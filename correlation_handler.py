"""
Module containing the class used to compute the correlation function
from same-event and mixed-event distributions.
"""

from ROOT import TH1F
from ROOT import TH2F
import logging
import numpy as np


class CorrelationHandler:
    """
    Class used to compute the correlation function
    from same-event and mixed-event k* distributions.
    Performs rebinning, normalisation and advanced operations.

    Parameters:
    ____________________________________________________

    inputSE may be a same-event k* (TH1F) or k* vs mult (TH2F) distribution
    inputME may be a same-event k* (TH1F) or k* vs mult (TH2F) distribution

    __se_k = same-event k* distribution (TH1F)
    __me_k = mixed-event k* distribution (TH1F)

    reweight method generates:
    __se_k = projects k* (TH1F) from k* vs mult (TH2F)
    __me_k = projects k* from k* vs mult and reweighs

    __se_mult = projects mult (TH1F) from k* vs mult (TH2F)
    __me_mult = projects mult from k* vs mult and reweighs

    __me_k_unw = unweighted k* projection
    __me_mult_unw = unweighted mult projection

    """

    def __init__(self, name, inputSE = None, inputME = None):
        # attributes
        self.__name = name
        self.__class = None

        self.__inputSE = None
        self.__inputME = None

        self.__se_k = None
        self.__me_k = None
        self.__cf = None

        if inputSE:
            if inputSE.ClassName() == 'TH2F':
                self.__class = 'TH2F'
                self.set_se(inputSE)
                self.__se_mult = None
                self.__me_mult = None

                self.__me_k_unw = None
                self.__me_mult_unw = None
                self.__cf_unw = None
            else:
                self.__class = 'TH1F'
                self.set_se(inputSE)
        else:
            print("input SE not set")
        if inputME:
            self.set_me(inputME)
        else:
            print("input ME not set")

    def set_se(self, inputSE = None):
        if self.__class == 'TH2F':
            self.__inputSE = inputSE.Clone(f'hSame_mult_{self.__name}')
        else:
            self.__inputSE = self.__se_k = inputSE.Clone(f'hSame_{self.__name}')

    def set_me(self, inputME = None):
        if self.__class == 'TH2F':
            self.__inputME = inputME.Clone(f'hMixed_mult_{self.__name}')
        else:
            self.__inputME = self.__me_k = inputME.Clone(f'hMixed_{self.__name}')

    def self_normalise(self):
        if self.__se_k:
            n_same = self.__se_k.GetNbinsX()
            integral_same = self.__se_k.Integral(1, n_same+1)
            self.__se_k.Scale(1/integral_same)
        else:
            logging.info('Same-event not set.')
        if self.__me_k:
            n_mixed = self.__me_k.GetNbinsX()
            integral_mixed = self.__me_k.Integral(1, n_mixed+1)
            self.__me_k.Scale(1/integral_mixed)
        else:
            logging.info('Mixed-event not set.')

    def move_to_MeV(self):
        if self.__se_k is None or self.__me_k is None:
            logging.error('Set same-event end/or mixed-event distribution first')
        else:
            if self.__cf:
                logging.info('After the change of units evaluate the CF again.')
            nbins = self.__se_k.GetNbinsX()
            low_edge = self.__se_k.GetBinLowEdge(1)
            up_edge = self.__se_k.GetBinLowEdge(nbins+1)
            self.__se_k.SetName('hSame_old')
            self.__me_k.SetName('hMixed_old')
            hSame_new = TH1F(
                f'hSame_{self.__name}', r';k* (MeV/#it{c}); Entries', nbins, low_edge*1000, up_edge*1000)
            hMixed_new = TH1F(
                f'hMixed_{self.__name}', r';k* (MeV/#it{c}); Entries', nbins, low_edge*1000, up_edge*1000)
            for i in range(0, nbins+2):
                hSame_new.SetBinContent(i, self.__se_k.GetBinContent(i))
                hSame_new.SetBinError(i, self.__se_k.GetBinError(i))
                hMixed_new.SetBinContent(i, self.__me_k.GetBinContent(i))
                hMixed_new.SetBinError(i, self.__me_k.GetBinError(i))
            self.__se_k = hSame_new
            self.__me_k = hMixed_new

    def move_to_GeV(self):
        if self.__se_k is None or self.__me_k is None:
            logging.error(
                'Set same-event end/or mixed-event distribution first')
        else:
            if self.__cf:
                logging.info(
                    'After the change of units evaluate the CF again.')
            nbins = self.__se_k.GetNbinsX()
            low_edge = self.__se_k.GetBinLowEdge(1)
            up_edge = self.__se_k.GetBinLowEdge(nbins+1)
            self.__se_k.SetName('hSame_old')
            self.__se_k.SetName('hMixed_old')
            hSame_new = TH1F(
                f'hSame_{self.__name}', r';k* (GeV/#it{c}); Entries', nbins, low_edge/1000, up_edge/1000)
            hMixed_new = TH1F(
                f'hMixed_{self.__name}', r';k* (GeV/#it{c}); Entries', nbins, low_edge/1000, up_edge/1000)
            for i in range(0, nbins+2):
                hSame_new.SetBinContent(i, self.__se_k.GetBinContent(i))
                hSame_new.SetBinError(i, self.__se_k.GetBinError(i))
                hMixed_new.SetBinContent(i, self.__me_k.GetBinContent(i))
                hMixed_new.SetBinError(i, self.__me_k.GetBinError(i))
            self.__se_k = hSame_new
            self.__me_k = hMixed_new

    def rebin(self, n_rebin=2) -> None:
        self.__se_k.Rebin(n_rebin)
        self.__me_k.Rebin(n_rebin)

    def make_correlation_function(self):
        self.__cf = self.__se_k.Clone(f'hCF_{self.__name}')
        self.__cf.Reset()
        self.__cf.GetYaxis().SetTitle('C(k*)')
        for i in range(1, self.__cf.GetNbinsX()):
            vSame = self.__se_k.GetBinContent(i)
            eSame = self.__se_k.GetBinError(i)
            vMixed = self.__me_k.GetBinContent(i)
            eMixed = self.__me_k.GetBinError(i)
            if vSame < 1.e-12 or vMixed < 1.e-12:
                continue
            else:
                vCF = vSame/vMixed
                eCF = np.sqrt((eMixed/vMixed)*(eMixed/vMixed) +
                              (eSame/vSame)*(eSame/vSame)) * vCF
                self.__cf.SetBinContent(i, vCF)
                self.__cf.SetBinError(i, eCF)

    def make_cf_unw(self):
        if not self.__me_k_unw:
            print("No unweighted mixed event distribution, run reweight() first!")
            return
        self.__cf_unw = self.__se_k.Clone(f'hCF_unw_{self.__name}')
        self.__cf_unw.Reset()
        self.__cf_unw.GetYaxis().SetTitle('C(k*)')
        for i in range(1, self.__cf_unw.GetNbinsX()):
            vSame1 = self.__se_k.GetBinContent(i)
            eSame1 = self.__se_k.GetBinError(i)
            vMixed1 = self.__me_k_unw.GetBinContent(i)
            eMixed1 = self.__me_k_unw.GetBinError(i)
            if vSame1 < 1.e-12 or vMixed1 < 1.e-12:
                continue
            else:
                vCF = vSame1 / vMixed1
                eCF = np.sqrt((eMixed1 / vMixed1)*(eMixed1 / vMixed1) +
                              (eSame1 / vSame1)*(eSame1 / vSame1)) * vCF
                self.__cf_unw.SetBinContent(i, vCF)
                self.__cf_unw.SetBinError(i, eCF)

    def normalise(self, low=0.6, high=1.) -> None:
        if self.__cf:
            low_bin = self.__cf.FindBin(low)
            low_edge = self.__cf.GetBinLowEdge(low_bin)
            high_bin = self.__cf.FindBin(high)
            high_edge = self.__cf.GetBinLowEdge(high_bin+1)
            cf_integral = self.__cf.Integral(low_bin, high_bin, 'width')
            self.__cf.Scale((high_edge-low_edge)/cf_integral)
        else:
            logging.error('Compute the correlation function first')

    def normalize_unw(self, low=0.6, high=1.):
        if self.__cf_unw:
            low_bin = self.__cf_unw.FindBin(low)
            low_edge = self.__cf_unw.GetBinLowEdge(low_bin)
            high_bin = self.__cf_unw.FindBin(high)
            high_edge = self.__cf_unw.GetBinLowEdge(high_bin+1)
            cf_integral = self.__cf_unw.Integral(low_bin, high_bin, 'width')
            self.__cf_unw.Scale((high_edge-low_edge)/cf_integral)
        else:
            logging.error('Compute the unweighted correlation function first')

    def reweight(self):
        if not self.__class == 'TH2F':
            print("Error in reweight(): Input is not a TH2F k* vs multiplicity distribution")
            return
        se_mult_integral = self.__inputSE.Integral()
        me_mult_integral = self.__inputME.Integral()

        se_k = self.__inputSE.ProjectionX("se_k", 1, self.__inputSE.GetNbinsY())
        se_mult = self.__inputSE.ProjectionY("se_mult", 1, self.__inputSE.GetNbinsX())

        me_k = self.__inputME.ProjectionX("me_k", 1, self.__inputME.GetNbinsY())
        me_k_unw = me_k.Clone("me_k_unw")

        me_mult = self.__inputME.ProjectionY("me_mult", 1, self.__inputME.GetNbinsX())
        me_mult_unw = me_mult.Clone("me_mult_unw")

        print(self.__inputSE.GetNbinsY())
        for ibin in range(1, self.__inputSE.GetNbinsY()):
            se_bin = self.__inputSE.ProjectionX("se_bin", ibin, ibin)
            se_bin.Sumw2()
            me_bin = self.__inputME.ProjectionX("me_bin", ibin, ibin)
            me_bin.Sumw2()
            se_ratio = se_bin.Integral() / se_mult_integral
            me_ratio = me_bin.Integral() / me_mult_integral

            print("Currently at bin ", +ibin)
            print(se_ratio)
            print(me_ratio)

            if (me_ratio > 0. and se_ratio > 0.):
                me_bin.Scale(se_ratio / me_ratio)
                me_mult.SetBinContent(ibin,me_bin.Integral())
                me_k.Add(me_bin,1)
            else:
                print("Attention: Multiplicity bin "+ibin+" has no entries. Please check!")

        self.__se_k = se_k
        self.__se_mult = se_mult
        self.__me_k = me_k
        self.__me_mult = me_mult
        self.__me_k_unw = me_k_unw
        self.__me_mult_unw = me_mult_unw


    ###   Getter functions   ###
    # Histogramms for same event distribution
    def get_se(self):
        return self.__inputSE

    def get_se_k(self):
        return self.__se_k

    def get_se_mult(self):
        return self.__se_mult

    # Histograms for the mixed event distribution
    def get_me(self):
        return self.__inputME

    def get_me_k(self):
        return self.__me_k

    def get_me_mult(self):
        return self.__me_mult

    def get_me_k_unw(self):
        return self.__me_k_unw

    def get_me_mult_unw(self):
        return self.__me_mult_unw

    # Correlation function
    def get_correlation_function(self):
        return self.__cf

    def get_cf_unw(self):
        return self.__cf_unw
