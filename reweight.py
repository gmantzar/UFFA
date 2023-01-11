from ROOT import TH1F
from ROOT import TH2F
import numpy as np

class Reweighter:
    """
    Class used to reweight the mixed event distribution
    based on the multiplicity distribution in the
    same event distribution

    Parameters:
    _____________________________________________________

    same_event_mult = k* vs. multiplicity for same event (TH2F)
    mixed_event_mult = k* vs. multiplicity for mixed event distribution (TH2F)

    Can calculate:

    same_event = k* distribution of the same event (Projection of same_event_mult)
    mixed_event_unw = mixed-event k* distribution (Projection of mixed_event_mult)
    mixed_event = weighed mixed-event k* distribution
    same_mult = multiplicity distribution of the same-event (Projection)
    mixed_mult_unw = unweighted multiplicity distribution of the mixed-events (Projection)
    mixed_mult = weighted multiplicity distribution of the mixed-events (Projection)

    --> The weighting reweights the mixed-event distribution such that its multiplicity
        matches that of the same-event distribution
    """

    def __init__(self,name, same_event_mult=None, mixed_event_mult=None) -> None:
        self.name = name
        if same_event_mult:
            self.set_same_event_mult(same_event_mult)
        else:
            self.same_event_mult = None
        if mixed_event_mult:
            self.set_mixed_event_mult(mixed_event_mult)
        else:
            self.mixed_event_mult = None
        self.same_event = None
        self.same_mult = None
        self.mixed_event = None
        self.mixed_event_unw = None
        self.mixed_mult = None
        self.mixed_mult_unw = None
        self.testhisto = None

    def get_testhisto(self):
        return self.testhisto

    #Histogramms for same event distribution
    def set_same_event_mult(self, same_event_mult=None) -> None:
        self.same_event_mult = same_event_mult.Clone(f'hSame_mult_{self.name}')

    def get_same_event_mult(self):
        return self.same_event_mult

    def get_same_event(self):
        return self.same_event

    def get_same_mult(self):
        return self.same_mult

    #histograms for the mixed event distribution
    def set_mixed_event_mult(self, mixed_event_mult=None) -> None:
        self.mixed_event_mult = mixed_event_mult.Clone(f'hMixed_mult_{self.name}')

    def get_mixed_event_mult(self):
        return self.mixed_event_mult

    def get_mixed_event(self):
        return self.mixed_event

    def get_mixed_event_unw(self):
        return self.mixed_event_unw

    def get_mixed_mult(self):
        return self.mixed_mult

    def get_mixed_mult_unw(self):
        return self.mixed_mult_unw

    def do_reweight(self):
        se_mult_integral = self.same_event_mult.Integral()
        me_mult_integral = self.mixed_event_mult.Integral()

        se = self.same_event_mult.ProjectionX("se",1,self.same_event_mult.GetNbinsY())
        #se.Reset("ICESM") #not necessary. We can just use the projection for the same event distribution
        se_mult = self.same_event_mult.ProjectionY("se_mult",1,self.same_event_mult.GetNbinsX())

        me_unw = self.mixed_event_mult.ProjectionX("me_unw",1,self.mixed_event_mult.GetNbinsY())
        me_mult_unw = self.mixed_event_mult.ProjectionY("me_mult_unw",1,self.mixed_event_mult.GetNbinsX())
        me = me_unw.Clone("me")
        me.Reset("ICESM")
        me_mult = me_mult_unw.Clone("me_mult")
        me_mult.Reset("ICESM")

        print(self.same_event_mult.GetNbinsY())
        for ibin in range(1,self.same_event_mult.GetNbinsY()):
        #for ibin in range(1,5):
            SE_i = self.same_event_mult.ProjectionX("SE_i",ibin,ibin)
            SE_i.Sumw2()
            ME_i = self.mixed_event_mult.ProjectionX("ME_i",ibin,ibin)
            ME_i.Sumw2()
            SEratio = SE_i.Integral()/se_mult_integral
            MEratio = ME_i.Integral()/me_mult_integral
            print("Currently at bin ",+ibin)
            print(SEratio)
            print(MEratio)

            if(MEratio>0. and SEratio>0.):
                ME_i.Scale(SEratio/MEratio)
                me_mult.SetBinContent(ibin,ME_i.Integral())
                me.Add(ME_i,1)
                #se.Add(SE_i,1) not necessary. The same event distribution does not change
            else:
                print("Attention: Multiplicity bin "+ibin+" has no entries. Please check")

        self.same_event=se
        self.same_mult=se_mult
        self.mixed_event_unw=me_unw
        self.mixed_mult_unw=me_mult_unw
        self.mixed_event=me
        self.mixed_mult=me_mult
