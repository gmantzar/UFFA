import ROOT
import FileUtils as FU
import FemtoDreamSaver as FDS
import FemtoDreamReader as FDR
import CorrelationHandler as CH
import TemplateFit as TF

def UFFA(settings):
    conf = config(settings)
    match conf['function']:
        case 'cf':
            UFFA_cf(conf)
        case 'tf':
            UFFA_tf(conf)
        case 'syst':
            if conf['htype'] == 'mtmult':
                UFFA_syst_3d(conf)
            else:
                UFFA_syst(conf)

# correlation function
def UFFA_cf(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
    ch = cf_handler(fdr, conf)
    fds = FDS.FemtoDreamSaver(conf, ch.get_histos())

# template fits
def UFFA_tf(settings):
    conf = config(settings)
    if conf['file']:
        fdr1 = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])
        dca_data = fdr1.get_dca()
    elif conf['data']:
        dca_data = conf['data']
    else:
        print('UFFA_tf: Missing input data!')
    if conf['templates']:
        if type(conf['templates']) == str:
            fdr2 = FDR.FemtoDreamReader(conf['templates'], conf['mcTDir'])
            dca_mcplots = fdr2.get_dca_mc()
        else:
            dca_mcplots = conf['templates']
    else:
        dca_mcplots = fdr1.get_dca_mc()

    fds = FDS.FemtoDreamSaver(settings)
    ofile = fds.getFile()
    if 'rebin' in conf and conf['rebin']:
        bins = [conf['bins'], conf['rebin']]
    else:
        bins = conf['bins']
    TF.TemplateFit(ofile, dca_data, dca_mcplots, conf['tftype'], conf['namelist'], conf['fitrange'], bins, conf['outDir'])

# systematics
def UFFA_syst(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

    # default cf
    ch = cf_handler(fdr, conf)
    cf, cf_unw = ch.get_cf()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

    cf_list = []
    if conf['rebin']:
        len_rebin = len(conf['rebin'])

    if conf['atype'] == 'int':                              # integrated
        ck, ck_rebin = cf[0]
        cf_list.append([ck, ck_rebin])
        syst = [[Systematics(ck), []]]                          # [[syst cf, [rebins]]]
        if conf['rebin']:
            for i in range(len_rebin):
                syst[0][1].append(Systematics(ck_rebin[i]))
    elif conf['atype'] == 'dif':                            # differential
        syst = []
        for n, [ck, ck_rebin] in enumerate(cf):
            cf_list.append([ck, ck_rebin])
            syst.append([Systematics(ck), []])                  # [[syst cf, [rebins]], [bin2...], ...]
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[n][1].append(Systematics(ck_rebin[i]))

    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th1
    file_dir = fdr.get_dir();
    niter = 1           # start at _Var01
    fdr.cd(0)           # class method of FileSaver to return to root of file
    while (fdr.cd(file_dir + f"_Var{niter:02d}")):
        # allows to include/exclude specific variations
        if conf['exclude'] and file_dir + f"_Var{niter:02d}" in conf['exclude']:
            niter += 1
            continue
        elif conf['include']:
            if fdr.get_dir() in conf['include']:
                pass
            else:
                niter += 1
                continue

        ch_var = cf_handler(fdr, conf)
        cf_var, cf_var_unw = ch_var.get_cf()
        for n, [ck_var, ck_var_rebin] in enumerate(cf_var):
            syst[n][0].AddVar(ck_var)
            if conf['rebin']:
                for i in range(len_rebin):
                    syst[n][1][i].addVar(ck_var_rebin[i])
        niter += 1
        del ch_var

    # generate th2 plots for systematics
    for n in range(len(syst)):
        syst[n][0].GenSyst()
        if conf['rebin']:
            for i in range(len_rebin):
                syst[n][1][i].GenSyst()

    syst_plots = []                                         # [[[cf, diff, syst, dev], [rebins]], [bin2...], ...]
    for n in range(len(syst)):
        syst_plots.append([syst[n][0].GetAll(), []])
        if conf['rebin']:
            for i in range(len_rebin):
                syst_plots[n][1].append(syst[n][1][i].GetAll())

    # generates the graphs with the systematic errors for the cf and the rebinned entries
    tgraphs = []
    for n, (hist, hist_rebin) in enumerate(cf_list):
        tgraphs.append([ROOT.TGraphErrors(), []])
        for i in range(1, hist.GetNbinsX() + 1):
            tgraphs[n][0].SetName("CF syst graph")
            tgraphs[n][0].SetPoint(i - 1, hist.GetBinCenter(i), hist.GetBinContent(i))
            tgraphs[n][0].SetPointError(i - 1, 0, syst_plots[n][0][2].GetBinContent(i))
        if conf['rebin']:
            for i in range(len_rebin):
                tgraphs[n][1].append(ROOT.TGraphErrors())
                for j in range(1, hist.GetNbinsX() + 1):
                    tgraphs[n][1][i].SetName("CF syst graph")
                    tgraphs[n][1][i].SetPoint(j - 1, hist_rebin[i].GetBinCenter(j), hist_rebin[i].GetBinContent(j))
                    tgraphs[n][1][i].SetPointError(j - 1, 0, syst_plots[n][1][j][2].GetBinContent(j))

    histos = (cf_list, syst_plots, tgraphs)
    fds = FDS.FemtoDreamSaver(conf, histos)

# systematics
def UFFA_syst_3d(settings):
    conf = config(settings)
    fdr = FDR.FemtoDreamReader(conf['fullpath'], conf['fileTDir'])

    # default cf
    ch = cf_handler(fdr, conf)
    histos = ch.get_cf_3d()                                # [[cf, [rebins]], [bin2...], ...], [[cf unw, [rebins]], [bin2...], ...]

    syst = []
    syst_plots = []

    if conf['rebin']:
        len_rebin = len(conf['rebin'])

    # [[[bin1-1 cf, [rebin cf]], [bin1-2 cf, [rebin cf]], ...], [[bin2-1 cf, [rebin cf]], [bin2-2 cf, [rebin cf]], ...], ...]
    for n, bin1 in enumerate(histos):
        syst.append([])
        for nn, [cf, cf_rebin] in enumerate(bin1):
            syst[n].append([Systematics(cf), []])
            if conf['rebin']:
                for nnn in range(len_rebin):
                    syst[n][nn][1].append(Systematics(cf_rebin[nnn]))

    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th1
    file_dir = fdr.get_dir();
    niter = 1           # start at _Var01
    fdr.cd(0)           # class method of FileSaver to return to root of file
    while (fdr.cd(file_dir + f"_Var{niter:02d}")):
        # allows to include/exclude specific variations
        if conf['exclude'] and file_dir + f"_Var{niter:02d}" in conf['exclude']:
            niter += 1
            continue
        elif conf['include']:
            if fdr.get_dir() in conf['include']:
                pass
            else:
                niter += 1
                continue

        ch_var = cf_handler(fdr, conf)
        histos_var = ch_var.get_cf_3d()

        for n, bin1 in enumerate(histos_var):
            for nn, [cf, cf_rebin] in enumerate(bin1):
                syst[n][nn][0].AddVar(cf)
                if conf['rebin']:
                    for nnn in range(len_rebin):
                        syst[n][nn][1][nnn].AddVar(cf_rebin[nnn])
        niter += 1
        del ch_var

    # generate th2 plots for systematics
    for n, bin1 in enumerate(syst):
        syst_plots.append([])
        for nn, bin2 in enumerate(bin1):
            syst[n][nn][0].GenSyst()
            syst_plots[n].append([syst[n][nn][0].GetAll(), []])
            if conf['rebin']:
                for nnn in range(len_rebin):
                    syst[n][nn][1][nnn].GenSyst()
                    syst_plots[n][nn][1].append(syst[n][nn][1][nnn].GetAll())

    # generates the graphs with the systematic errors for the cf and the rebinned entries
    tgraphs = []
    for n, bin1 in enumerate(histos):
        tgraphs.append([])
        for nn, [hist, hist_rebin] in enumerate(bin1):
            tgraphs[n].append([ROOT.TGraphErrors(), []])
            for nnn in range(1, hist.GetNbinsX() + 1):
                tgraphs[n][nn][0].SetName("CF syst graph")
                tgraphs[n][nn][0].SetPoint(nnn - 1, hist.GetBinCenter(nnn), hist.GetBinContent(nnn))
                tgraphs[n][nn][0].SetPointError(nnn - 1, 0, syst_plots[n][nn][0][2].GetBinContent(nnn))
            if conf['rebin']:
                for nnn in range(len_rebin):
                    tgraphs[n][nn][1].append(ROOT.TGraphErrors())
                    for nnnn in range(1, hist.GetNbinsX() + 1):
                        tgraphs[n][nn][1][nnn].SetName("CF syst graph")
                        tgraphs[n][nn][1][nnn].SetPoint(nnnn - 1, hist_rebin[nnn].GetBinCenter(nnnn), hist_rebin[nnn].GetBinContent(nnnn))
                        tgraphs[n][nn][1][nnn].SetPointError(nnnn - 1, 0, syst_plots[n][nn][1][nnn][2].GetBinContent(nnnn))

    all_histos = (histos, syst_plots, tgraphs)
    fds = FDS.FemtoDreamSaver(conf, all_histos)

# class that returns the systematics of a cf
# add variations with AddVar(var) before calling GenSyst()
# GetAll() returns [th2 cf, th2 difference, th1 systematics, th1 std dev]
class Systematics():
    counter = 0
    #ybins = 128
    ybins = 256
    #ybins = 512
    #ybins = 1024
    #ybins = 2048
    def __init__(self, cf):
        self._cf = cf
        self._xaxis = cf.GetXaxis()
        self._xbins = cf.GetNbinsX()
        if Systematics.counter == 0:
            text = ""
        else:
            text = " " + str(Systematics.counter)
        self._var = ROOT.TH2D("CF th2" + text, "CF th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins, 0, 10)
        self._dif = ROOT.TH2D("diff th2" + text, "diff th2", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax(), Systematics.ybins*2, -10, 10)
        self._sys = ROOT.TH1D("syst th1" + text, "syst th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        self._dev = ROOT.TH1D("dev th1" + text, "dev th1", self._xbins, self._xaxis.GetXmin(), self._xaxis.GetXmax())
        Systematics.counter = Systematics.counter + 1

    def AddVar(self, cf_var):
        for i in range(1, self._xbins + 1):
            self._var.SetBinContent(i, self._var.GetYaxis().FindBin(cf_var.GetBinContent(i)), 1)

    def GenSyst(self):
        for i in range(1, self._xbins + 1):
            cf_proj = self._var.ProjectionY("cf xbin" + str(i), i, i)
            dev = cf_proj.GetStdDev()
            self._dev.SetBinContent(i, dev)
            cont_def = self._cf.GetBinContent(i)
            for j in range(Systematics.ybins):
                var_min = cf_proj.GetBinCenter(cf_proj.FindFirstBinAbove(0))
                var_max = cf_proj.GetBinCenter(cf_proj.FindLastBinAbove(0))
                #self._dif.Fill(cf_proj.GetBinCenter(i), var_min)                # fill th2 with difference to default
                #self._dif.Fill(cf_proj.GetBinCenter(i), var_max)                # fill th2 with difference to default
                self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_min), 1)
                self._dif.SetBinContent(i, self._dif.GetYaxis().FindBin(var_max), 1)
            dif_proj = self._dif.ProjectionY("diff xbin" + str(i), i, i)
            proj_min = dif_proj.GetBinCenter(dif_proj.FindFirstBinAbove(0))     # minimum value of difference
            proj_max = dif_proj.GetBinCenter(dif_proj.FindLastBinAbove(0))      # maximum value of difference
            self._sys.SetBinContent(i, (proj_max - proj_min) / (12**0.5))       # assume a square distribution
        self._var.SetDirectory(0)
        self._dif.SetDirectory(0)
        self._sys.SetDirectory(0)
        self._dev.SetDirectory(0)

    def SetBinning(self, n):
        Systematics.ybins = n

    def GetVar(self):
        return self._var

    def GetDiff(self):
        return self._dif

    def GetSyst(self):
        return self._sys

    def GetDev(self):
        return self._dev

    def GetAll(self):
        return [self._var, self._dif, self._sys, self._dev]

# class that handles the retrieving of histos and computing of correlation functions
class cf_handler():
    def __init__(self, FileReader, conf):
        self._file  = FileReader
        self._pair  = conf['pair']
        self._atype = conf['atype']             # analysis type
        self._htype = conf['htype']             # histo type
        self._mc    = conf['mc']                # bool monte carlo data
        self._bins  = conf['bins']              # bin range for differential
        self._diff3d = conf['diff3d']           # which axis to split first in a 3D analysis
        self._binsdiff3d = conf['binsdiff3d']   # bin range for the first differential split in case of a 3D analysis
        self._rebin = conf['rebin']             # rebin factors for all se, me, cf plots
        self._norm  = conf['normalize']         # normalization range
        self._se = None
        self._me = None
        self._se_mc = None
        self._me_mc = None
        self._event = None
        self._tracks = None
        self._tracks_mc = None
        self._v0 = None
        self._get_histos()

    # retrieves histos from the provided file reader
    def _get_histos(self):
        if self._htype == 'k':        # TH1 kstar
            self._se, self._me = self._file.get_kstar()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kstar_mc()
        elif self._htype == 'mult':      # TH2 k-mult
            self._se, self._me = self._file.get_kmult()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmult_mc()
        elif self._htype == 'mt':      # TH2 k-mt
            self._se, self._me = self._file.get_kmt()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmt_mc()
        elif self._htype == 'mtmult':
            self._se, self._me = self._file.get_kmtmult()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmtmult_mc()
        self._event = self._file.get_event()
        self._tracks = self._file.get_tracks()
        if self._mc:
            self._tracks_mc = self._file.get_tracks_mc()
        if self._pair == 'pl':
            self._v0 = self._file.get_v0()

    # computes the cf for integrated or differential analysis and for mc data
    # and returns the histos for all the different options
    # [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]
    def get_histos(self):
        histos = []
        histos_mc = []
        histos_unw = []
        histos_unw_mc = []
        if self._atype == 'int':        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
            if self._mc:
                histos_mc, histos_unw_mc = getIntegrated(self._se_mc, self._me_mc, self._htype, self._rebin, self._norm)
        elif self._atype == 'dif':      # differential analysis
            if self._htype == 'mtmult': # 3D differantial analysis
                histos = getDifferential3D(self._se, self._me, self._diff3d, self._binsdiff3d, self._htype, self._bins, self._rebin, self._norm)
            else:
                histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)
                if self._mc:
                    histos_mc = getDifferential(self._se_mc, self._me_mc, self._htype, self._bins, self._rebin, self._norm)

        return [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc, self._v0]

    # returns a list of cf and their rebinned version
    # [[cf, [rebin 1, rebin 2, ...]], [bin2...], ...] same for unweighted if integrated analysis
    def get_cf(self):
        histos = []
        histos_unw = []
        cf_list = []
        cf_list_unw = []
        if self._atype == 'int':        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
            if self._htype == 'mult':
                cf_list_unw.append(histos_unw[1])
                cf_list_unw.append([])
        elif self._atype == 'dif':      # differential analysis
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)

        cf_list.append([histos[1][2], []])                          # cf, for differential 1st bin
        if self._rebin:
            for n in range(len(self._rebin)):
                cf_list[0][1].append(histos[1][3][n][2])            # rebinned cf
                if self._atype == 'int' and self._htype == 'mult':
                    cf_list_unw[1].append(histos_unw[2][n][1])      # rebinned unw cf for integrated
        if self._atype == 'dif':
            for n in range(1, len(self._bins) - 1):
                cf_list.append([histos[n + 1][2], []])
                if self._rebin:
                    for m in range(len(self._rebin)):
                        cf_list[n + 1][1].append(histos[n + 1][3][n][2])    # rebinned cf appended to rebin list
        return [cf_list, cf_list_unw]

# returns all the cf's for a 3D mt/mult histo
# [[[bin1-1 cf, [rebin cf]], [bin1-2 cf, [rebin cf]], ...], [[bin2-1 cf, [rebin cf]], [bin2-2 cf, [rebin cf]], ...], ...]
    def get_cf_3d(self):
        cf_list = []

        histos = getDifferential3D(self._se, self._me, self._diff3d, self._binsdiff3d, self._htype, self._bins, self._rebin, self._norm)
        for n in range(1, len(histos)):
            cf_list.append([])
            for nn in range(1, len(histos[n])):
                cf_list[n - 1].append([histos[n][nn][2], []])

        return cf_list

# splits th2 in section based on provided bins
def getBinRangeHistos(iSE, iME, bins):
    xAxis = iSE.GetXaxis()
    yAxis = iSE.GetYaxis()

    if type(bins) == list:
        limits = []
        for mt_value in bins:
            bin_value = yAxis.FindBin(float(mt_value))
            limits.append(bin_value)
    else:
        print("Error in getBinRangeHistos: bin input \"" + str(bins) + "\" not a list of ranges!")
        exit()

    histos = []
    for n in range(1, len(limits)):
        name = "[%.2f-%.2f)" % (yAxis.GetBinLowEdge(limits[n - 1]), yAxis.GetBinLowEdge(limits[n]))
        se = iSE.ProjectionX("se_k", limits[n - 1], limits[n])
        me = iME.ProjectionX("me_k", limits[n - 1], limits[n])
        histos.append([name, se.Clone(), me.Clone()])

    return histos

# splits th3 in section based on provided bins
def getBinRangeHistos3D(iSE, iME, diff1, binsdiff1):
    if diff1 == 'mt':
        diffAxisSE = iSE.GetYaxis()
        diffAxisME = iME.GetYaxis()
        projOpt = "zx"
    elif diff1 == 'mult':
        diffAxisSE = iSE.GetZaxis()
        diffAxisME = iME.GetZaxis()
        projOpt = "yx"
    else:
        print("Error in getBinRangeHistos: diff1 axis not known. Please choose either 'mT' or 'mult'")
        exit()

    if type(binsdiff1) == list:
        limits = []
        for diff_value in binsdiff1:
            bin_value = diffAxisSE.FindBin(float(diff_value))
            limits.append(bin_value)
    else:
        print("Error in getBinRangeHistos: bin input \"" + str(binsdiff1) + "\" not a list of ranges!")
        exit()

    histos = []
    for n in range(1, len(limits)):
        name = "[%.2f-%.2f)" % (diffAxisSE.GetBinLowEdge(limits[n - 1]), diffAxisSE.GetBinLowEdge(limits[n]))
        diffAxisSE.SetRange(limits[n - 1], limits[n])
        diffAxisME.SetRange(limits[n - 1], limits[n])
        se = iSE.Project3D("SE_"+projOpt+"_"+name)
        me = iME.Project3D("ME_"+projOpt+"_"+name)
        histos.append([se.Clone(), me.Clone()])

    return histos

# calculate purity
def getPurity(hPt, hPt_mc):
    ratio = hPt.Clone("hPt_ratio")
    ratio.Divide(hPt_mc)
    return ratio

# helper function for the CF
# output: [se, me, cf]
def getCorrelation(se, me, name, conf, norm = None):
    ch = CH.CorrelationHandler(name, se, me)
    #ch.normalize()
    ch.make_cf()
    minmax = norm if norm else [0.24, 0.34]
    ch.normalize_cf(minmax[0], minmax[1])
    se = ch.get_se().Clone("SE")
    me = ch.get_me().Clone("ME")
    cf = ch.get_cf().Clone("CF")
    se.SetTitle(conf)
    me.SetTitle(conf)
    cf.SetTitle(conf)
    del ch
    return [se, me, cf]

# returns [[iSE, iME], [se, me, cf]] for a list of mt or mult ranges
# [[iSE, iME], [se, me, cf, [rebin: [...], [...], ...], [bin 2 [rebin]], ...]
def getDifferential(iSE, iME, htype, bins, rebin, norm):
    histos = []
    if htype == 'mult':
        conf = "mult: "
        histos.append([iSE.Clone("SE kmult"), iME.Clone("ME kmult")])
    elif htype == 'mt':
        conf = "mt: "
        histos.append([iSE.Clone("SE kmT"), iME.Clone("ME kmT")])
    else:
        print("getDifferential: no kmT or kmult input!")
        exit()

    mt_histos = getBinRangeHistos(iSE, iME, bins)
    for n, [name, se, me] in enumerate(mt_histos, 1):
        histos.append(getCorrelation(se, me, name, conf + name, norm))
        histos[n].append([])
        if rebin:       # append a list of rebinned [se, me, cf] in the original [se, me, cf]
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                histos[n][3].append(getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf, norm))
    return histos

# [[iSE, iME], [[1st proj SE, 1st proj ME], [se, me, cf, [rebin], [bin 2 [rebin]]]], ...]
def getDifferential3D(iSE, iME, diff3d, binsdiff3d, htype, bins, rebin, norm):
    histos = []
    histos.append([iSE.Clone("SE kmTmult"), iME.Clone("ME kmTmult")])

    diff3d_histos = getBinRangeHistos3D(iSE, iME, diff3d, binsdiff3d)

    htypeSplit2 = ""
    if diff3d == 'kmult':
        htypeSplit2 = "kmt"
    elif diff3d == 'mt':
        htypeSplit2 = "mult"
    for se, me in diff3d_histos:
        histos.append(getDifferential(se, me, htypeSplit2, bins, rebin, norm))

    return histos

# returns a list of [[iSE, iME], [se, me, cf]] for rel pair k* input
# or reweights and returns ([[iSE, iME], [se, me, cf, [rebin]]], [me_unw, cf_unw, [rebin]]) for kmult
# and [[iSE, iME], [se, me, cf, [rebin]]] for kmT
def getIntegrated(iSE, iME, htype, rebin, norm):
    histos = []
    histos_unw = []
    if htype == 'k':      # k* input
        histos.append([iSE.Clone("SE kstar"), iME.Clone("ME kstar")])
        se = iSE
        me = iME
    elif htype == 'mult':    # kmult input
        histos.append([iSE.Clone("SE kmult"), iME.Clone("ME kmult")])
        hReweight = reweight(iSE, iME)
        se          = hReweight[0]
        me          = hReweight[1]
        me_unw      = hReweight[2]
        se_mult     = hReweight[3]
        me_mult     = hReweight[4]
        me_mult_unw = hReweight[5]
        histos[0].append(se_mult.Clone("SE mult"))
        histos[0].append(me_mult.Clone("ME mult"))
        histos[0].append(me_mult_unw.Clone("ME mult unw"))
    elif htype == 'mt':    # kmT input
        histos.append([iSE.Clone("SE kmT"), iME.Clone("ME kmT")])
        hReweight = reweight(iSE, iME)
        se = hReweight[0]
        me = hReweight[2]   # unweighted me, i.e. normal me projection of the kmT histo

    histos.append(getCorrelation(se, me, "cf", "", norm))
    if rebin:               # append rebinned histos to list of histos
        histos_rebin = []
        for factor in rebin:
            se_rebin = rebin_hist(se, factor)
            me_rebin = rebin_hist(me, factor)
            rebin_conf = " rebin: " + str(factor)
            histos_rebin.append(getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm))
        histos[1].append(histos_rebin)
    if htype == 'mult':      # 2nd list with unweighted histos
        se, me, cf = getCorrelation(se, me_unw, "cf_unw", "unweighted", norm)
        histos_unw.append(me.Clone("ME unw"))
        histos_unw.append(cf.Clone("CF unw"))
        histos_unw.append([])
        if rebin:           # append rebinned histos to list of histos
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                se_rebin, me_rebin, cf_rebin = getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm)
                histos_unw[2].append([me_rebin.Clone("ME unw"), cf_rebin.Clone("CF unw")])

    return histos, histos_unw

# projects and reweights se and me from kmult histos
# returns [se, me, me unweighted, se mult, me mult, me mult unweighted]
def reweight(iSE, iME):
    se_k = iSE.ProjectionX("se_k")
    me_k = iSE.ProjectionX("me_k")
    se_int = se_k.Integral()
    me_int = me_k.Integral()

    se_mult = iSE.ProjectionY("se_mult")
    me_mult = iME.ProjectionY("me_mult")

    me_k_unw = iME.ProjectionX("me_k_unw")
    me_mult_unw = iME.ProjectionY("me_mult_unw")

    me_k.Reset("ICESM")
    me_mult.Reset("ICESM")

    for n in range(iSE.GetNbinsY()):
        se_n = iSE.ProjectionX("se_bin", n, n)
        me_n = iME.ProjectionX("me_bin", n, n)

        se_ratio = se_n.Integral() / se_int
        me_ratio = me_n.Integral() / me_int

        if me_ratio > 0. and se_ratio > 0.:
            me_n.Scale(se_ratio / me_ratio)
            me_mult.SetBinContent(n, me_n.Integral())
            me_k.Add(me_n)

    return [se_k, me_k, me_k_unw, se_mult, me_mult, me_mult_unw]

# returns rebinned copy of histo
def rebin_hist(input_histo, binning):
    histo = input_histo.Clone()
    histo = histo.Rebin(binning)
    return histo

# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int or type(rebin) == str:
        rebin = [rebin]
    elif type(rebin) != list:
        return None
    rebin_list.extend(rebin)
    return rebin_list

# returns list with the configured settings
# [path, file name, infile dir, analysis type, hist type, mc, bin range, rebin factors]
def config(dic_conf):
    dic = {
            "function":     None,
            "pair":         None,
            "path":         "",
            "file":         None,
            "fullpath":     None,
            "fileTDir":     "",
            "newfile":      None,
            "mc":           None,
            "mcTDir":       "",
            "outDir":       "",
            "rename":       None,
            "bins":         None,
            "diff3d":       "",
            "binsdiff3d":   None,
            "rebin":        None,
            "atype":        None,
            "htype":        None,
            "tftype":       None,
            "data":         None,
            "templates":    None,
            "namelist":     None,
            "fitrange":     None,
            "normalize":    None,
            "include":      None,
            "exclude":      None,
            "debug":        False
        }

    k_keys      = ['k', 'kstar', '1']
    mult_keys   = ['mult', 'kmult', '2']
    mt_keys     = ['mt', 'kmt', '3']
    mtmult_keys = ['mtmult','kmtmult', '4']
    int_keys = ['int', 'integrated', '1']
    dif_keys = ['diff', 'dif', 'differential', '2']

    # function to be used
    if 'function' in dic_conf:
        dic['function'] = dic_conf['function']

    # type of particle pair
    if 'pair' in dic_conf:
        if dic_conf['pair']:
            dic['pair'] = dic_conf['pair'].lower()

    # input directory
    if 'path' in dic_conf:
        dic['path'] = dic_conf['path']

    # file name, file directory
    if 'file' in dic_conf:
        path_name = dic_conf['file'].rsplit('/', 1)
        if len(path_name) == 1:
            dic['file']  = path_name[0]
        else:
            dic['path'] = FU.path_expand(path_name[0]) + '/'
            dic['file']  = path_name[1]
    dic['fullpath'] = dic['path'] + dic['file']

    # file TDir/TList
    if 'fileTDir' in dic_conf:
        if dic_conf['fileTDir'] == "":
            dic['fileTDir'] = "femto-dream-pair-task-track-track"
        elif dic_conf['fileTDir'][0] == '_':
            dic['fileTDir'] = "femto-dream-pair-task-track-track" + dic_conf['fileTDir']
        else:
            dic['fileTDir'] = dic_conf['fileTDir']

    if 'outDir' in dic_conf:
        if dic_conf['outDir'] != "" and dic_conf['outDir']:
            dic['outDir'] = FU.path_expand(dic_conf['outDir'])
            if ROOT.gSystem.AccessPathName(dic['outDir']):
                print("output directory \"" + dic['outDir'] + "\" does not exist!")
                exit()
    else:
        dic['outDir'] = dic['path']

    # rename output file
    if 'rename' in dic_conf:
        dic['rename'] = dic_conf['rename']

    # create file
    if 'newfile' in dic_conf:
        if dic_conf['newfile'] in [1, "new"]:
            dic['newfile'] = "new"
        if dic_conf['newfile'] in [2, "recreate"]:
            dic['newfile'] = "recreate"
        if dic_conf['newfile'] in [3, "update"]:
            dic['newfile'] = "update"

    # analysis type
    if 'atype' in dic_conf:
        atype = dic_conf['atype']
        if type(atype) == str:
            atype = atype.lower()
        if atype in int_keys:
            dic['atype'] = 'int'
        elif atype in dif_keys:
            dic['atype'] = 'dif'

    # histogram type
    if 'htype' in dic_conf:
        htype = dic_conf['htype']
        if type(htype) == str:
            htype = htype.lower()
        if htype in k_keys:
            dic['htype'] = 'k'
        elif htype in mult_keys:
            dic['htype'] = 'mult'
        elif htype in mt_keys:
            dic['htype'] = 'mt'
        elif htype in mtmult_keys:
            dic['htype'] = 'mtmult'

    # template fit type
    if 'tftype' in dic_conf:
        if dic_conf['tftype']:
            tftype = dic_conf['tftype'].lower()
            if tftype == 'dca':
                dic['tftype'] = 'dca'
            elif tftype == 'cpa':
                dic['tftype'] = 'cpa'

    # mc data file
    if 'mc' in dic_conf:
        dic['mc'] = dic_conf['mc']

    # mc data file TDir/TList
    if 'mcTDir' in dic_conf:
        dic['mcTDir'] = dic_conf['mcTDir']

    # bin ranges
    if 'bins' in dic_conf:
        dic['bins'] = dic_conf['bins']

   # which axis to be used for the first split in a 3D analysis
    if 'diff3d' in dic_conf:
        diff3d = dic_conf['diff3d']
        if type(diff3d) == str:
            diff3d = diff3d.lower()
        if diff3d in mult_keys:
            dic['diff3d'] = 'mult'
        elif diff3d in mt_keys:
            dic['diff3d'] = 'mt'

    # binning of the first differential split in case of a 3D analysis;
    # if used, 'bins' will be used for the binning of the second split
    if 'binsdiff3d' in dic_conf:
        dic['binsdiff3d'] = dic_conf['binsdiff3d']

    # rebin factor/s
    if 'rebin' in dic_conf:
        dic['rebin'] = bin2list(dic_conf['rebin'])

    # general purpose data entry
    if 'data' in dic_conf:
        dic['data'] = dic_conf['data']

    # template file/plots
    if 'templates' in dic_conf:
        dic['templates'] = dic_conf['templates']

    # names for template file/plots
    if 'namelist' in dic_conf:
        dic['namelist'] = dic_conf['namelist']

    # fitrange
    if 'fitrange' in dic_conf:
        dic['fitrange'] = dic_conf['fitrange']

    # normalize cf in range
    if 'normalize' in dic_conf:
        dic['normalize'] = dic_conf['normalize']

    # normalize cf in range
    if 'include' in dic_conf:
        dic['include'] = bin2list(dic_conf['include'])

    # normalize cf in range
    if 'exclude' in dic_conf:
        dic['exclude'] = bin2list(dic_conf['exclude'])

    # normalize cf in range
    if 'debug' in dic_conf:
        dic['debug'] = bool(dic_conf['debug'])

    return dic

# debug class
class cDEBUG():
    def __init__(self):
        self.n = 0
        self.x = 2
        self.y = 1
        self.N = 2
        self.h = 1024
        self.w = 768
        self.i = 0
        self.canvas = None

    def setSize(self, h, w):
        self.h = h
        self.w = w

    def divide(self, x, y):
        self.x = x
        self.y = y
        self.N = x * y

    def makeCanvas(self):
        self.canvas = ROOT.TCanvas("c", "c")
        self.canvas.SetCanvasSize(self.h, self.w)
        self.canvas.Divide(self.x, self.y)
        self.canvas.SetWindowSize(self.w + (self.w - self.canvas.GetWw()), self.h + (self.h - self.canvas.GetWh()));

    def add(self, histo, n = None):
        if n:
            self.canvas.cd(n)
            histo.Draw()
        elif self.n < self.N:
            self.n += 1
            self.canvas.cd(self.n)
            histo.Draw()
        else:
            return None

    def update(self):
        self.canvas.Update()

    def print(self):
        self.canvas.SaveAs("cDEBUG_" + str(self.i) + ".pdf")
        self.i += 1

