import ROOT
import FileUtils as FU
import FemtoDreamSaver as FDS
import FemtoDreamReader as FDR
import CorrelationHandler as CH
import TemplateFit as TF

def UFFA(settings):
    match settings['function']:
        case 'cf':
            UFFA_cf(settings)
        case 'tf':
            UFFA_tf(settings)

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
    bins = [conf['bins'], conf['rebin']] if conf['rebin'] else conf['bins']
    TF.TemplateFit(ofile, dca_data, dca_mcplots, conf['tftype'], conf['namelist'], conf['fitrange'], bins, conf['outDir'])

def UFFA_syst(fname, fdir, new_file, atype, htype, mc = None, bins = None, rebin = None, dirOut = None):
    conf = config(dirOut, fname, fdir, new_file, atype, htype, mc, bins, rebin)
    fdr = FDR.FemtoDreamReader(fname, fdir)

    # default cf
    ch = cf_handler(fdr, conf)
    cf, cf_unw = ch.get_cf()
    cf, cf_rebin = cf

    # systematics
    syst = Systematics(cf)

    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th2
    file_dir = fdr.get_dir()
    while (fdr.cd(file_dir + "_var" + n)):
        ch_var = cf_handler(fdr, conf)
        cf_var, cf_var_rebin = ch_var.get_cf()[0]

        # add variation
        syst.AddVar(cf_var)

        n += 1
        del ch_var

    # generate th2 plots for systematics
    syst.GenSyst()

    # save output
    ofile = ROOT.TFile("UFFA_syst_" + conf[1], "recreate")
    root = ofile.mkdir(conf[2])
    root.cd()

    cf.Write()
    h2cf_var.Write()
    h2cf_dif.Write()
    cf_sys.Write()
    cf_dev.Write()

    for i in range(len(rebin)):
        new_dir = root.mkdir("rebin: " + rebin[i])
        new_dir.cd()

        cf_rebin[i].Write()
        h2cf_var_rebin[i].Write()
        h2cf_dif_rebin[i].Write()
        cf_sys_rebin[i].Write()
        cf_dev_rebin[i].Write()

# class that returns the systematics of a cf
# add variations with AddVar(var) before calling GenSyst()
class Systematics():
    ybins = 512
    def __init__(self, cf):
        self._cf = cf
        self._xaxis = cf.GetXaxis()
        self._xbins = cf.GetNbinsX()
        self._var = ROOT.TH2D("CF", "CF", xbins, xaxis.GetXmin(), xaxis.GetXmax(), Systematics.ybins, 0, 3)
        self._dif = ROOT.TH2D("diff", "diff", xbins, xaxis.GetXmin(), xaxis.GetXmax(), Systematics.ybins*2, -3, 3)
        self._sys = ROOT.TH1D("syst", "syst", xbins, xaxis.GetXmin(), xaxis.GetXmax())
        self._dev = ROOT.TH1D("dev", "dev", xbins, xaxis.GetXmin(), xaxis.GetXmax())

    def AddVar(self, cf_var):
        for i in range(cf_var.GetEntries()):
            self._var.Fill(i + 1, cf_var.GetBinContent(i))

    def GenSyst(self):
        for i in range(self.xbins):
            cf_proj = self._var.ProjectionY("cf xbin" + str(i + 1), i + 1, i + 1)
            dev = cf_proj.GetStdDev()
            self._dev.SetBinContent(i + 1, dev)
            cont_def = self._cf.GetBinContent(i)
            for j in range(ybins):
                cont_var = cf_proj.GetBinContent(j)
                self._dif.Fill(i + 1, cont_def - cont_var)                      # fill th2 with difference to default
            dif_proj = self._dif.ProjectionY("diff xbin" + str(i + 1), i + 1, i + 1)
            proj_min = dif_proj.GetBinContent(dif_proj.FindFirstBinAbove(0))    # minimum value of difference
            proj_max = dif_proj.GetBinContent(dif_proj.FindLastBinAbove(0))     # maximum value of difference
            self._sys.SetBinContent(i + 1, (proj_max - proj_min) / (12**0.5))   # assume a square distribution

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
        self._atype = conf['atype']         # analysis type
        self._htype = conf['htype']         # histo type
        self._mc    = conf['mc']            # bool monte carlo data
        self._bins  = conf['bins']          # bin range for differential
        self._rebin = conf['rebin']         # rebin factors for all se, me, cf plots
        self._norm  = conf['normalize']     # normalization range
        self._se = None
        self._me = None
        self._se_mc = None
        self._me_mc = None
        self._event = None
        self._tracks = None
        self._tracks_mc = None
        self._get_histos()

    # retrieves histos from the provided file reader
    def _get_histos(self):
        if self._htype == 1:        # TH1 kstar
            self._se, self._me = self._file.get_kstar()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kstar_mc()
        elif self._htype == 2:      # TH2 k-mult
            self._se, self._me = self._file.get_kmult()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmult_mc()
        elif self._htype == 3:      # TH2 k-mt
            self._se, self._me = self._file.get_kmt()
            if self._mc:
                self._se_mc, self._me_mc = self._file.get_kmt_mc()
        try:
            self._event = self._file.get_event()        # event histos
            self._tracks = self._file.get_tracks()      # track histos
            if self._mc:
                self._tracks_mc = self._file.get_tracks_mc()
        except:
            pass

    # computes the cf for integrated or differential analysis and for mc data
    # and returns the histos for all the different options
    # [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]
    def get_histos(self):
        histos = []
        histos_mc = []
        histos_unw = []
        histos_unw_mc = []
        if self._atype == 1:        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
            if self._mc:
                histos_mc, histos_unw_mc = getIntegrated(self._se_mc, self._me_mc, self._htype, self._rebin, self._norm)
        elif self._atype == 2:      # differential analysis
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)
            if self._mc:
                histos_mc = getDifferential(self._se_mc, self._me_mc, self._htype, self._bins, self._rebin, self._norm)

        return [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]

    # returns a list of cf and their rebinned version
    # [cf, [rebin 1, rebin 2, ...]], same for unweighted if integrated analysis
    def get_cf(self):
        histos = []
        histos_unw = []
        cf_list = []
        cf_list_unw = []
        if self._atype == 1:        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin, self._norm)
            cf_list_unw.append(histos_unw[0][1])
            cf_list_unw.append([])
        elif self._atype == 2:      # differential analysis
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin, self._norm)

        cf_list.append(histos[1][2])
        cf_list.append([])
        if self._rebin:
            for n in range(len(self._rebin)):
                cf = histos[2][n][2]
                cf.SetName("CF rebin: " + self._rebin[n])
                cf_list[1].append(cf)
                if self._atype == 1:
                    cf_unw = histos_unw[1][n][1]
                    cf_unw.SetName("CF unw rebin: " + self._rebin[n])
                    cf_list_unw[1].append(cf_unw)

        return [cf_list, cf_list_unw]

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
            "rebin":        None,
            "atype":        None,
            "htype":        None,
            "tftype":       None,
            "data":         None,
            "templates":    None,
            "namelist":     None,
            "fitrange":     None,
            "normalize":    None,
            "debug":        False
        }

    k_keys    = ['k', 'kstar', '1']
    mult_keys = ['mult', 'kmult', '2']
    mt_keys   = ['mt', 'kmt', '3']
    int_keys  = ['int', 'integrated', '1']
    dif_keys  = ['diff', 'differential', '2']

    # function to be used
    if 'function' in dic_conf:
        dic['function'] = dic_conf['function']

    # type of particle pair
    if 'pair' in dic_conf:
        dic['pair'] = dic_conf['pair']

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
            dic['atype'] = 1
        elif atype in dif_keys:
            dic['atype'] = 2

    # histogram type
    if 'htype' in dic_conf:
        htype = dic_conf['htype']
        if type(htype) == str:
            htype = htype.lower()
        if htype in k_keys:
            dic['htype'] = 1
        elif htype in mult_keys:
            dic['htype'] = 2
        elif htype in mt_keys:
            dic['htype'] = 3

    # template fit type
    if 'tftype' in dic_conf:
        tftype = dic_conf['tftype'].lower()
        if tftype == 'dca':
            dic['tftype'] = 1
        elif tftype == 'cpa':
            dic['tftype'] = 2

    # mc data file
    if 'mc' in dic_conf:
        dic['mc'] = dic_conf['mc']

    # mc data file TDir/TList
    if 'mcTDir' in dic_conf:
        dic['mcTDir'] = dic_conf['mcTDir']

    # bin ranges
    if 'bins' in dic_conf:
        dic['bins'] = dic_conf['bins']

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
    if 'debug' in dic_conf:
        dic['debug'] = bool(dic_conf['debug'])

    return dic

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

# calculate purity
def getPurity(hPt, hPt_mc):
    ratio = hPt.Clone("hPt_ratio")
    ratio.Divide(hPt_mc)
    return ratio

# helper function for the CF
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
# and a list of rebinned [se, me, cf] appended to the firt list for rebin factors
# [[iSE, iME], [bin 1], [bin 1 rebin], [bin 2], [bin 2 rebin]...]
def getDifferential(se, me, htype, bins, rebin, norm):
    histos = []
    if htype == 2:
        conf = "mult: "
        histos.append([se.Clone("SE kmult"), me.Clone("ME kmult")])
    elif htype == 3:
        conf = "mt: "
        histos.append([se.Clone("SE kmT"), me.Clone("ME kmT")])
    else:
        print("getDifferential: no kmT or kmult input!")
        exit()

    mt_histos = getBinRangeHistos(se, me, bins)
    for name, se, me in mt_histos:
        histos.append(getCorrelation(se, me, name, conf + name, norm))
        if rebin:       # append a list of rebinned [se, me, cf] after the original [se, me, cf]
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                histos_rebin.append(getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf, norm))
            histos.append(histos_rebin)
    return histos

# returns a list of [[iSE, iME], [se, me, cf]] for rel pair k* input
# or reweights and returns ([[iSE, iME], [se, me, cf]], [[me_unw, cf_unw]]) for kmult
# and [[iSE, iME], [se, me, cf]] for kmT
def getIntegrated(se, me, htype, rebin, norm):
    histos = []
    histos_unw = []
    if htype == 1:      # k* input
        histos.append([se.Clone("SE kstar"), me.Clone("ME kstar")])
    elif htype == 2:    # kmult input
        histos.append([se.Clone("SE kmult"), me.Clone("ME kmult")])
        hReweight = reweight(se, me)
        se          = hReweight[0]
        me          = hReweight[1]
        me_unw      = hReweight[2]
        se_mult     = hReweight[3]
        me_mult     = hReweight[4]
        me_mult_unw = hReweight[5]
        histos[0].append(se_mult.Clone("SE mult"))
        histos[0].append(me_mult.Clone("ME mult"))
        histos[0].append(me_mult_unw.Clone("ME mult unw"))
    elif htype == 3:    # kmT input
        histos.append([se.Clone("SE kmT"), me.Clone("ME kmT")])
        hReweight = reweight(se, me)
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
        histos.append(histos_rebin)
    if htype == 2:      # 2nd list with unweighted histos
        se, me, cf = getCorrelation(se, me_unw, "cf_unw", "unweighted", norm)
        histos_unw.append([me.Clone("ME unw"), cf.Clone("CF unw")])
        if rebin:           # append rebinned histos to list of histos
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                se_rebin, me_rebin, cf_rebin = getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf, norm)
                histos_rebin.append([me_rebin.Clone("ME unw"), cf_rebin.Clone("CF unw")])
            histos_unw.append(histos_rebin)

    return histos, histos_unw

# projects and reweights se and me from kmult histos
# returns [se, me, me unweighted, se mult, me mult, me mult unweighted]
def reweight(iSE, iME):
    se_int = iSE.Integral()
    me_int = iME.Integral()

    se_k = iSE.ProjectionX("se_k", 1, iSE.GetNbinsY())
    me_k = iSE.ProjectionX("me_k", 1, iME.GetNbinsY())
    me_k.Reset()

    se_mult = iSE.ProjectionY("se_mult", 1, iSE.GetNbinsX())
    me_mult = iME.ProjectionY("me_mult", 1, iME.GetNbinsX())
    me_mult.Reset()

    me_k_unw = iME.ProjectionX("me_k_unw", 1, iME.GetNbinsY())
    me_mult_unw = iME.ProjectionY("me_mult_unw", 1, iME.GetNbinsX())

    for n in range(1, iSE.GetNbinsY()):
        se_n = iSE.ProjectionX("se_bin", n, n)
        me_n = iME.ProjectionX("me_bin", n, n)
        #se_n.Sumw2()
        #me_n.Sumw2()
        se_ratio = se_n.Integral() / se_int
        me_ratio = me_n.Integral() / me_int

        if me_ratio > 0. and se_ratio > 0.:
            me_n.Scale(se_ratio / me_ratio)
            me_mult.SetBinContent(n, me_n.Integral())
            me_k.Add(me_n, 1)

    return [se_k, me_k, me_k_unw, se_mult, me_mult, me_mult_unw]

# returns rebinned copy of histo
def rebin_hist(input_histo, binning):
    histo = input_histo.Clone()
    histo = histo.Rebin(binning)
    return histo

# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int:
        rebin = [rebin]
    elif type(rebin) != list:
        return None
    rebin_list.extend(rebin)
    return rebin_list
