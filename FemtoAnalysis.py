import ROOT
import FemtoDreamReader as FDR
import CorrelationHandler as CH
from math import ceil

def UFFA_pp(path, fname, fdir, new_file, anal_type, hist_type, mc = False, bins = False, rebin = False):
    conf = config(path, fname, fdir, new_file, anal_type, hist_type, mc, bins, rebin)
    fdr = FDR.FemtoDreamReader(conf[0]+conf[1], conf[2])

    # block for the correlation calculation
    ch = cf_handler(fdr, conf)
    histos = ch.get_histos()
    #######################################

    # file saver
    fds = FemtoDreamSaver(fname, histos, conf)
    ############

# class that handles the retrieving of histos and computing of correlation functions
class cf_handler():
    def __init__(self, FileReader, conf):
        self._file = FileReader
        self._atype = conf[3]   # analysis type
        self._htype = conf[4]   # histo type
        self._mc    = conf[5]   # bool monte carlo data
        self._bins  = conf[6]   # bin range for differential
        self._rebin = conf[7]   # rebin factors for all se, me, cf plots
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
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin)
            if self._mc:
                histos_mc, histos_unw_mc = getIntegrated(self._se_mc, self._me_mc, self._htype, self._rebin)
        elif self._atype == 2:      # differential analysis
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin)
            if self._mc:
                histos_mc = getDifferential(self._se_mc, self._me_mc, self._htype, self._bins, self._rebin)

        return [histos, histos_unw, histos_mc, histos_unw_mc, self._event, self._tracks, self._tracks_mc]

# class that handles the saving of the histograms in the correct file structure
class FemtoDreamSaver():
    last_edit = None

    def __init__(self, ofile, histos, conf):
        self._ofile = ofile     # output file name
        self._histos = histos   # all histos from the cf_handler
        self._opath = conf[0]   # path to file
        self._oname = conf[1]   # name of file
        self._odir  = conf[2]   # TDirectory name
        self._atype = conf[3]   # analysis type
        self._htype = conf[4]   # histo type
        self._mc    = conf[5]   # bool monte carlo data
        self._bins  = conf[6]   # bin range for differential
        self._rebin = conf[7]   # rebin factors for all se, me, cf plots
        self._newfile = conf[8] # "new", "recreate" or "update"
        self._save_histos()

    # write all histos from the list 'histos'
    def _write(self, histos):
        for hist in histos:
            hist.Write()

    # make directory in dir_root with name dir_name and save all histos from the list 'histos'
    def _mkdir_write(self, dir_root, dir_name, histos):
        dir_new = dir_root.mkdir(dir_name)
        dir_new.cd()
        try:
            self._write(histos)
        except:
            pass
        dir_root.cd()
        del dir_new

    # function that looks if 'file' exists and if yes it append -n, where n depends on if another file was already created
    def _file_exists(self, file):
        if not ROOT.gSystem.AccessPathName(file):
            name, ext = file.rsplit('.')
            digit = 1
            while not ROOT.gSystem.AccessPathName(name + '-' + str(digit) + '.' + ext):
                digit += 1
            file = name + '-' + str(digit) + '.' + ext
        return file

    # function that saves all the histos in the correct file format
    # Obacht! very chunky function!
    def _save_histos(self):
        # for "new" file: rename output if the output file already exists
        if self._newfile == "new":
            new_name = self._file_exists(self._opath + "UFFA_" + self._ofile)
            print(new_name + " created!")
            self._ofile = new_name
            ofile_name = self._opath + self._ofile
            FemtoDreamSaver.last_edit = self._ofile
        else:
            ofile_name = self._opath + "UFFA_" + self._ofile

        if FemtoDreamSaver.last_edit:
            ofile_name = FemtoDreamSaver.last_edit
        ofile = ROOT.TFile(ofile_name, self._newfile)

        if self._atype == 1:        # integrated
            hist_std = self._histos[0]      # [[iSE, iME], [se, me, cf], [se rebin1, me rebin1, cf rebin1], ...]
            hist_unw = self._histos[1]      # [[me unw, cf unw], [me unw rebin1, cf unw rebin1], ...]
            hist_in  = hist_std[0]          # [iSE, iME]
            hist_smc = hist_std[1]          # [se, me, cf]
            if self._mc:            # monte carlo
                hist_std_mc = self._histos[2]
                hist_in_mc  = hist_std_mc[0]
                hist_smc_mc = hist_std_mc[1]
                if self._htype == 2:
                    hist_unw_mc = self._histos[3]
        elif self._atype == 2:      # differential
            hist_std = self._histos[0]      # [[iSE, iME], [se, me, cf], [[se rebin1, me rebin1, cf rebin1], ...], ...]
            hist_in  = hist_std[0]          # [iSE, iME]
            if self._mc:            # monte carlo
                hist_std_mc = self._histos[2]
                hist_in_mc  = hist_std_mc[0]

        hist_event = self._histos[4]
        hist_track = self._histos[5]
        if self._mc:
            hist_track_mc = self._histos[6]
            hist_pur = getPurity(hist_track[0], hist_track_mc[0])

        if self._odir == "" or self._odir == "femto-dream-pair-task-track-track":
            dir_root = ofile.mkdir("_std")
        else:
            dir_root = ofile.mkdir(self._odir)
        dir_root.cd()

        if self._atype == 1:            # integrated
            self._write(hist_in)        # iSE, iME
            self._write(hist_smc)       # se, me, cf
            if self._htype == 2:
                self._write(hist_unw[0])   # me unw, cf unw
            if self._rebin:             # all rebinned se, me, cf, me unw, cf unw, in individual directories
                for n in range(len(self._rebin)):
                    dir_rebin = dir_root.mkdir("rebin: " + str(self._rebin[n]))
                    dir_rebin.cd()
                    self._write(hist_std[2][n])        # [[iSE, iME], [no rebin], [1st rebin], [2nd rebin], ...]
                    if self._htype == 2:
                        self._write(hist_unw[1][n])    # [[no rebin], [1st rebin], [2nd rebin], ...]
                    dir_root.cd()
                    del dir_rebin
            self._mkdir_write(dir_root, "Event", hist_event)
            self._mkdir_write(dir_root, "Tracks_one", hist_track)
            if self._mc:
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self._write(hist_in_mc)
                self._write(hist_smc_mc)
                if self._htype == 2:
                    self._write(hist_unw_mc[0])
                if self._rebin:         # rebin directories in main directory
                    for n in range(len(self._rebin)):
                        dir_rebin = dir_mc.mkdir("rebin: " + str(self._rebin[n]))
                        dir_rebin.cd()
                        self._write(hist_std_mc[2][n])
                        if self._htype == 2:
                            self._write(hist_unw_mc[1][n])
                        dir_mc.cd()
                        del dir_rebin
                self._mkdir_write(dir_mc, "Tracks_oneMC", hist_track_mc)
                dir_root.cd()
                hist_pur.Write()

        elif self._atype == 2:          # differential
            self._write(hist_in)
            for n in range(len(self._bins) - 1):    # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                dir_bin = dir_root.mkdir("bin: " + str(n + 1))
                dir_bin.cd()
                if self._rebin:         # rebin directories inside each mt/mult bin directory
                    self._write(hist_std[1 + 2*n])
                    for m in range(len(self._rebin)):
                        dir_rebin = dir_bin.mkdir("rebin: " + str(self._rebin[m]))
                        dir_rebin.cd()
                        self._write(hist_std[2 + 2*n][m])
                        dir_bin.cd()
                        del dir_rebin
                else:
                    self._write(hist_std[1 + n])
                dir_root.cd()
                del dir_bin
            self._mkdir_write(dir_root, "Event", hist_event)
            self._mkdir_write(dir_root, "Tracks_one", hist_track)
            if self._mc:                # monte carlo directory
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self._write(hist_in_mc)
                for n in range(len(self._bins) - 1):    # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                    dir_bin = dir_mc.mkdir("bin: " + str(n + 1))
                    dir_bin.cd()
                    if self._rebin:     # rebin directories inside each mt/mult bin directory
                        self._write(hist_std_mc[1 + 2*n])
                        for m in range(len(self._rebin)):
                            dir_rebin = dir_bin.mkdir("rebin: " + str(self._rebin[m]))
                            dir_rebin.cd()
                            self._write(hist_std[2 + 2*n][m])
                            dir_bin.cd()
                            del dir_rebin
                    else:
                        self._write(hist_std[1 + n])
                    dir_mc.cd()
                    del dir_bin
                dir_root.cd()
                del dir_mc

# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int:
        rebin = [rebin]
    elif type(rebin) != list and type(rebin) != bool:
        print("rebin factor is not an int or a list!")
        exit()
    if rebin:
        rebin_list.extend(rebin)
    return rebin_list

# returns list with the configured settings
# [path, file name, infile dir, analysis type, hist type, mc, bin range, rebin factors]
def config(path, fname, fdir, new_file, anal_type, hist_type, mc, bins, rebin):
    k_keys = ['k', 'kstar', '1']
    mult_keys = ['mult', 'kmult', '2']
    mt_keys = ['mt', 'kmt', '3']
    int_keys = ['int', 'integrated', '1']
    dif_keys = ['diff', 'differential', '2']

    ipath = path            # path to file
    iname = fname           # file name
    idir = fdir             # TDirectory in file

    if new_file not in [1, 2, 3]:
        print("\nInput error:\tWrong 'new_file' option!\nOptions:\n\t1 -> \"new\"\n\t2 -> \"recreate\"\n\t3 -> \"update\"\n")
        exit()

    if new_file == 1:
        new_file = "new"
    elif new_file == 2:
        new_file = "recreate"
    elif new_file == 3:
        new_file = "update"
    else:
        print("Input error: This ain't it chief!")

    # analysis type
    anal_type = str(anal_type).lower()
    if anal_type in int_keys:
        anal_type = 1
    elif anal_type in dif_keys:
        anal_type = 2
    else:
        print("\nInput error:\tWrong analysis type!\nOptions:\n\tintegrated: \t'int', 'integrated', 1 \n\tdifferential: \t'dif', 'differential', 2\n")
        exit()

    # histogram type
    hist_type = str(hist_type).lower()
    if hist_type in k_keys:
        hist_type = 1
    elif hist_type in mult_keys:
        hist_type = 2
    elif hist_type in mt_keys:
        hist_type = 3
    else:
        print("\nInput error:\tWrong histo type!\nOptions:\n\tTH1 k*: \t'k', 'kstar', 1 \n\tTH2 k-mult: \t'mult', 'kmult', 2 \n\tTH2 k-mt: \t'mt', 'kmt', 3\n")
        exit()
    mc = bool(mc)
    rebin = bin2list(rebin)

    return [ipath, iname, idir, anal_type, hist_type, mc, bins, rebin, new_file]

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
        name = "%.2f-%.2f" % (yAxis.GetBinCenter(limits[n - 1]), yAxis.GetBinCenter(limits[n]))
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
def getCorrelation(se, me, name, conf):
    ch = CH.CorrelationHandler(name, se, me)
    #ch.normalize()
    ch.make_cf()
    ch.normalize_cf(0.24, 0.34)
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
def getDifferential(se, me, hist_type, bins, rebin):
    histos = []
    if hist_type == 2:
        conf = "mult: "
        histos.append([se.Clone("SE kmult"), me.Clone("ME kmult")])
    elif hist_type == 3:
        conf = "mt: "
        histos.append([se.Clone("SE kmT"), me.Clone("ME kmT")])
    else:
        print("getDifferential: no kmT or kmult input!")
        exit()

    mt_histos = getBinRangeHistos(se, me, bins)
    for name, se, me in mt_histos:
        histos.append(getCorrelation(se, me, name, conf + name))
        if rebin:       # append a list of rebinned [se, me, cf] after the original [se, me, cf]
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                histos_rebin.append(getCorrelation(se_rebin, me_rebin, name, conf + name + rebin_conf))
            histos.append(histos_rebin)
    return histos

# returns a list of [[iSE, iME], [se, me, cf]] for rel pair k* input
# or reweights and returns ([[iSE, iME], [se, me, cf]], [me_unw, cf_unw]) for kmult
# and [[iSE, iME], [se, me, cf]] for kmT
def getIntegrated(se, me, hist_type, rebin):
    histos = []
    histos_unw = []
    if hist_type == 1:      # k* input
        histos.append([se.Clone("SE kstar"), me.Clone("ME kstar")])
    elif hist_type == 2:    # kmult input
        histos.append([se.Clone("SE kmult"), me.Clone("ME kmult")])
        hReweight = reweight(se, me)
        se = hReweight[0]
        me = hReweight[1]
        me_unw = hReweight[2]
    elif hist_type == 3:    # kmT input
        histos.append([se.Clone("SE kmT"), me.Clone("ME kmT")])
        hReweight = reweight(se, me)
        se = hReweight[0]
        me = hReweight[2]   # unweighted me, i.e. normal me projection of the kmT histo

    histos.append(getCorrelation(se, me, "cf", ""))
    if rebin:               # append rebinned histos to list of histos
        histos_rebin = []
        for factor in rebin:
            se_rebin = rebin_hist(se, factor)
            me_rebin = rebin_hist(me, factor)
            rebin_conf = " rebin: " + str(factor)
            histos_rebin.append(getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf))
        histos.append(histos_rebin)
    if hist_type == 2:      # 2nd list with unweighted histos
        se, me, cf = getCorrelation(se, me_unw, "cf_unw", "unweighted")
        histos_unw.append([me.Clone("ME unw"), cf.Clone("CF unw")])
        if rebin:           # append rebinned histos to list of histos
            histos_rebin = []
            for factor in rebin:
                se_rebin = rebin_hist(se, factor)
                me_rebin = rebin_hist(me, factor)
                rebin_conf = " rebin: " + str(factor)
                se_rebin, me_rebin, cf_rebin = getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf)
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
