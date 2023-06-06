import ROOT
import FemtoDreamReader as FDR
import CorrelationHandler as CH
import v2TemplateFit as TF
import time

def UFFA_pp(fname, fdir, new_file, atype, htype, mc = None, bins = None, rebin = None, dirOut = None):
    conf = config(dirOut, fname, fdir, new_file, atype, htype, mc, bins, rebin)
    # file reader
    fdr = FDR.FemtoDreamReader(fname, fdir)
    # correlation calculation
    ch = cf_handler(fdr, conf)
    # file saver
    fds = FemtoDreamSaver(ch.get_histos(), conf)

def UFFA_pp_syst(fname, fdir, new_file, atype, htype, mc = None, bins = None, rebin = None, dirOut = None):
    conf = config(dirOut, fname, fdir, new_file, atype, htype, mc, bins, rebin)
    fdr = FDR.FemtoDreamReader(fname, fdir)

    se = fdr.get_se()
    xaxis = se.GetXaxis()
    xbins = se.GetNbinsX()
    ybins = 512

    # default cf
    ch = cf_handler(fdr, conf)
    cf, cf_rebin = ch.get_cf()

    # histograms
    h2cf_var = ROOT.TH2D("CF", "CF", xbins, xaxis.GetXmin(), xaxis.GetXmax(), ybins, 0, 3)
    h2cf_dif = ROOT.TH2D("CF diff", "CF diff", xbins, xaxis.GetXmin(), xaxis.GetXmax(), ybins*2, -3, 3)
    cf_sys = ROOT.TH1D("syst", "syst", xbins, xaxis.GetXmin(), xaxis.GetXmax())
    cf_dev = ROOT.TH1D("std dev", "std dev", xbins, xaxis.GetXmin(), xaxis.GetXmax())
    # lists for the rebinned versions
    h2cf_var_rebin = []
    h2cf_dif_rebin = []
    cf_sys_rebin = []
    cf_dev_rebin = []
    if rebin:
        for i in range(len(rebin)):
            h2cf_var_rebin.append(ROOT.TH2D("cf rebin: " + rebin[i], "cf rebin: " + rebin[i], xbins, xaxis.GetXmin(), xaxis.GetXmax(), ybins, 0, 3))
            h2cf_dif_rebin.append(ROOT.TH2D("cf diff rebin: " + rebin[i], "cf diff rebin: " + rebin[i], xbins, xaxis.GetXmin(), xaxis.GetXmax(), ybins*2, -3, 3))
            cf_sys_rebin.append(ROOT.TH1D("syst rebin: " + rebin[i], "syst rebin: " + rebin[i], xbins, xaxis.GetXmin(), xaxis.GetXmax()))
            cf_dev_rebin.append(ROOT.TH1D("std dev rebin: " + rebin[i], "std dev rebin: " + rebin[i], xbins, xaxis.GetXmin(), xaxis.GetXmax()))
            h2cf_var_rebin.RebinX(rebin[i])
            h2cf_dif_rebin.RebinX(rebin[i])
            cf_sys_rebin.Rebin(rebin[i])
            cf_dev_rebin.Rebin(rebin[i])

    n = 0
    file_dir = fdr.get_dir()
    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th2
    while (fdr.cd(file_dir + "_var" + n)):
        ch_var = cf_handler(fdr, conf)
        cf_var, cf_var_rebin = ch_var.get_cf()
        del ch_var

        # fill the th2 with each cf
        for i in range(cf_var.GetEntries()):
            h2cf_var.Fill(i + 1, cf_var.GetBinContent(i))
        if rebin:
            for i in range(len(rebin)):
                for j in range(cf_var_rebin[i]):
                    h2cf_var_rebin[i].Fill(j + 1, cf_var_rebin[i].GetBinContent(j))
        n += 1

    # loop over each x bin and get the content
    for i in range(xbins):
        value_def = cf.GetBinContent(i)
        cf_proj = h2cf_var.ProjectionY("cf bin: " + str(i + 1), i + 1, i + 1)
        dev = cf_proj.GetStdDev()
        cf_dev.SetBinContent(i + 1, dev)
        for j in range(ybins):
            value_var = cf_proj.GetBinContent(j)
            h2cf_dif.Fill(i + 1, value_def - value_var)                     # fill th2 with difference to default
        dif_proj = h2cf_dif.ProjectionY("diff bin: " + str(i + 1), i + 1, i + 1)
        proj_min = dif_proj.GetBinContent(dif_proj.FindFirstBinAbove(0))    # minimum value of difference
        proj_max = dif_proj.GetBinContent(dif_proj.FindLastBinAbove(0))     # maximum value of difference
        cf_sys.SetBinContent(i + 1, (proj_max - proj_min) / (12**0.5))      # assume a square distribution
        # same loop for all the rebinned versions
        if rebin:
            for j in range(len(rebin)):
                value_def = cf_rebin[j]
                cf_proj_rebin = h2cf_var_rebin[j].Projection("cf rebin: " + rebin[j] + " bin: " + str(i + 1), i + 1, i + 1)
                dev = cf_proj_rebin.GetStdDev()
                cf_dev_rebin[j].SetBinContent(i + 1, dev)
                for k in range(ybins):
                    value_var = cf_proj_rebin.GetBinContent(k)
                    h2cf_dif_rebin[j].Fill(i + 1, value_def - value_var)
                dif_proj_rebin = h2cf_dif_rebin[j].Projection("diff rebin: " + rebin[j] + " bin: " + str(i + 1), i + 1, i + 1)
                proj_min = dif_proj_rebin.GetBinContent(dif_proj_rebin.FindFirstBinAbove(0))
                proj_max = dif_proj_rebin.GetBinContent(dif_proj_rebin.FindLastBinAbove(0))
                cf_sys_rebin[j].SetBinContent(i + 1, (proj_max - proj_min) / (12**0.5))

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

def UFFA_pp_syst2(fname, fdir, new_file, atype, htype, mc = None, bins = None, rebin = None, dirOut = None):
    conf = config(dirOut, fname, fdir, new_file, atype, htype, mc, bins, rebin)
    fdr = FDR.FemtoDreamReader(fname, fdir)

    # default cf
    ch = cf_handler(fdr, conf)
    cf, cf_rebin = ch.get_cf()

    h2cf_var_rebin = []
    h2cf_dif_rebin = []
    cf_sys_rebin = []
    cf_dev_rebin = []

    n = 0
    file_dir = fdr.get_dir()
    # loop over data variations in file and calculate the cf for each
    # which is then saved in a th2 from which the systematic error is computed and saved in a th2
    while (fdr.cd(file_dir + "_var" + n)):
        ch_var = cf_handler(fdr, conf)
        cf_var, cf_var_rebin = ch_var.get_cf()
        del ch_var

        # fill the th2 with each cf
        for i in range(cf_var.GetEntries()):
            h2cf_var.Fill(i + 1, cf_var.GetBinContent(i))
        if rebin:
            for i in range(len(rebin)):
                for j in range(cf_var_rebin[i]):
                    h2cf_var_rebin[i].Fill(j + 1, cf_var_rebin[i].GetBinContent(j))
        n += 1

    # loop over each x bin and get the content
    for i in range(xbins):
        value_def = cf.GetBinContent(i)
        cf_proj = h2cf_var.ProjectionY("cf bin: " + str(i + 1), i + 1, i + 1)
        dev = cf_proj.GetStdDev()
        cf_dev.SetBinContent(i + 1, dev)
        for j in range(ybins):
            value_var = cf_proj.GetBinContent(j)
            h2cf_dif.Fill(i + 1, value_def - value_var)                     # fill th2 with difference to default
        dif_proj = h2cf_dif.ProjectionY("diff bin: " + str(i + 1), i + 1, i + 1)
        proj_min = dif_proj.GetBinContent(dif_proj.FindFirstBinAbove(0))    # minimum value of difference
        proj_max = dif_proj.GetBinContent(dif_proj.FindLastBinAbove(0))     # maximum value of difference
        cf_sys.SetBinContent(i + 1, (proj_max - proj_min) / (12**0.5))      # assume a square distribution
        # same loop for all the rebinned versions
        if rebin:
            for j in range(len(rebin)):
                value_def = cf_rebin[j]
                cf_proj_rebin = h2cf_var_rebin[j].Projection("cf rebin: " + rebin[j] + " bin: " + str(i + 1), i + 1, i + 1)
                dev = cf_proj_rebin.GetStdDev()
                cf_dev_rebin[j].SetBinContent(i + 1, dev)
                for k in range(ybins):
                    value_var = cf_proj_rebin.GetBinContent(k)
                    h2cf_dif_rebin[j].Fill(i + 1, value_def - value_var)
                dif_proj_rebin = h2cf_dif_rebin[j].Projection("diff rebin: " + rebin[j] + " bin: " + str(i + 1), i + 1, i + 1)
                proj_min = dif_proj_rebin.GetBinContent(dif_proj_rebin.FindFirstBinAbove(0))
                proj_max = dif_proj_rebin.GetBinContent(dif_proj_rebin.FindLastBinAbove(0))
                cf_sys_rebin[j].SetBinContent(i + 1, (proj_max - proj_min) / (12**0.5))

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
            cf_proj = self._var.ProjectionY("cf xbin: " + str(i + 1), i + 1, i + 1)
            dev = cf_proj.GetStdDev()
            self._dev.SetBinContent(i + 1, dev)
            cont_def = self._cf.GetBinContent(i)
            for j in range(ybins):
                cont_var = cf_proj.GetBinContent(j)
                self._dif.Fill(i + 1, cont_def - cont_var)                      # fill th2 with difference to default
            dif_proj = self._dif.ProjectionY("diff xbin: " + str(i + 1), i + 1, i + 1)
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

TF_conf = {
        "type": None,
        "fit_range": None,
        "pt_bins": None,
        "names": None,
        "data": None,
        "data_dir": None,
        "mc": None,
        "mc_dir": None,
        "rename": None,
        "dirOut": None
        }

def TemplateFit(dcacpa, fit_range, pt_bins, dca_mcplots_names, data_file, data_dir,
                mc_file = None, mc_dir = None, fname = None, dirOut = None):
    if type(data_file) == str:
        fdr1 = FDR.FemtoDreamReader(data_file, data_dir)
        dca_data = fdr1.get_dca()
        if not dca_data:
            dca_data = fdr1.get_histo("Tracks_one/DCAxy")
    else:
        dca_data = data_file
    if mc_file:
        if type(mc_file) == str:
            fdr2 = FDR.FemtoDreamReader(mc_file, mc_dir)
            dca_mcplots = fdr2.get_dca_mc()
        else:
            dca_mcplots = mc_file
    else:
        dca_mcplots = fdr1.get_dca_mc()

    if not fname:
        fname = data_file
    TF.TemplateFit(fname, dca_data, dca_mcplots, dcacpa, dca_mcplots_names, fit_range, pt_bins, dirOut)

def tf_binning(pt_bins):
    if type(pt_bins) == int:
        pass


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

    # returns a list of cf and their rebinned version
    # [cf, [rebin 1, rebin 2, ...]], same for unweighted if integrated analysis
    def get_cf(self):
        histos = []
        histos_unw = []
        cf_list = []
        cf_list_unw = []
        if self._atype == 1:        # integrated analysis
            histos, histos_unw = getIntegrated(self._se, self._me, self._htype, self._rebin)
            cf_list_unw.append(histos_unw[0][1])
            cf_list_unw.append([])
        elif self._atype == 2:      # differential analysis
            histos = getDifferential(self._se, self._me, self._htype, self._bins, self._rebin)

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

# class that handles the saving of the histograms in the correct file structure
class FemtoDreamSaver():
    last_edit = None

    def __init__(self, histos, conf):
        self._histos = histos   # all histos from the cf_handler
        self._ipath = conf[0]   # path to file
        self._oname = conf[1]   # name of file
        self._idir  = conf[2]   # TDirectory name
        self._atype = conf[3]   # analysis type
        self._htype = conf[4]   # histo type
        self._mc    = conf[5]   # bool monte carlo data
        self._bins  = conf[6]   # bin range for differential
        self._rebin = conf[7]   # rebin factors for all se, me, cf plots
        self._nfile = conf[8]   # "new", "recreate" or "update"
        self._opath = conf[9]   # output directory
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
        if self._nfile == "new":
            new_name = self._file_exists(self._opath + "UFFA_" + self._oname)
            print(new_name + " created!")
            self._oname = new_name
            ofile_name = self._opath + self._oname
            FemtoDreamSaver.last_edit = self._oname
        else:
            ofile_name = self._opath + "UFFA_" + self._oname

        if FemtoDreamSaver.last_edit:
            ofile_name = FemtoDreamSaver.last_edit
        ofile = ROOT.TFile(ofile_name, self._nfile)

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

        default = "femto-dream-pair-task-track-track"
        if self._idir == "" or self._idir == default:
            dir_root = ofile.mkdir(default + "_std")
        elif self._idir[0] == '_':
            dir_root = ofile.mkdir(default + self._idir)
        else:
            dir_root = ofile.mkdir(self._idir)
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
                self._mkdir_write(dir_root, "Tracks_one_MC", hist_track_mc)
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
                    del dir_bin
                del dir_mc
                self._mkdir_write(dir_root, "Tracks_one_MC", hist_track_mc)

# generates list with rebin factors
def bin2list(rebin):
    rebin_list = []
    if type(rebin) == int:
        rebin = [rebin]
    elif type(rebin) != list:
        print("rebin factor is not an int or a list!")
        exit()
    rebin_list.extend(rebin)
    return rebin_list

# returns list with the configured settings
# [path, file name, infile dir, analysis type, hist type, mc, bin range, rebin factors]
def config(dirOut = None, fname = None, fdir = None, new_file = None,
           atype = None, htype = None, mc = None, bins = None, rebin = None):
    k_keys = ['k', 'kstar', '1']
    mult_keys = ['mult', 'kmult', '2']
    mt_keys = ['mt', 'kmt', '3']
    int_keys = ['int', 'integrated', '1']
    dif_keys = ['diff', 'differential', '2']

    path_name = fname.rsplit('/', 1)
    if len(path_name) == 1:
        ipath = ""
        iname = path_name[0]
    else:
        ipath = path_name[0] + '/'
        iname = path_name[1]

    if ROOT.gSystem.AccessPathName(ipath + iname):
        print("file \"" + ipath + iname + "\" not found!")
        exit()

    if not fdir:
        idir = ""
    else:
        idir = fdir             # TDirectory in file

    if not dirOut:
        dirOut = ipath
    if ROOT.gSystem.AccessPathName(dirOut):
        print("output directory \"" + dirOut + "\" does not exist!")
        exit()

    if new_file == 1:
        new_file = "new"
    elif new_file == 2:
        new_file = "recreate"
    elif new_file == 3:
        new_file = "update"
    else:
        print("\nInput error:\tWrong 'new_file' option!\nOptions:\n\t1 -> \"new\"\n\t2 -> \"recreate\"\n\t3 -> \"update\"\n")
        exit()

    # analysis type
    atype = str(atype).lower()
    if atype in int_keys:
        atype = 1
    elif atype in dif_keys:
        atype = 2

    # histogram type
    htype = str(htype).lower()
    if htype in k_keys:
        htype = 1
    elif htype in mult_keys:
        htype = 2
    elif htype in mt_keys:
        htype = 3

    if mc:
        mc = True

    if rebin:
        rebin = bin2list(rebin)
    else:
        rebin = None

    return [ipath, iname, idir, atype, htype, mc, bins, rebin, new_file, dirOut]

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
# [[iSE, iME], [bin 1], [bin 1 rebin], [bin 2], [bin 2 rebin]...]
def getDifferential(se, me, htype, bins, rebin):
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
# or reweights and returns ([[iSE, iME], [se, me, cf]], [[me_unw, cf_unw]]) for kmult
# and [[iSE, iME], [se, me, cf]] for kmT
def getIntegrated(se, me, htype, rebin):
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

    histos.append(getCorrelation(se, me, "cf", ""))
    if rebin:               # append rebinned histos to list of histos
        histos_rebin = []
        for factor in rebin:
            se_rebin = rebin_hist(se, factor)
            me_rebin = rebin_hist(me, factor)
            rebin_conf = " rebin: " + str(factor)
            histos_rebin.append(getCorrelation(se_rebin, me_rebin, "rebin: " + str(factor), rebin_conf))
        histos.append(histos_rebin)
    if htype == 2:      # 2nd list with unweighted histos
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
