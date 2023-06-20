import ROOT
import FileUtils as FU
import FileSaver as FS

# class that handles the saving of the histograms in the correct file structure
class FemtoDreamSaver(FS.FileSaver):
    last_edit = None

    def __init__(self, settings, histos = None):
        FS.FileSaver.__init__(self)
        self._histos = histos
        self._func  = settings['function']          # save structure for 'function'
        self._oname = settings['rename'] if settings['rename'] \
                else  settings['file']              # output file name
        self._nfile = settings['newfile']           # "new", "recreate" or "update"
        self._ipath = settings['path']              # input directory
        self._opath = settings['outDir']            # output directory
        self._idir  = settings['fileTDir']          # file TDirectory
        self._atype = settings['atype']             # analysis type
        self._htype = settings['htype']             # histo type
        self._mc    = settings['mc']                # bool monte carlo data
        self._bins  = settings['bins']              # bin range for differential
        self._rebin = settings['rebin']             # rebin factors for all se, me, cf plots
        self._debug = settings['debug']
        FS.FileSaver.setDebug(self._debug)
        if self._func == 'cf':
            self._save_cf()
        elif self._func == 'tf':
            self._save_tf()
        else:
            self._save_cf()

    # create the output file
    def _create_output(self, prepend):
        if self._nfile == 'new':
            new_name = FU.file_exists(self._opath + prepend + self._oname)
            FemtoDreamSaver.last_edit = new_name
        elif self._nfile == 'recreate':
            new_name = self._opath + prepend + self._oname
            FemtoDreamSaver.last_edit = None
        elif self._nfile == 'update':
            if FemtoDreamSaver.last_edit:
                new_name = FemtoDreamSaver.last_edit
            else:
                new_name = self._opath + prepend + self._oname
        self.touch(new_name, self._nfile)

    # function that saves all the histos in the correct file format
    # saver for UFFA_cf()
    def _save_cf(self):
        # create the output file
        self._create_output("UFFA_")

        # histos setup
        hist_std = self._histos[0]
        hist_in = hist_std[0]               # [iSE, iME]
        hist_smc = hist_std[1:]             # [[se, me, cf, [rebin]], ...]
        if self._mc:
            hist_std_mc = self._histos[2]
            hist_in_mc  = hist_std_mc[0]    # [iSE, iME] <- MC
            hist_smc_mc = hist_std_mc[1:]   # [[se, me, cf, [rebin]], ...] <- MC
        if self._htype == 'mult':
            hist_unw = self._histos[1]      # [me unw, cf unw, [rebin]]
            if self._mc:
                hist_unw_mc = self._histos[3]

        hist_event = self._histos[4]
        hist_track = self._histos[5]
        if self._mc:
            hist_track_mc = self._histos[6]
            hist_pur = getPurity(hist_track[0], hist_track_mc[0])

        default = "femto-dream-pair-task-track-track"
        if self._idir == default:
            dir_root = self._file.mkdir(default + "_std")
        else:
            dir_root = self._file.mkdir(self._idir)
        dir_root.cd()

        if self._atype == 'int':                # integrated
            self.write(hist_in)                           # iSE, iME
            self.write(hist_smc[0][:-1])                  # se, me, cf
            if self._htype == 'mult':
                self.write(hist_unw[:-1])                 # me unw, cf unw
            if self._rebin:                         # all rebinned se, me, cf, me unw, cf unw, in individual directories
                for n in range(len(self._rebin)):
                    dir_rebin = dir_root.mkdir("rebin: " + str(self._rebin[n]))
                    dir_rebin.cd()
                    self.write(hist_smc[0][3][n])         # [[se, me, cf, [rebins]]]
                    if self._htype == 'mult':
                        self.write(hist_unw[2][n])        # [me unw, cf unw, [rebins]]
                    dir_root.cd()
                    del dir_rebin
            self.writeInDir(dir_root, "Event", hist_event)
            self.writeInDir(dir_root, "Tracks_one", hist_track)
            if self._mc:
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self.write(hist_in_mc)                    # [iSE, iME]
                self.write(hist_smc_mc[0][:-1])           # [[se, me, cf, [rebins]]
                if self._htype == 'mult':
                    self.write(hist_unw_mc[:-1])          # [me unw, cf unw, [rebins]]
                if self._rebin:                     # rebin directories in main directory
                    for n in range(len(self._rebin)):
                        dir_rebin = dir_mc.mkdir("rebin: " + str(self._rebin[n]))
                        dir_rebin.cd()
                        self.write(hist_smc_mc[0][3][n])
                        if self._htype == 'mult':
                            self.write(hist_unw_mc[2][n])
                        dir_mc.cd()
                        del dir_rebin
                self.writeInDir(dir_root, "Tracks_one_MC", hist_track_mc)
                hist_pur.Write()
        elif self._atype == 'dif':              # differential
            self.write(hist_in)                           # [iSE, iME]
            for n in range(len(self._bins) - 1):        # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                dir_bin = dir_root.mkdir("bin: " + str(n + 1))
                dir_bin.cd()
                if self._rebin:                     # rebin directories inside each mt/mult bin directory
                    self.write(hist_smc[n][:-1])          # [[se, me, cf, [rebins]], ...]
                    for m in range(len(self._rebin)):
                        dir_rebin = dir_bin.mkdir("rebin: " + str(self._rebin[m]))
                        dir_rebin.cd()
                        self.write(hist_smc[n][3][m])     # [[se, me, cf, [rebins]], ...]
                        dir_bin.cd()
                        del dir_rebin
                dir_root.cd()
                del dir_bin
            self.writeInDir(dir_root, "Event", hist_event)
            self.writeInDir(dir_root, "Tracks_one", hist_track)
            if self._mc:                            # monte carlo directory
                dir_mc = dir_root.mkdir("mc")
                dir_mc.cd()
                self.write(hist_in_mc)                    # [iSE, iME]
                for n in range(len(self._bins) - 1):    # list: [1, 2, 3, 4] -> ranges: [1-2, 2-3, 3-4]
                    dir_bin = dir_mc.mkdir("bin: " + str(n + 1))
                    dir_bin.cd()
                    if self._rebin:                 # rebin directories inside each mt/mult bin directory
                        self.write(hist_smc_mc[n][:-1])
                        for m in range(len(self._rebin)):
                            dir_rebin = dir_bin.mkdir("rebin: " + str(self._rebin[m]))
                            dir_rebin.cd()
                            self.write(hist_smc_mc[n][3][m])
                            dir_bin.cd()
                            del dir_rebin
                    del dir_bin
                del dir_mc
                self.writeInDir(dir_root, "Tracks_one_MC", hist_track_mc)

    # saver for UFFA_tf()
    def _save_tf(self):
        # create the output file
        self._create_output("TemplateFit_")
