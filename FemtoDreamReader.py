import ROOT
import FileReader as FR

class FemtoDreamReader(FR.FileReader):
    def __init__(self, ifile, directory = None):
        if directory == "":
            directory = "femto-dream-pair-task-track-track"
        elif directory[0] == '_':
            directory = "femto-dream-pair-task-track-track" + directory
        FR.FileReader.__init__(self, ifile, directory)

    ### Getter Functions ###
    def get_pt(self):
        return self.get_histo("Tracks_one/hPt")

    def get_pt_mc(self):
        return self.get_histo("Tracks_one_MC/hPt_ReconNoFake")

    def get_dca(self):
        return self.get_histo("Tracks_one/hDCAxy")

    def get_dca_mc(self):
        return [self.get_histo("Tracks_one_MC/hDCAxy_Primary"), \
                self.get_histo("Tracks_one_MC/hDCAxy_DaughterLambda"), \
                self.get_histo("Tracks_one_MC/hDCAxy_DaughterSigmaplus"), \
                self.get_histo("Tracks_one_MC/hDCAxy_Material")]

    def get_zvxt(self):
        return self.get_histo("Event/zvtxhist")

    def get_event(self):
        return self.get_histos("Event")

    def get_tracks(self):
        return self.get_histos("Tracks_one")

    def get_tracks_mc(self):
        return self.get_histos("Tracks_one_MC")

    def get_se(self):
        return self.get_histo("SameEvent/relPairDist")

    def get_me(self):
        return self.get_histo("MixedEvent/relPairDist")

    def get_kstar(self):
        return self.get_histo("SameEvent/relPairDist"), self.get_histo("MixedEvent/relPairDist")

    def get_kstar_mc(self):
        return self.get_histo("SameEvent_MC/relPairDist"), self.get_histo("MixedEvent_MC/relPairDist")

    def get_kmt(self):
        return self.get_histo("SameEvent/relPairkstarmT"), self.get_histo("MixedEvent/relPairkstarmT")

    def get_kmt_mc(self):
        return self.get_histo("SameEvent_MC/relPairkstarmT"), self.get_histo("MixedEvent_MC/relPairkstarmT")

    def get_kmult(self):
        return self.get_histo("SameEvent/relPairkstarMult"), self.get_histo("MixedEvent/relPairkstarMult")

    def get_kmult_mc(self):
        return self.get_histo("SameEvent_MC/relPairkstarMult"), self.get_histo("MixedEvent_MC/relPairkstarMult")
