import ROOT
import FileReader as FR

class FemtoDreamReader(FR.FileReader):
    def __init__(self, ifile, directory = None):
        if directory == "":
            directory = "femto-dream-pair-task-track-track"
        elif directory[0] == '_':
            directory =="femto-dream-pair-task-track-track" + directory
        FR.FileReader.__init__(self, ifile, directory)

    ### Getter Functions ###
    def get_pt(self):
        return self.get_histo("hPt", "Tracks_one")

    def get_pt_mc(self):
        return self.get_histo("hPt", "Tracks_one_MC")

    def get_dca(self):
        return self.get_histo("hDCAxy", "Tracks_one")

    def get_zvxt(self):
        return self.get_histo("zvtxhist", "Event")

    def get_event(self):
        return self.get_histos("Event")

    def get_tracks(self):
        return self.get_histos("Tracks_one")

    def get_tracks_mc(self):
        return self.get_histos("Tracks_one_MC")

    def get_se(self):
        return self.get_histo("relPairDist", "SameEvent")

    def get_me(self):
        return self.get_histo("relPairDist", "MixedEvent")

    def get_kstar(self):
        return self.get_histo("relPairDist", "SameEvent"), self.get_histo("relPairDist", "MixedEvent")

    def get_kstar_mc(self):
        return self.get_histo("relPairDist", "SameEvent_MC"), self.get_histo("relPairDist", "MixedEvent_MC")

    def get_kmt(self):
        return self.get_histo("relPairkstarmT", "SameEvent"), self.get_histo("relPairkstarmT", "MixedEvent")

    def get_kmt_mc(self):
        return self.get_histo("relPairkstarmT", "SameEvent_MC"), self.get_histo("relPairkstarmT", "MixedEvent_MC")

    def get_kmult(self):
        return self.get_histo("relPairkstarMult", "SameEvent"), self.get_histo("relPairkstarMult", "MixedEvent")

    def get_kmult_mc(self):
        return self.get_histo("relPairkstarMult", "SameEvent_MC"), self.get_histo("relPairkstarMult", "MixedEvent_MC")
