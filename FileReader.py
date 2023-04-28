import ROOT

class FileReader:

    def __init__(self, ifile, directory = None):
        self._ifile = ROOT.TFile(ifile, "read")
        self._dir = self._ifile.GetDirectory(directory) if directory else self._ifile

    # find and set the current working directory
    def _set_dir(self, dir_name):
        if self._dir.GetDirectory(dir_name):
            self._dir = self._dir.GetDirectory(dir_name)
        elif self.__file.GetDirectory(dir_name):
            self._dir = self.__file.GetDirectory(dir_name)
        else:
            print("TDirectory \""+dir_name+"\" not found!")

    # find and return a TDirectory by name
    def _set_localdir(self, dir_name):
        if self._dir.GetDirectory(dir_name):
            directory = self._dir.GetDirectory(dir_name)
        elif self._ifile.GetDirectory(dir_name):
            directory = self._ifile.GetDirectory(dir_name)
        else:
            return None
        return directory

    # retrieves a histogram by name in the current directory
    # or if given, from the full path or a subdirectory
    def get_histo(self, histo_name, dir_name = None):
        if dir_name:
            directory = self._set_localdir(dir_name)
            if not directory:
                print("get_histo: TDirectory '"+dir_name+"' not found!")
                return
            histo = directory.Get(histo_name)
        else:
            histo = self._dir.Get(histo_name)
        if not histo:
            print("get_histo: \""+histo_name+"\" not found!")
            return
        histo.SetDirectory(0)
        return histo

    # function to retrieve all histograms in a directory as a list
    def get_histos(self, dir_name = None):
        if dir_name:
            directory = self._set_localdir(dir_name)
        histos = []
        if directory:
            lobj = directory.GetListOfKeys()
            lobj_ent = lobj.GetEntries()
            lnk = lobj.FirstLink()
            for n in range(lobj_ent):
                obj = lnk.GetObject().ReadObj()
                if obj.InheritsFrom(ROOT.TH1.Class()):
                    histos.append(lnk.GetObject().ReadObj())
                lnk = lnk.Next()
            for hist in histos:
                hist.SetDirectory(0)
        else:
            print("get_histos: TDirectory '"+dir_name+"' not found!")
            return
        return histos

    # function to retrieve a full directory as a list
    def get_dir(self, dir_name = None):
        if dir_name:
            directory = self._set_localdir(dir_name)
        dir_content = []
        if directory:
            lobj = directory.GetListOfKeys()
            lobj_ent = lobj.GetEntries()
            lnk = lobj.FirstLink()
            for n in range(lobj_ent):
                obj = lnk.GetObject().ReadObj()
                dir_content.append(lnk.GetObject().ReadObj())
                lnk = lnk.Next()
            for entry in dir_content:
                entry.SetDirectory(0)
        else:
            print("get_histos: TDirectory '"+dir_name+"' not found!")
            return
        return dir_content

    # function to chance directory
    def cd(self, dir_name = None):
        # 0 is the root directory of the file
        if dir_name:
            if dir_name == 0:
                self._dir = self._ifile
            else:
                self._set_dir(dir_name)
        else:
            if self._dir.GetMotherDir():
                self._dir = self._dir.GetMotherDir()

    # list directory content
    def ls(self):
        self._dir.ls()

    # list file content
    def fls(self):
        self._ifile.ls()

    # print working directory
    def pwd(self):
        self._dir.pwd()
