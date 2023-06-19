import ROOT
import FileUtils as FU

class FileReader:
    DEBUG = False
    def __init__(self, ifile, directory = None):
        ifile = FU.path_expand(ifile)
        if not ROOT.gSystem.AccessPathName(ifile):
            self._ifile = ROOT.TFile(ifile, "read")
        else:
            print("File \"" + ifile + "\" not found!")
        self._tree  = [self._ifile]
        self._wdir  = self._ifile
        if directory:
            directory = FU.path_fix(directory)
            self._set_dir(directory)

    # find object in dir_obj
    def _find_obj(self, obj_name, dir_obj):
        obj = None
        try:                # try to find in TDirectory
            obj = dir_obj.Get(obj_name)
        except:
            pass
        if not obj:         # try to find in TList
            try:
                obj = dir_obj.FindObject(obj_name)
            except:
                pass
        return obj

    # set directory
    def _set_dir(self, dir_name):
        dir_name = FU.path_fix(dir_name)
        name_list = dir_name.rsplit('/')
        for name in name_list:
            dir_new = self._find_obj(name, self._wdir)              # find in current directory
            if not dir_new:
                dir_new = self._find_obj(name, self._ifile)         # find in root directory
                if dir_new:
                    self._tree = [self._ifile]                      # if found in root dir, reset working directory path
            if not dir_new:
                if FileReader.DEBUG:
                    print("\n\tDirectory \"" + name + "\" not found!\n")
                    self.ls()
                return None
            self._tree.append(dir_new)                              # append dir to dir path
            self._wdir = dir_new
            if dir_new.Class() == ROOT.TList.Class():
                continue
            dir_new.cd()

    # find and return object
    def _get_obj(self, obj_name):
        obj_name = FU.path_fix(obj_name)
        name_list = obj_name.rsplit('/')
        obj = self._find_obj(name_list[0], self._wdir)
        if not obj:
            obj = self._find_obj(name_list[0], self._ifile)
        if not obj:
            if FileReader.DEBUG:
                print("Object \"" + name_list[0] + "\" not found!")
            return None
        for name in name_list[1:]:
            obj = self._find_obj(name, obj)
            if not obj:
                if FileReader.DEBUG:
                    print("Object \"" + name + "\" not found!")
                return None
        return obj

    # retrieves a histogram by name in the current directory
    # or if given, from the full path or a subdirectory
    def get_histo(self, histo_name, dir_name = None):
        dir_name = FU.path_fix(dir_name)
        if not dir_name or dir_name == "":
            return self._get_obj(histo_name)
        return self._get_obj(dir_name + '/' + histo_name)

    def getHisto(self, histo_name, dir_name = None):
        dir_name = FU.path_fix(dir_name)
        if not dir_name or dir_name == "":
            return self._get_obj(histo_name)
        return self._get_obj(dir_name + '/' + histo_name)

    # function to retrieve all histograms in a directory as a list
    def get_histos(self, dir_name = None):
        if dir_name:
            dir_name = FU.path_fix(dir_name)             # fix the path
            directory = self._get_obj(dir_name)             # set directory
        else:
            directory = self._wdir                          # find in current directory
        histos = []
        if directory.Class() == ROOT.TList.Class():         # for TList
            lnk = directory.FirstLink()
            lobj_ent = directory.GetEntries()
        else:                                               # for TDirectory
            lobj = directory.GetListOfKeys()
            lobj_ent = lobj.GetEntries()
            lnk = lobj.FirstLink()
        for n in range(lobj_ent):                           # iterate over entries
            if directory.Class() == ROOT.TList.Class():
                obj = lnk.GetObject()
            else:
                obj = lnk.GetObject().ReadObj()
            if obj.InheritsFrom(ROOT.TH1.Class()):          # check if its a derivative of TH1
                histos.append(obj)
            lnk = lnk.Next()
        for hist in histos:
            hist.SetDirectory(0)
        return histos

    def getHistos(self, dir_name = None):
        if dir_name:
            dir_name = FU.path_fix(dir_name)             # fix the path
            directory = self._get_obj(dir_name)             # set directory
        else:
            directory = self._wdir                          # find in current directory
        histos = []
        if directory.Class() == ROOT.TList.Class():         # for TList
            lnk = directory.FirstLink()
            lobj_ent = directory.GetEntries()
        else:                                               # for TDirectory
            lobj = directory.GetListOfKeys()
            lobj_ent = lobj.GetEntries()
            lnk = lobj.FirstLink()
        for n in range(lobj_ent):                           # iterate over entries
            if directory.Class() == ROOT.TList.Class():
                obj = lnk.GetObject()
            else:
                obj = lnk.GetObject().ReadObj()
            if obj.InheritsFrom(ROOT.TH1.Class()):          # check if its a derivative of TH1
                histos.append(obj)
            lnk = lnk.Next()
        for hist in histos:
            hist.SetDirectory(0)
        return histos

    # function to retrieve a full directory as a list
    def get_dir(self, dir_name = None):
        if dir_name:
            dir_name = FU.path_fix(dir_name)
            directory = self._get_obj(dir_name)
        else:
            directory = self._wdir
        dir_content = []
        if directory.Class() == ROOT.TList.Class():
            lnk = directory.FirstLink
            lobj_ent = directory.GetEntries()
        else:
            lobj = directory.GetListOfKeys()
            lobj_ent = lobj.GetEntries()
            lnk = lobj.FirstLink()
        for n in range(lobj_ent):
            if directory.Class() == ROOT.TList.Class():
                obj = lnk.GetObject()
            else:
                obj = lnk.GetObject().ReadObj()
            dir_content.append(obj)
            lnk = lnk.Next()
        for entry in dir_content:
            entry.SetDirectory(0)
        return dir_content

    # function to change directory
    def cd(self, dir_name = None):
        if not dir_name and dir_name != 0:
            if len(self._tree) == 1:
                if FileReader.DEBUG:
                    print("Already in root directory!")
                    self._wdir.ls()
            else:
                self._tree.pop()                        # remove entry for working directory list
                self._wdir = self._tree[-1]
        else:
            # 0 is the root directory of the file
            if dir_name == 0:
                self._wdir = self._ifile
                self._tree = [self._ifile]              # reset working directory list
            else:
                self._set_dir(dir_name)

    # list directory content
    def ls(self):
        self._wdir.ls()

    # list file content
    def fls(self):
        self._ifile.ls()

    # return working directory
    def pwd(self):
        pwd = ""
        if len(self._tree) == 1:
            return self._ifile.GetName()
        for obj in self._tree[1:]:
            pwd += obj.GetName() + "/"
        if FileReader.DEBUG:
            print("pwd: \"" + pwd + "\"")
        return pwd

    # return current directory
    def get_dir(self):
        if FileReader.DEBUG:
            print("dir: \"" + self._wdir.GetName() + "\"")
        return self._wdir.GetName()

    # return file
    def get_file(self):
        if FileReader.DEBUG:
            print("file: \"" + self._ifile.GetName() + "\"")
        return self._ifile.GetName()

    # toggle debug output or set with 'option'
    def setDebug(option = None):
        if type(option) == bool:
            FileReader.DEBUG = option
            return
        FileReader.DEBUG = not FileReader.DEBUG

