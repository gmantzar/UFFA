import ROOT

# expand ~ to home directory
def path_expand(ifile):
    if ifile[0] == '~':
        home = ROOT.gSystem.GetHomeDirectory()
        home += ifile[1:]
        return home
    else:
        return ifile

# fix dir path, i.e. remove preceding and trailing /
def path_fix(dir_name):
    if dir_name:
        if dir_name[-1] == '/':
            dir_name = dir_name[:-1]
        if dir_name[0] == '/':
            dir_name = dir_name[1:]
    else:
        dir_name = ""
    return dir_name

# function that looks if 'file' exists and if yes it append -n, where n depends on if another file was already created
def file_exists(file):
    file = path_expand(file)
    if not ROOT.gSystem.AccessPathName(file):
        name, ext = file.rsplit('.')
        digit = 1
        while not ROOT.gSystem.AccessPathName(name + '-' + str(digit) + '.' + ext):
            digit += 1
        file = name + '-' + str(digit) + '.' + ext
    return file
