import FemtoAnalysis as FA
import FemtoDreamReader as FDR

#########################################################################################################
#                                                                                                       #
# UFFA_pp(path, fname, fdir, new_file, atype, htype, mc = False, bins = False, rebin = False, dirOut):  #
# new_file:                                                                                             #
#       1  ->  "new"                                                                                    #
#       2  ->  "recreate"                                                                               #
#       3  ->  "update"                                                                                 #
# atype:                                                                                                #
#       Integrated    ->  "Int",  "Integrated",   1                                                     #
#       Differential  ->  "Diff", "Differential", 2                                                     #
# htype:                                                                                                #
#       kstar  ->  "kstar", "k",    1                                                                   #
#       kmult  ->  "kmult", "mult", 2                                                                   #
#       kmt    ->  "kmt",   "mt",   3                                                                   #
# mc:                                                                                                   #
#       False (default)                                                                                 #
#       True                                                                                            #
# bins:                                                                                                 #
#       False (default)                                                                                 #
#       mult (atype = 2):                                                                               #
#               [0, 20, 40, 60]                                                                         #
#       mt (atype = 3):                                                                                 #
#               [0.5, 1.5, 2.5, 4.5]                                                                    #
# rebin:                                                                                                #
#       False (default)                                                                                 #
#       int   ->  e.g. 5                                                                                #
#       list  ->  [2, 3, 7]                                                                             #
# dirOut:                                                                                               #
#       default -> dirIn                                                                                #
#                                                                                                       #
#########################################################################################################

ipath = "root_input/"
opath = "root_output/"

file = "AnalysisResults.root"
#file = "AnalysisResults_22m_pass2_no-offset.root"
#file = "AnalysisResults_18l.root"
#file = "input_CFOutput_pp.root"
#file = "Task-MC-LHC21k6-protons-full-AnalysisResults.root"
#file = "input_CFOutput_pp.root"
#file = "newest_AnalysisResults.root"
#file = "run2dca-UFFAinput.root"

templates = "AnalysisResults_forTemplates.root"
filename0 = "AnalysisResults_22f.root"
filename1 = "AnalysisResults_22q.root"
filename2 = "AnalysisResults_22o.root"
filename3 = "AnalysisResults_LHC22f_pass3.root"
filename4 = "AnalysisResults_LHC22o_pass3_small.root"
filename5 = "AnalysisResults_22qpass3.root"
filename6 = "AnalysisResults_LHC21k6_pp.root"

fdir = "femto-dream-pair-task-track-track_base-tempFit"
Out = ""

dca_names0 = ['prim', 'lam', 'sig', 'mat']
#dca_names1 = ['prim', 'lam', 'sig', 'mat', 'sec', 'fake']
dca_names1 = ['prim', 'lam', 'sig', 'mat', 'fake']
dca_names2 = ['prim', 'mat', 'notprim', 'fake']


f1 = FDR.FemtoDreamReader(templates, "")
#mc_plots = [f1.get_histo("Tracks_one_MC/hDCAxy_Primary"),
        #f1.get_histo("Tracks_one_MC/hDCAxy_Material"), \
        #f1.get_histo("Tracks_one_MC/hDCAxy_NotPrimary"), \
        #f1.get_histo("Tracks_one_MC/hDCAxy_Fake")]

mc_plots0 = f1.get_dca_mc()
#mc_plots0.append(f1.get_histo("Tracks_one_MC/hDCAxy_Daughter"))
mc_plots0.append(f1.get_histo("Tracks_one_MC/hDCAxy_Fake"))

#FA.TemplateFit("dca", 0.1, 15, dca_names1, filename6, "", mc_plots0)
#FA.TemplateFit("dca", 0.1, [0.5, 1.0, 2.0, 4], dca_names1, filename6, "", mc_plots0)
#FA.TemplateFit("dca", 0.1, [[0.5, 1.0, 2.0, 4], 2], dca_names1, filename6, "", mc_plots0)
#FA.TemplateFit("dca", dca_names0, filename1, "", filename0, "_base-tempFit")

#FA.UFFA_pp(file, "", 1, 1, 1, dirOut = opath)
#FA.UFFA_pp(file, "", 1, 2, 2, True, [0, 20, 40, 60], dirOut = opath)
#FA.UFFA_pp(ipath + file, "", 2, 1, 2, False, False, False, dirOut = opath)
#FA.UFFA_pp(file, "_apap", 3, 1, 2, False, False, False, opath)
#FA.UFFA_pp(file, "", 1, 1, 1, True, False, False)
#FA.UFFA_pp(file, "", 1, 1, 1, True, False, 5)
#FA.UFFA_pp(file, "", 1, 1, 1, True, False, [2, 3, 5])
#FA.UFFA_pp(file, "", 1, 2, 2, False, [0, 20, 40, 60], False)
#FA.UFFA_pp(file, "", 1, 2, 2, False, [0, 20, 40, 60], 3)
#FA.UFFA_pp(file, "", 1, 2, 2, False, [0, 20, 40, 60], [2, 3, 5])
#FA.UFFA_pp(file, "", 1, 2, 3, False, [0.5, 1.5, 2.5, 4.5], False)
#FA.UFFA_pp(file, "", 1, 2, 3, True, [0.5, 1.5, 2.5, 4.5], 3)
#FA.UFFA_pp(file, "", 1, 2, 3, False, [0.5, 1.5, 2.5, 4.5], [2, 3, 5])



