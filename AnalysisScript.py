import FemtoAnalysis as FA

################################################################################################
#                                                                                              #
# UFFA_pp(path, fname, fdir, new_file, atype, htype, mc = False, bins = False, rebin = False): #
# new_file:                                                                                    #
#       1  ->  "new"                                                                           #
#       2  ->  "recreate"                                                                      #
#       3  ->  "update"                                                                        #
# atype:                                                                                       #
#       Integrated    ->  "Int",  "Integrated",   1                                            #
#       Differential  ->  "Diff", "Differential", 2                                            #
# htype:                                                                                       #
#       kstar  ->  "kstar", "k",    1                                                          #
#       kmult  ->  "kmult", "mult", 2                                                          #
#       kmt    ->  "kmt",   "mt",   3                                                          #
# mc:                                                                                          #
#       False (default)                                                                        #
#       True                                                                                   #
# bins:                                                                                        #
#       False (default)                                                                        #
#       mult (atype = 2):                                                                      #
#               [0, 20, 40, 60]                                                                #
#       mt (atype = 3):                                                                        #
#               [0.5, 1.5, 2.5, 4.5]                                                           #
# rebin:                                                                                       #
#       False (default)                                                                        #
#       int   ->  e.g. 5                                                                       #
#       list  ->  [2, 3, 7]                                                                    #
#                                                                                              #
################################################################################################

path = ""

#file = "AnalysisResults_22m_pass2_no-offset.root"
#file = "AnalysisResults.root"
#file = "AnalysisResults_18l.root"
#file = "input_CFOutput_pp.root"
file = "Task-MC-LHC21k6-protons-full-AnalysisResults.root"

#FA.UFFA_pp(path, file, "", 1, 1, 1, False, False, False)
#FA.UFFA_pp(path, file, "", 1, 1, 1, True, False, False)
#FA.UFFA_pp(path, file, "", 1, 1, 1, True, False, 5)
#FA.UFFA_pp(path, file, "", 1, 1, 1, True, False, [2, 3, 5])
#FA.UFFA_pp(path, file, "", 1, 2, 2, True, [0, 20, 40, 60], False)
#FA.UFFA_pp(path, file, "", 1, 2, 2, True, [0, 20, 40, 60], 3)
#FA.UFFA_pp(path, file, "", 1, 2, 2, True, [0, 20, 40, 60], [2, 3, 5])
#FA.UFFA_pp(path, file, "", 1, 2, 3, True, [0.5, 1.5, 2.5, 4.5], False)
#FA.UFFA_pp(path, file, "", 1, 2, 3, True, [0.5, 1.5, 2.5, 4.5], 3)
#FA.UFFA_pp(path, file, "", 1, 2, 3, True, [0.5, 1.5, 2.5, 4.5], [2, 3, 5])



