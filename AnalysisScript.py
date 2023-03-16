import FemtoAnalysis as FemtoAnalysis


# Here, I use the saveHistogramms function for every configuration and file I want to analyse
# Note, that the bool variable newfile is true for the base configuration and has to be set to false
# for the variations within the same dataset in order to save the output in the same file


########## LHC18l ##########
thepath = "~/work/projects/pp_python/"
#thefilename = "AnalysisResults_22m_pass2_no-offset.root"
#thefilename = "AnalysisResults.root"
thefilename = "Merged_AnalysisResults.root"
#thefilename = "AnalysisResults_18l.root"
#thefilename = "input_CFOutput_pp.root"
#thefilename2 = "Task-MC-LHC21k6-protons-full-AnalysisResults.root"
# saveHistogramms(filepath, filename, new_output, TDir_name, hist_type, monte_carlo, binning, rebin)
#FemtoAnalysis.saveHistogramms(thepath, thefilename, True, "", 5, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, True, "", "mult", False, [0.5, 1.5, 2.5, 4.5])
FemtoAnalysis.saveHistogramms(thepath, thefilename, False, "_apap", "mult", False, [0.5, 1.5, 2.5, 4.5])
#FemtoAnalysis.saveHistogramms(thepath, thefilename2, True, "", "mult", True)
#FemtoAnalysis.saveHistogramms(thepath, thefilename, False, "_ap-base", "mult", False, [2, 4])
