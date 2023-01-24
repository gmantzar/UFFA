import FemtoAnalysis as FemtoAnalysis


# Here, I use the saveHistogramms function for every configuration and file I want to analyse
# Note, that the bool variable newfile is true for the base configuration and has to be set to false
# for the variations within the same dataset in order to save the output in the same file


########## LHC18l ##########
thepath = "~/work/projects/pp_python/"
#thefilename = "AnalysisResults_18l.root"
thefilename = "AnalysisResults_22m_pass2_no-offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, True, "", 'th2')
FemtoAnalysis.saveHistogramms(thepath, thefilename, False, "_ap-base", 'mult', [2, 4])
#FemtoAnalysis.saveHistogramms(thepath, thefilename, "_clusters100", 2, False)
#FemtoAnalysis.saveHistogramms(thepath, thefilename, "_eta07", 2, False)
#FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-clusters100", 2, False)
#FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-eta07", 2, False)
