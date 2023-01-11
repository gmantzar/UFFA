import FemtoAnalysis


# Here, I use the saveHistogramms function for every configuration and file I want to analyse
# Note, that the bool variable newfile is true for the base configuration and has to be set to false
# for the variations within the same dataset in order to save the output in the same file


########## LHC22m_pass1_subset ##########

#LHC22m_pass1_subset_offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22m_pass1_subset/offset/"
thefilename = "AnalysisResults_22m_pass1_subset_offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 2, False)

#LHC22m_pass1_subset_no-offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22m_pass1_subset/no-offset/"
thefilename = "AnalysisResults_22m_pass1_subset_no-offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)

########## LHC22m_pass2_subset ##########

#LHC22m_pass1_subset_offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22m_pass2_subset/"
thefilename = "AnalysisResults_22m_pass2_subset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 2, False)

########## LHC22r_pass1_subset ##########

#LHC22r_pass1_subset_offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22r_pass1_subset/offset/"
thefilename = "AnalysisResults_22r_pass1_subset_offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 5, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 5, False)

#LHC22r_pass1_subset_no-offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22r_pass1_subset/no-offset/"
thefilename = "AnalysisResults_22r_pass1_subset_no-offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 5, True)


########## LHC22p_pass1 ##########

#LHC22p_pass1_offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22p_pass1/offset/"
thefilename = "AnalysisResults_22p_pass1_offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 2, False)

#LHC22p_pass1_no-offset
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22p_pass1/no-offset/"
thefilename = "AnalysisResults_22p_pass1_no-offset.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)

#LHC22f_pass2
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC22f_pass2/"
thefilename = "AnalysisResults_22f_pass2.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 5, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_dca35", 5, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 5, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-dca35", 5, False)

########## LHC18l ##########
thepath = "~/femto/pprun3/TrainOutput/run20221230/LHC18l/"
thefilename = "AnalysisResults_18l.root"
FemtoAnalysis.saveHistogramms(thepath, thefilename, "", 2, True)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_clusters100", 2, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_eta07", 2, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-base", 2, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-clusters100", 2, False)
FemtoAnalysis.saveHistogramms(thepath, thefilename, "_ap-eta07", 2, False)
