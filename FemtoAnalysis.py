import ROOT
import correlation_handler as gentlefemto


def getCorrelationFunction(infilename, config, hist_type, rebinfactor):

    """
    Function to get the same event (se) and mixed event (me) distribution and calculate the
    correlation function using the correlation_handler class.

    Parameters:

    infilename = path and name to the AnalysisResults.root file, where the same event and mixed event distrubution is

    config = configuration name of the subwagon. Use an empty string ("") for the standard configuration.
             use an underscore before the name of another configuration (e.g. "_ap-base")

    rebinfactor = rebinning factor to be applied ONLY to the correlation function within the correlation_handler class

    Output:

    array with the three histograms in the following order: se, me, cf
    """
    infile = ROOT.TFile(infilename,"read")
    directory = infile.GetDirectory("femto-dream-pair-task-track-track"+config)
    if not directory:
        print("TDirectory femto-dream-pair-task-track-track" + config + " does not exists!")
        exit()

    if isinstance(rebinfactor, int):
        rebinfactor = [rebinfactor]

    th2 = ['TH2', 'th2', 'mult', 'kmult']
    if hist_type in th2:
        se = directory.Get("SameEvent/relPairkstarMult")
        me = directory.Get("MixedEvent/relPairkstarMult")
    else:
        se = directory.Get("SameEvent/relPairDist")
        me = directory.Get("MixedEvent/relPairDist")

    se.SetDirectory(0)
    me.SetDirectory(0)

    infile.Close()

    ch = gentlefemto.CorrelationHandler("cf", se, me)
    if hist_type in th2:
        ch.reweight()

    histos = []
    for factor in rebinfactor:
        ch.rebin(factor)
        ch.make_cf()
        ch.normalize(0.24, 0.34)
        if hist_type in th2:
            ch.make_cf_unw()
            ch.normalize_unw(0.24, 0.34)

        h1 = [
                ["hSE rebin: "+str(factor), ch.get_se_k().Clone()],
                ["hME rebin: "+str(factor), ch.get_se_k().Clone()],
                ["hCF rebin: "+str(factor), ch.get_cf().Clone()],
        ]
        histos.extend(h1)
        if hist_type in th2:
            h2 = [
                    ["hME_unw rebin: "+str(factor), ch.get_me_k_unw().Clone()],
                    ["hCF_unw rebin: "+str(factor), ch.get_cf_unw().Clone()],
            ]
            histos.extend(h2)

        ch.resetbin()

    return histos

def getSingleParticlePlots(infilename, config):
    """
    Function to get the histograms of the single particle properties, i.e.
     - transverse momentum pT
     - pseudorapidity eta
     - angle phi
     - transverse DCA

    Parameters:

    infilename = path and name to the AnalysisResults.root file, where the same event and mixed event distrubution is

    config = configuration name of the subwagon. Use an empty string ("") for the standard configuration.
             use an underscore before the name of another configuration (e.g. "_ap-base")

    Output:

    array with the four histograms in the following order: hPt, hEta, hPhi, hDCAxy
    """

    directory = "femto-dream-pair-task-track-track"+config

    infile = ROOT.TFile(infilename, "read")
    histos = []


    histos.append(["hPt", infile.Get(directory+"/Tracks_one/hPt")])
    histos.append(["hEta", infile.Get(directory+"/Tracks_one/hEta")])
    histos.append(["hPhi", infile.Get(directory+"/Tracks_one/hPhi")])
    histos.append(["hDCAxy", infile.Get(directory+"/Tracks_one/hDCAxy")])
    for name, hist in histos:
        hist.SetDirectory(0)

    infile.Close()

    return histos

def getEventHistos(infilename, config):
    """
    Function to get the histograms of the event properties, i.e.
     - position distribution of the primary vertex in z direction
     - amplitude measured in the V0 detector (used as a multiplicity estimator)
     - number of tracks contributing to finding the primary vertex (i.e. the multiplicity at the primary vertex)

    Parameters:

    infilename = path and name to the AnalysisResults.root file, where the same event and mixed event distrubution is

    config = configuration name of the subwagon. Use an empty string ("") for the standard configuration.
             use an underscore before the name of another configuration (e.g. "_ap-base")

    Output:

    array with the three histograms in the following order: zvtxhist, MultV0M, MultNTracksPV
    """

    directory = "femto-dream-pair-task-track-track"+config

    infile = ROOT.TFile(infilename, "read")
    histos = []

    histos.append(["zvtxhist", infile.Get(directory+"/Event/zvtxhist")])
    histos.append(["MultV0M", infile.Get(directory+"/Event/MultV0M")])
    histos.append(["MultNTracksPV", infile.Get(directory+"/Event/MultNTracksPV")])
    for name, hist in histos:
        hist.SetDirectory(0)

    infile.Close()

    return histos

def saveHistogramms(path, filename, newfile, config, hist_type = 'TH1F', rebin = 2):
    """
    Function to save the histogramms of an "AnalysisResults.root" file to a new file.
    Its name is the name of the infput file with the prefix "GF-output_".
    It saves the output of the specified subwagon to an individual directory withing the file.
    If you want the output of more substrains within the same output file, you have to call the function
    multiple times and ster the parameter "newfile" to true only for the first time

    Parameters:

    path = path to the AnalysisResults.root file and the path, where the new root file will be stored

    filename = filename of the AnalysisResults.root file

    config = configuration name of the subwagon. Use an empty string ("") for the standard configuration.
             use an underscore before the name of another configuration (e.g. "_ap-base")

    rebin = rebinning factor to be applied ONLY to the correlation function within the correlation_handler class

    newfile = bool value to create a new file or not.
              -> for true the outputfile is called with the option "recreate"
              -> for false the outputfile is created with the option "update"

    Output:

    no output.. just saving the new root file
    """

    infile = path+filename

    #Get the histogramms from the input file
    correlationHistos = getCorrelationFunction(infile, config, hist_type, rebin)
    trackHistos = getSingleParticlePlots(infile, config)
    eventHistos = getEventHistos(infile, config)

    histos = []
    histos.extend(correlationHistos)
    histos.extend(trackHistos)
    histos.extend(eventHistos)
    for n in range(len(histos)):
        name = histos[n][0]
        if name == "hPt":
            nTPC = histos[n][1].Integral(0, histos[n][1].FindBin(0.75))
        if name == "zvtxhist":
            nEvents = histos[n][1].GetEntries()


    #print here some interesting numbers like number of pairs in a certain k^* region, number of events, etc.
    print(" ")
    print("### Getting file "+filename+" with the configuration "+config+" ###")

    print("Number of Events:")
    print(histos[0][1].GetEntries())

    print("nTPC/nEvents")
    print(nTPC/nEvents)

    print("Number of pairs in the low k* region")
    print(histos[0][1].Integral(1, histos[0][1].FindBin(0.2)))



    if newfile:
        fileaction = "recreate"
    else:
        fileaction = "update"

    outfile = path+"UFFA_"+filename
    output = ROOT.TFile(outfile,fileaction)
    if config=="":
        dir = output.mkdir("_std")
    else:
        dir = output.mkdir(config)

    dir.cd()

    for i in range(len(histos)):
        histos[i][1].Write(histos[i][0])

    output.Close()
