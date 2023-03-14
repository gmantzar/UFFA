import ROOT
import correlation_handler as gentlefemto


def getCorrelationFunction(infilename, config, hist_type, monte_carlo, rebinfactor):

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

    # set first bin factor to no rebinning and add other factors
    nobin = [1]
    if rebinfactor:
        if isinstance(rebinfactor, int):
            rebinfactor = [rebinfactor]
        nobin.extend(rebinfactor)
    rebinfactor = nobin

    th2 = ['th2f', 'th2', 'mult', 'kmult']
    hist_th2 = True if hist_type.lower() in th2 else False

    if hist_th2:
        if monte_carlo:
            se = directory.Get("SameEventMC/relPairkstarMult")
            me = directory.Get("MixedEventMC/relPairkstarMult")
        else:
            se = directory.Get("SameEvent/relPairkstarMult")
            me = directory.Get("MixedEvent/relPairkstarMult")
    else:
        if monte_carlo:
            se = directory.Get("SameEventMC/relPairDist")
            me = directory.Get("MixedEventMC/relPairDist")
        else:
            se = directory.Get("SameEvent/relPairDist")
            me = directory.Get("MixedEvent/relPairDist")

    se.SetDirectory(0)
    me.SetDirectory(0)

    infile.Close()

    histos = []
    ch = gentlefemto.CorrelationHandler("cf", se, me)
    if hist_th2:
        ch.reweight()
        hMult = [
                ["hSE_mult", ch.get_se().Clone()],
                ["hME_mult", ch.get_me().Clone()],
        ]
        histos.extend(hMult)

    for factor in rebinfactor:
        ch.rebin(factor)
        ch.make_cf()
        ch.normalize(0.24, 0.34)
        if hist_th2:
            ch.make_cf_unw()
            ch.normalize_unw(0.24, 0.34)

        conf_text = ""
        if factor > 1:
            conf_text += " rebin: " + str(factor)
        h1 = [
                ["hSE" + conf_text, ch.get_se_k().Clone()],
                ["hME" + conf_text, ch.get_me_k().Clone()],
                ["hCk" + conf_text, ch.get_cf().Clone()],
        ]
        histos.extend(h1)
        if hist_th2:
            h2 = [
                    ["hME_unw"+conf_text, ch.get_me_k_unw().Clone()],
                    ["hCk_unw"+conf_text, ch.get_cf_unw().Clone()],
            ]
            histos.extend(h2)
        ch.resetbin()

    return histos

def getSingleParticlePlots(infilename, config, monte_carlo = False):
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

    infile = ROOT.TFile(infilename, "read")
    directory = infile.GetDirectory("femto-dream-pair-task-track-track"+config)

    # list of TKey's in subdir
    subdir = directory.GetDirectory("Tracks_one")
    try:
        lobj = subdir.GetListOfKeys()
    except:
        print("Directory \"Tracks_one\" empty!")
        return
    lobj_ent = lobj.GetEntries()
    lnk = lobj.FirstLink()

    # loops over linked list of objects and saves histos
    histos = []
    for entry in range(lobj_ent):
        histos.append([lnk.GetObject().GetName(), lnk.GetObject().ReadObj()])
        lnk = lnk.Next()

    for name, hist in histos:
        hist.SetDirectory(0)

    if monte_carlo:
        subdir2 = directory.GetDirectory("Tracks_oneMC")
        lobjMC = subdir2.GetListOfKeys()
        lobj_entMC = lobjMC.GetEntries()
        lnkMC = lobjMC.FirstLink()
        histos2 = []
        for entry in range(lobj_entMC):
            histos2.append([lnkMC.GetObject().GetName(), lnkMC.GetObject().ReadObj()])
            lnkMC = lnkMC.Next()
        for name, hist in histos2:
            hist.SetDirectory(0)

    infile.Close()

    return (histos, histos2) if monte_carlo else histos

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

    infile = ROOT.TFile(infilename, "read")
    directory = infile.GetDirectory("femto-dream-pair-task-track-track"+config)
    subdir = directory.GetDirectory("Event")

    # list of TKey's in subdir
    try:
        lobj = subdir.GetListOfKeys()
    except:
        print("Directory \"Event\" empty!")
        return
    lobj_ent = lobj.GetEntries()
    lnk = lobj.FirstLink()

    # loops over linked list of objects and saves histos
    histos = []
    for entry in range(lobj_ent):
        histos.append([lnk.GetObject().GetName(), lnk.GetObject().ReadObj()])
        lnk = lnk.Next()

    for name, histo in histos:
        histo.SetDirectory(0)

    infile.Close()

    return histos

def saveHistogramms(path, filename, newfile, config, hist_type = 'TH1F', monte_carlo = False, rebin = None):
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

    # Get the histogramms from the input file
    correlationHistos = getCorrelationFunction(infile, config, hist_type, 0, rebin)

    if monte_carlo:
        hCorr_mc = getCorrelationFunction(infile, config, hist_type, 1, rebin)
        hTrack, hTrackMC = getSingleParticlePlots(infile, config, monte_carlo)
        hCorr_mc.extend(hTrackMC)
        for n in range(len(hCorr_mc)):
            name = hCorr_mc[n][0]
            if name == "hPt":
                hPt_mc = hCorr_mc[n][1]
    else:
        hTrack = getSingleParticlePlots(infile, config, monte_carlo)

    eventHistos = getEventHistos(infile, config)

    histos = []
    histos.extend(correlationHistos)
    try:
        histos.extend(hTrack)
        histos.extend(eventHistos)
    except:
        pass

    for n in range(len(histos)):
        name = histos[n][0]
        if name == "hPt":
            TPC_Int = histos[n][1].Integral(0, histos[n][1].FindBin(0.75))
            hPt = histos[n][1]
        if name == "zvtxhist":
            zvtx_Events = histos[n][1].GetEntries()


    # print here some interesting numbers like number of pairs in a certain k^* region, number of events, etc.
    print("\n### Getting file", filename, "###")
    print("Number of Events:\n", histos[0][1].GetEntries())
    try:
        print("nTPC/nEvents:\n", TPC_Int/zvtx_Events)
        print("Number of pairs in the low k* region:\n", histos[0][1].Integral(1, histos[0][1].FindBin(0.2)))
    except:
        pass


    fileaction = "recreate" if newfile else "update"

    outfile = path+"UFFA_"+filename
    output = ROOT.TFile(outfile, fileaction)

    if config == "":
        dir = output.mkdir("_std")
    else:
        dir = output.mkdir(config)

    dir.cd()

    for i in range(len(histos)):
        histos[i][1].Write(histos[i][0])

    # Subdir in same TDirectory for MC data
    if monte_carlo:
        hCorr_mc.extend(["hPt_ratio", ComputePurity(hPt, hPt_mc)])
        dir = dir.mkdir("MC")
        dir.cd()
        for i in range(len(hCorr_mc)):
            hCorr_mc[i][1].Write(hCorr_mc[i][0])

    output.Close()

def ComputePurity(hPt, hPt_mc):
    ratio = hPt.Clone("hPt_ratio")
    ratio.Divide(hPt_mc)
    return ratio.Clone()

