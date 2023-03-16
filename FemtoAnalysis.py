import ROOT
import correlation_handler as gentlefemto
from math import ceil


def getCorrelationFunction(infilename, config, hist_type, monte_carlo, binning, rebinfactor):

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
        #if isinstance(rebinfactor, int):
        if type(rebinfactor) == int:
            rebinfactor = [rebinfactor]
        nobin.extend(rebinfactor)
    rebinfactor = nobin

    mult_keywords = ['mult', 'kmult']
    mt_keywords = ['mt', 'mtk', 'kmt']
    hist_bin = False
    hist_mult = False
    hist_mt = False

    if type(binning) == int or type(binning) == list:
        hist_bin = True
        mt_bins = binning
    if hist_type.lower() in mult_keywords:
        hist_mult = True
        conf = "mult: "
    elif hist_type.lower() in mt_keywords:
        hist_mt = True
        conf = "mt: "

    if hist_mult:
        if monte_carlo:
            se = directory.Get("SameEventMC/relPairkstarMult")
            me = directory.Get("MixedEventMC/relPairkstarMult")
        else:
            se = directory.Get("SameEvent/relPairkstarMult")
            me = directory.Get("MixedEvent/relPairkstarMult")
    elif hist_mt:
        if monte_carlo:
            se = directory.Get("SameEventMC/relPairkstarmT")
            me = directory.Get("MixedEventMC/relPairkstarmT")
        else:
            se = directory.Get("SameEvent/relPairkstarmT")
            me = directory.Get("MixedEvent/relPairkstarmT")
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

    # correlation handler
    # mt differential
    if hist_bin:
        mt_histos = getBinRangeHistos(se, me, mt_bins)
        histo_lol = []
        for name, histo in mt_histos:
            ch = gentlefemto.CorrelationHandler("cf", histo[0], histo[1])
            ch.make_cf()
            ch.normalize(0.24, 0.34)
            histo_lol.extend([[conf + name, ch.get_se().Clone(), ch.get_me().Clone(), ch.get_cf().Clone()]])
            del ch
        return histo_lol
    # mt integrated
    else:
        ch = gentlefemto.CorrelationHandler("cf", se, me)
        histos = []

        if hist_mult:
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
            if hist_mult:
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
            if hist_mult:
                h2 = [
                        ["hME_unw"+conf_text, ch.get_me_k_unw().Clone()],
                        ["hCk_unw"+conf_text, ch.get_cf_unw().Clone()],
                ]
                histos.extend(h2)
            ch.resetbin()

    return histos

def getBinRangeHistos(inputSE, inputME, bins):

    inputSE.SetDirectory(0)
    inputME.SetDirectory(0)

    xAxis = inputSE.GetXaxis()
    yAxis = inputSE.GetYaxis()

    # generate list of bins from amount of mt bins or list of mt ranges
    if type(bins) == int:
        limits = [binLow]
        binLow = 1
        binHigh = inputSE.FindLastBinAbove(0, 2)    # last bin in Y axis with more than 0 counts
        binWidth = ceil((binHigh - binLow) / bins)
        for n in range(bins):
            limits.append((n + 1) * binWidth)
    elif type(bins) == list:
        limits = []
        for mt_value in bins:
            bin_value = yAxis.FindBin(float(mt_value))
            limits.append(bin_value)
    else:
        print("Error in getBinRangeHistos: bin input \"" + str(bins) + "\"")

    histos = []
    for n in range(1, len(limits)):
        name = "%.2f-%.2f" % (yAxis.GetBinCenter(limits[n - 1]), yAxis.GetBinCenter(limits[n]))
        se = inputSE.ProjectionX("se_k", limits[n - 1], limits[n])
        me = inputME.ProjectionX("me_k", limits[n - 1], limits[n])
        se.SetDirectory(0)
        me.SetDirectory(0)
        histos.append([name, [se.Clone(), me.Clone()]])

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
    if monte_carlo:
        subdir = directory.GetDirectory("Tracks_oneMC")
    else:
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

def getPurity(hPt, hPt_mc):
    ratio = hPt.Clone("hPt_ratio")
    ratio.Divide(hPt_mc)
    return ratio.Clone()

def saveHistogramms(path, filename, newfile, config, hist_type = 'TH1F', monte_carlo = False, binning = None, rebin = None):
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
    hCorr = getCorrelationFunction(infile, config, hist_type, 0, binning, rebin)
    hTrack = getSingleParticlePlots(infile, config, 0)
    hEvent = getEventHistos(infile, config)

    if monte_carlo:
        hCorrMC = getCorrelationFunction(infile, config, hist_type, 1, binning, rebin)
        hTrackMC = getSingleParticlePlots(infile, config, 1)

        for name, histo in (hCorrMC + hTrackMC):
            if name == "hPt":
                hPt_mc = histo

    histos = []
    histos.extend(hCorr)
    try:
        histos.extend(hTrack)
        histos.extend(hEvent)
        for name, histo in (hTrack + hEvent):
            if name == "hPt":
                TPC_Int = histo.Integral(0, histo.FindBin(0.75))
                hPt = histo
            if name == "zvtxhist":
                zvtx_Events = histo.GetEntries()
    except:
        pass

    # save file
    fileaction = "recreate" if newfile else "update"
    outfile = path+"UFFA_"+filename
    output = ROOT.TFile(outfile, fileaction)

    if config == "":
        dir = output.mkdir("_std")
    else:
        dir = output.mkdir(config)
    dir.cd()

    if binning:
        for name, se, me, cf in hCorr:
            mt_dir = dir.mkdir(name)
            mt_dir.cd()
            se.Write("hSE " + name)
            me.Write("hME " + name)
            cf.Write("hCF " + name)
            dir.cd()
        for name, histo in (hTrack + hEvent):
            histo.Write(name)
    else:
        for name, histo in (hCorr + hTrack + hEvent):
            histo.Write(name)

    # Subdir in same TDirectory for MC data
    if monte_carlo:
        dir_mc = dir.mkdir("MC")
        dir_mc.cd()
        for name, histo in (hCorrMC + hTrackMC):
            histo.Write(name)
        hPurity = getPurity(hPt, hPt_mc)
        hPurity.Write("hPt_ratio")

    output.Close()

