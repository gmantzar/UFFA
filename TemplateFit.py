import ROOT
import FemtoDreamSaver as FDS
import FemtoAnalysis as FA
import array as arr

class ftotal():
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
        self.axis = data.GetXaxis()
    def __call__(self, arr, par):
        total = 0
        nbin = self.axis.FindBin(arr[0])
        for n in range(len(self.histos)):
            total += par[n]*self.histos[n].GetBinContent(nbin)
        return total

def TemplateFit(fname, dca_data, dca_mcplots, dcacpa, dca_mcplots_names, fit_range, pt_bins, dirOut):
    # Output file

    if type(fname) == str:
        ofile = ROOT.TFile(dirOut + "TemplateFit_" + fname, "recreate")
    else:
        ofile = fname

    xAxis = dca_data.GetXaxis()
    yAxis = dca_data.GetYaxis()
    # setup for the bin ranges for the projections
    # output: [(bin1, bin2), (bin2, bin3), ...]
    if type(pt_bins) == int:                                            # input: bins
        pt_ent = pt_bins
        pt_range = []
        for n in range(pt_ent):
            pt_range.append((n + 1, n + 2))
        pt_rebin = None
    elif type(pt_bins) == list:
        if type(pt_bins[0]) == list and type(pt_bins[1]) == list:       # input: [[bin range], [rebin factors]]
            pt_ent = len(pt_bins[0]) - 1
            pt_range = []
            for i in range(pt_ent):
                bin_value1 = xAxis.FindBin(pt_bins[0][i])
                bin_value2 = xAxis.FindBin(pt_bins[0][i + 1])
                pt_range.append((bin_value1, bin_value2))
            if (len(pt_bins[1]) > 1):
                pt_rebin = pt_bins[1]*pt_ent
            else:
                pt_rebin = pt_bins[1]
        elif type(pt_bins[0]) == list and type(pt_bins[1]) == int:      # input: [[bin range], rebin factor]
            pt_ent = len(pt_bins[0]) - 1
            pt_range = []
            for i in range(pt_ent):
                bin_value1 = xAxis.FindBin(pt_bins[0][i])
                bin_value2 = xAxis.FindBin(pt_bins[0][i + 1])
                pt_range.append((bin_value1, bin_value2))
            pt_rebin = [pt_bins[1]]*pt_ent
        else:                                                           # input: [bin range]
            pt_ent = len(pt_bins) - 1
            pt_range = []
            for i in range(pt_ent):
                bin_value1 = xAxis.FindBin(pt_bins[i])
                bin_value2 = xAxis.FindBin(pt_bins[i + 1])
                pt_range.append((bin_value1, bin_value2))
            pt_rebin = None
    else:
        print("TemplateFit.py: pt_bins not an int or list of ranges: " + pt_bins)
        exit()

    print(pt_bins)
    print(pt_rebin)

    # titles of the graphs
    pt_names = []
    for n in range(pt_ent):
        name = "p_{T} range: [%.3f-%.3f) GeV" % (xAxis.GetBinLowEdge(pt_range[n][0]), xAxis.GetBinLowEdge(pt_range[n][1]))
        if pt_rebin:
            name += " rebin: " + str(pt_rebin[n])
        pt_names.append(name)

    dca_ent = len(dca_mcplots)
    if type(dcacpa) == str:
        if dcacpa.lower() == "dca":
            dcacpa = 'dca'
        elif dcacpa.lower() == "cpa":
            dcacpa = 'cpa'
        else:
            print("Error in plot type: \"dca\" or \"cpa\"")

    # dca cut
    CPAcut = 0.99
    if dcacpa == 'dca':
        fitmax = fit_range
        fitmin = -fitmax
    else:
        fitmin = 0.9
        fitmax = 1.0

    # graph initialization
    gChi = ROOT.TGraph(pt_ent - 1)
    m2ent = [0]*dca_ent
    dataEntries = []
    mcEntries = []
    parDCA_mc = []
    gDCA_mc = []
    for i in range(dca_ent):
        parDCA_mc.append([])
        gDCA_mc.append(ROOT.TGraph(pt_ent - 1))
        gDCA_mc[i].SetName(dca_mcplots_names[i])
        gDCA_mc[i].SetTitle(dca_mcplots_names[i])
        gDCA_mc[i].SetLineWidth(2)
        gDCA_mc[i].SetLineColor(i + 3)
        gDCA_mc[i].SetMarkerStyle(21)
        gDCA_mc[i].SetMarkerSize(1)
        gDCA_mc[i].SetMarkerColor(i + 3)
        gDCA_mc[i].GetXaxis().SetLabelSize(0.05)
        gDCA_mc[i].GetYaxis().SetLabelSize(0.05)
        gDCA_mc[i].GetXaxis().SetTitleSize(0.05)
        gDCA_mc[i].GetXaxis().SetTitle("<p_{T}> (GeV)")

    # main loop for fitting
    for n in range(pt_ent):
        # canvas for fitting
        canvas = ROOT.TCanvas("canvas_" + str(n + 1), "canvas_" + str(n + 1))
        canvas.SetCanvasSize(1024, 768)
        canvas.cd()
        ROOT.gPad.SetLogy()

        # data
        data = dca_data.ProjectionY("hDCAxy_" + str(n + 1), pt_range[n][0], pt_range[n][1] - 1)
        if pt_rebin and pt_rebin[n] != 1:
            data.Rebin(pt_rebin[n])
        data.SetAxisRange(fitmin, fitmax)
        if dcacpa == "cpa":
            data_int = data.Integral(data.FindBin(CPAcut), data.FindBin(1))
        else:
            data_int = data.Integral(data.FindBin(fitmin), data.FindBin(fitmax))
        dataEntries.append(data_int)
        if data.Integral():
            data.Scale(1. / data.Integral())

        # MC templates
        hDCA_mc = []
        for i in range(dca_ent):
            hDCA_mc.append(dca_mcplots[i].ProjectionY(dca_mcplots_names[i] + '_' + str(n + 1), pt_range[n][0], pt_range[n][1] - 1))
            if pt_rebin and pt_rebin[n] != 1:
                hDCA_mc[i].Rebin(pt_rebin[n])
            hDCA_mc[i].SetAxisRange(fitmin, fitmax)
            hDCA_mc[i].SetTitle(pt_names[n])
            if hDCA_mc[i].Integral():
                hDCA_mc[i].Scale(1. / hDCA_mc[i].Integral())

        # fit function
        adj = ftotal(data, hDCA_mc)
        ftot = ROOT.TF1("ftot", adj, fitmin, fitmax, dca_ent)

        for i in range(dca_ent):
            ftot.SetParameter(i, 1. / dca_ent)
            ftot.SetParLimits(i, 0., 1.)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", fitmin, fitmax)
        data.SetLineColor(1)
        data.SetLineWidth(2)
        data.Draw("hist")

        for i in range(dca_ent):
            hDCA_mc[i].Scale(ftot.GetParameter(i))
            hDCA_mc[i].SetLineColor(i + 3)
            hDCA_mc[i].SetLineWidth(2)
            hDCA_mc[i].Draw("same")

        htot = hDCA_mc[0].Clone("htot_" + str(n + 1))
        for i in range(1, dca_ent):
            htot.Add(hDCA_mc[i])
        htot.SetLineColor(2)
        htot.SetLineWidth(2)
        htot.Draw("same")

        # legend
        legend = ROOT.TLegend(0.15, 0.40, 0.35, 0.8)
        legend.SetTextSize(0.05)
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(1001)

        legend.AddEntry(data, "total", "l")
        legend.AddEntry(htot, "fit", "l")
        for i in range(dca_ent):
            legend.AddEntry(hDCA_mc[i], dca_mcplots_names[i], "l")
        legend.Draw()

        # write to file
        ofile.cd()
        for i in range(dca_ent):
            hDCA_mc[i].Write()
        canvas.Write()

        # fractions
        if dcacpa == "cpa":
            htot_int = htot.Integral(htot.FindBin(CPAcut), htot.FindBin(1))
        else:
            htot_int = htot.Integral(htot.FindBin(fitmin), htot.FindBin(fitmax))
        mcEntries.append(htot_int)
        if dcacpa == "cpa":
            for i in range(dca_ent):
                parDCA_mc[i].append(hDCA_mc[i].Integral(hDCA_mc[i].FindBin(CPAcut), hDCA_mc[i].FindBin(1)) / mcEntries[n])
        else:
            for i in range(dca_ent):
                tmp = hDCA_mc[i].Integral(hDCA_mc[i].FindBin(fitmin), hDCA_mc[i].FindBin(fitmax))
                if mcEntries[n]:
                    tmp = tmp / mcEntries[n]
                parDCA_mc[i].append(tmp)

        # fill x and y lists for graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1])) / 2
        for i in range(dca_ent):
            gDCA_mc[i].SetPoint(n, pt_avg, parDCA_mc[i][n])
        gChi.SetPoint(n, pt_avg, ftot.GetChisquare() / (ftot.GetNDF() + dca_ent))

        canvas.Close()
        del canvas

    # canvas for chi^2
    c2 = ROOT.TCanvas("c2", "c2", 0, 0, 1024, 768)
    c2.cd()
    gChi.SetLineColor(1)
    gChi.Draw()

    # calculate final pT weighted result
    m2tot = 0
    for i in range(pt_ent):
        m2tot += dataEntries[i]
        for j in range(dca_ent):
            m2ent[j] += dataEntries[i]*parDCA_mc[j][i]

    data_bins = dca_data.GetXaxis().GetNbins()
    BinRanges_pT = arr.array('d', [0]*(data_bins + 1))
    BinRanges_pT[data_bins] = dca_data.GetXaxis().GetBinUpEdge(data_bins)
    pT_Weights = ROOT.TH1F("pT_Weights", "pT_Weights", pt_ent, BinRanges_pT)
    for i in range(pt_ent - 1):
        pT_Weights.SetBinContent(i, dataEntries[i] / m2tot)

    # chi2 graph
    gChi.SetLineWidth(2)
    gChi.GetXaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetTitle("chi2/ndf")
    gChi.GetXaxis().SetTitle("<p_{T}> (GeV)")
    gChi.GetXaxis().SetTitleSize(0.05)
    gChi.GetYaxis().SetTitleSize(0.05)

    # canvas for all fractions
    fractions = ROOT.TCanvas("fractions", "fractions", 1024, 768)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    fractions.cd()
    fractions.SetFillColor(0)
    fractions.SetBorderMode(0)
    fractions.SetBorderSize(2)
    fractions.SetFrameBorderMode(0)
    fractions.SetMargin(0.15, 0.05, 0.2, 0.05)

    gEmpty = gDCA_mc[0]
    gEmpty.SetTitle("")
    gEmpty.SetFillStyle(1000)
    gEmpty.SetMinimum(0.)
    gEmpty.SetMaximum(1.0)
    gx = gEmpty.GetXaxis()
    gy = gEmpty.GetYaxis()
    gx.SetLimits(xAxis.GetBinLowEdge(pt_range[0][0] - 1), xAxis.GetBinLowEdge(pt_range[-1][1] + 1))
    gx.SetTitle("<p_{T}> (GeV)")
    gx.SetLabelFont(42)
    gx.SetLabelSize(0.06)
    gx.SetTitleSize(0.07)
    gx.SetTitleOffset(1)
    gx.SetTitleFont(42)
    gy.SetTitle("fractions")
    gy.SetLabelFont(42)
    gy.SetLabelSize(0.06)
    gy.SetTitleSize(0.07)
    gy.SetTitleFont(42)
    gEmpty.Draw("alp")
    for i in range(1, dca_ent):
        gDCA_mc[i].Draw("lp same")
    for i in range(pt_ent):
        line = ROOT.TLine(xAxis.GetBinLowEdge(pt_range[i][0]), parDCA_mc[0][i], xAxis.GetBinLowEdge(pt_range[i][1]), parDCA_mc[0][i])
        line.DrawClone("same")

    primAvg = ROOT.TF1("primAvg", "[0]", 0, 5)
    primAvg.SetParameter(0, m2ent[0] / m2tot)
    primAvg.SetLineColor(12)
    primAvg.SetLineStyle(2)
    primAvg.Draw("same")

    for i in range(dca_ent):
        print("%-12s %.6f" % (dca_mcplots_names[i], m2ent[i] / m2tot))

    ofile.cd()
    for i in range(dca_ent):
        gDCA_mc[i].Write()
    gChi.Write("chi2/ndf")
    fractions.Write("all")
    primAvg.Write()
    pT_Weights.Write()

    c2.Close()
    fractions.Close()
    del c2, fractions, htot, ftot, adj

