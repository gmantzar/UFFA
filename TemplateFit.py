import ROOT
import FemtoDreamReader as FDR
import FemtoAnalysis as FA
import array as arr

class ftotal:
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
    def __call__(self, arr, par):
        total = 0
        nbin = self.data.GetXaxis().FindBin(arr[0])
        for n in range(len(self.histos)):
            total += par[n]*self.histos[n].GetBinContent(nbin)
        return total

def TemplateFit(fname, dca_data, dca_mcplots, dcacpa, dca_mcplots_names, pt_bins, dirOut = None):
    # Output file
    if dirOut and dirOut != "" and dirOut[-1] != '/':
        dirOut += '/'
    ofile = ROOT.TFile(dirOut + "TemplateFit_" + fname, "recreate")

    dca_ent = len(dca_mcplots)
    if type(dcacpa) == str:
        if dcacpa.lower() == "dca":
            dcacpa = 1
        elif dcacpa.lower() == "cpa":
            dcacpa = 2
        else:
            print("Error in plot type: \"dca\" or \"cpa\"")

    # dca cut
    if dcacpa == 1:
        DCAcut = 0.1
        fitmin = -0.3
        fitmax = 0.3
    else:
        DCAcut = 0.99
        fitmin = 0.9
        fitmax = 1.0

    gChi = ROOT.TGraph(pt_bins - 1)
    m2ent = [0]*dca_ent
    dataEntries = []
    mcEntries = []
    parDCA_mc = []
    gDCA_mc = []
    for i in range(dca_ent):
        parDCA_mc.append([])
        gDCA_mc.append(ROOT.TGraph(pt_bins - 1))
        gDCA_mc[i].SetName(dca_mcplots_names[i])
        gDCA_mc[i].SetLineWidth(2)
        gDCA_mc[i].GetXaxis().SetLabelSize(0.05)
        gDCA_mc[i].GetYaxis().SetLabelSize(0.05)
        gDCA_mc[i].GetXaxis().SetTitleSize(0.05)
        gDCA_mc[i].GetXaxis().SetTitle("p (GeV)")
        gDCA_mc[i].SetLineColor(i + 3)

    for n in range(1, pt_bins):
        # canvas for fitting
        canvas = ROOT.TCanvas("canvas_" + str(n), "canvas_" + str(n), 0, 0, 650, 550)
        canvas.cd()
        ROOT.gPad.SetLogy()

        # data
        data = dca_data.ProjectionY("DCAxy_" + str(n), n, n)
        data.Scale(1./data.Integral())
        dataEntries.append(data.Integral(data.FindBin(DCAcut), data.FindBin(1)))

        # MC templates
        hDCA_mc = []
        for i in range(dca_ent):
            hDCA_mc.append(dca_mcplots[i].ProjectionY(dca_mcplots_names[i] + '_' + str(n), n, n))
            hDCA_mc[i].Scale(1./hDCA_mc[i].Integral())

        # fitting
        adj = ftotal(data, hDCA_mc)
        ftot = ROOT.TF1("ftot", adj, fitmin, fitmax, dca_ent)

        for i in range(dca_ent):
            ftot.SetParameter(i, 0.5)
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

        htot = hDCA_mc[0].Clone("htot_" + str(n))
        for i in range(1, dca_ent):
            htot.Add(hDCA_mc[i])
        htot.SetLineColor(2)
        htot.SetLineWidth(2)
        htot.Draw("same")

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

        # ???
        mcEntries.append(htot.Integral(htot.FindBin(DCAcut), htot.FindBin(1)))
        for i in range(dca_ent):
            parDCA_mc[i].append(hDCA_mc[i].Integral(hDCA_mc[i].FindBin(DCAcut), hDCA_mc[i].FindBin(1)) / mcEntries[n - 1])

        # fill x and y lists for graphs
        for i in range(dca_ent):
            gDCA_mc[i].SetPoint(n - 1, dca_data.GetXaxis().GetBinCenter(n), parDCA_mc[i][n - 1])
        gChi.SetPoint(n - 1, dca_data.GetXaxis().GetBinCenter(n), ftot.GetChisquare() / (ftot.GetNDF() + 1))

    # canvas for chi^2
    c2 = ROOT.TCanvas("c2", "c2", 0, 0, 650, 550)
    c2.cd()
    gChi.SetLineColor(1)
    gChi.Draw()

    # calculate final pT weighted result
    m2tot = 0
    for n in range(pt_bins - 1):
        m2tot += dataEntries[n]
        for i in range(dca_ent):
            m2ent[i] += dataEntries[n]*parDCA_mc[i][n]

    data_bins = dca_data.GetXaxis().GetNbins()
    BinRanges_pT = arr.array('d', [0]*(data_bins + 1))
    BinRanges_pT[data_bins] = dca_data.GetXaxis().GetBinUpEdge(data_bins)
    pT_Weights = ROOT.TH1F("pT_Weights", "pT_Weights", pt_bins, BinRanges_pT)
    for n in range(pt_bins - 1):
        pT_Weights.SetBinContent(n + 1, dataEntries[n] / m2tot)

    # what even is this
    d12 = ROOT.TCanvas("d12", "d12", 0, 0, 650, 550)
    d12.cd()
    gChi.SetLineWidth(2)
    gChi.GetXaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetTitle("chi2/ndf")
    gChi.GetXaxis().SetTitle("p (GeV)")
    gChi.GetXaxis().SetTitleSize(0.05)
    gChi.GetYaxis().SetTitleSize(0.05)

    gDCA_mc[0].GetYaxis().SetRangeUser(0.,1.)
    gDCA_mc[0].Draw("")
    for i in range(1, dca_ent):
        gDCA_mc[i].Draw("same")

    primAvg = ROOT.TF1("primAvg", "[0]", 0, 5)
    primAvg.SetParameter(0, m2ent[0] / m2tot)
    primAvg.SetLineColor(12)
    primAvg.SetLineStyle(2)
    primAvg.Draw("same")

    ofile.cd()
    for i in range(dca_ent):
        gDCA_mc[i].Write()
    gChi.Write("chi2/ndf")
    d12.Write("all")
    pT_Weights.Write()


