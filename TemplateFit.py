import ROOT
import FemtoDreamReader as FDR
import FemtoAnalysis as FA

class ftotal:
    def __call__(self, arr, par):
        nbin = data.GetXaxis().FindBin(arr[0])
        prim_adj = par[0]*prim.GetBinContent(nbin)
        lam_adj = par[1]*prim.GetBinContent(nbin)
        sig_adj = par[2]*prim.GetBinContent(nbin)
        mat_adj = par[3]*prim.GetBinContent(nbin)
        fake_adj = par[4]*prim.GetBinContent(nbin)
        return prim_adj + lam_adj + sig_adj + mat_adj + fake_adj

# def TemplateFit(dirIn, fname, nfile, htype, bins, dirOut = None):
    #conf = FA.config(dirIn, dirOut, fname, False, nfile, False, htype, False, bins, False)
    #fdr = FDR.FemtoDreamReader(conf[0] + conf[1], conf[2])

def TemplateFit(filename):
    # Output file
    ofile = ROOT.TFile("TemplateFit_"+filename, "recreate")
    fdr = FDR.FemtoDreamReader(filename, "")

    # get data
    DCAxy_data = fdr.get_dca()

    # get MC templates
    DCAxy_prim  = fdr.get_histo("Tracks_oneMC/hDCAxy_Primary")
    DCAxy_lam   = fdr.get_histo("Tracks_oneMC/hDCAxy_DaughterLambda")
    DCAxy_sig   = fdr.get_histo("Tracks_oneMC/hDCAxy_DaughterSigmaplus")
    DCAxy_mat   = fdr.get_histo("Tracks_oneMC/hDCAxy_Material")
    DCAxy_fake  = fdr.get_histo("Tracks_oneMC/hDCAxy_Fake")

    # fitting function
    adj = ftotal()
    ftot = ROOT.TF1("ftot", adj, 0.9, 1, 4)

    # number of pt bins, might later be a function parameter
    Num_pT_bins = 8

    # x and y lists for graphs
    x_prim, y_prim = [], []
    x_lam,  y_lam  = [], []
    x_sig,  y_sig  = [], []
    x_mat,  y_mat  = [], []
    x_fake, y_fake = [], []
    x_chi,  y_chi  = [], []

    # dca cut
    DCAcut = 0.1
    CPAcut = 0.1

    dataEntries = []
    mcEntries = []
    parPrim = []
    parLam = []
    parSig = []
    parMat = []
    parFake = []

    # TGraphs
    gPrim = ROOT.TGraph(Num_pT_bins - 1)
    gLam = ROOT.TGraph(Num_pT_bins - 1)
    gSig = ROOT.TGraph(Num_pT_bins - 1)
    gMat = ROOT.TGraph(Num_pT_bins - 1)
    gFake = ROOT.TGraph(Num_pT_bins - 1)
    gChi = ROOT.TGraph(Num_pT_bins - 1)

    gPrim.SetName("Primary")
    gLam.SetName("Lambda")
    gSig.SetName("Sigma+")
    gMat.SetName("Material")
    gFake.SetName("Fake")
    gChi.SetName("Chi")

    for n in range(1, Num_pT_bins):
        # canvas for fitting
        c = ROOT.TCanvas("Can1_" + str(n), "Can1_" + str(n), 0, 0, 650, 550)
        c.cd()
        ROOT.gPad.SetLogy()

        # data
        data = DCAxy_data.ProjectionY("DCAxy_" + str(n), n, n)
        dataEntries.append(data.Integral(data.FindBin(DCAcut), data.FindBin(1)))

        momLow  = DCAxy_data.GetXaxis().GetBinLowEdge(n)
        momHigh = DCAxy_data.GetXaxis().GetBinLowEdge(n) + DCAxy_data.GetXaxis().GetBinWidth(n)

        data.Scale(1./data.Integral())

        # MC templates
        prim = DCAxy_prim.ProjectionY("prim_" + str(n), n, n)
        lam  = DCAxy_lam.ProjectionY("lam_" + str(n), n, n)
        sig  = DCAxy_sig.ProjectionY("sig_" + str(n), n, n)
        mat  = DCAxy_mat.ProjectionY("mat_" + str(n), n, n)
        fake = DCAxy_fake.ProjectionY("fake_" + str(n), n, n)

        # normalize
        try:
            prim.Scale(1./prim.Integral())
            lam.Scale(1./lam.Integral())
            sig.Scale(1./sig.Integral())
            mat.Scale(1./mat.Integral())
            fake.Scale(1./fake.Integral())
        except:
            pass

        # fitting
        ftot.SetParameters(0.6, 0.15, 0.1, 0.05, 0.1)
        ftot.SetParLimits(0, 0., 1.)
        ftot.SetParLimits(1, 0., 1.)
        ftot.SetParLimits(2, 0., 1.)
        ftot.SetParLimits(3, 0., 1.)
        ftot.SetParLimits(4, 0., 1.)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", 0.9, 1)
        data.SetLineColor(1)
        data.SetLineWidth(2)
        data.Draw("hist")

        prim.Scale(ftot.GetParameter(0))
        prim.SetLineColor(2)
        prim.SetLineWidth(2)
        prim.Draw("same")

        lam.Scale(ftot.GetParameter(1))
        lam.SetLineColor(4)
        lam.SetLineWidth(2)
        lam.Draw("same")

        sig.Scale(ftot.GetParameter(2))
        sig.SetLineColor(6)
        sig.SetLineWidth(2)
        sig.Draw("same")

        mat.Scale(ftot.GetParameter(3))
        mat.SetLineColor(3)
        mat.SetLineWidth(2)
        mat.Draw("same")

        fake.Scale(ftot.GetParameter(4))
        fake.SetLineColor(95)
        fake.SetLineWidth(2)
        fake.Draw("same")

        htot = prim.Clone("htot_" + str(n))
        htot.Add(lam)
        htot.Add(sig)
        htot.Add(mat)
        htot.Add(fake)
        htot.SetLineColor(2)
        htot.SetLineWidth(2)
        htot.Draw("same")

        # write to file
        ofile.cd()
        prim.Write()
        lam.Write()
        sig.Write()
        mat.Write()
        fake.Write()

        # ???
        mcEntries.append(htot.Integral(htot.FindBin(DCAcut), htot.FindBin(1)))
        if mcEntries[n - 1] == 0:
            continue
        parPrim.append(prim.Integral(prim.FindBin(DCAcut), prim.FindBin(1)) / mcEntries[n - 1])
        parLam.append(lam.Integral(lam.FindBin(DCAcut), lam.FindBin(1)) / mcEntries[n - 1])
        parSig.append(sig.Integral(sig.FindBin(DCAcut), sig.FindBin(1)) / mcEntries[n - 1])
        parMat.append(mat.Integral(mat.FindBin(DCAcut), mat.FindBin(1)) / mcEntries[n - 1])
        parFake.append(fake.Integral(fake.FindBin(DCAcut), fake.FindBin(1)) / mcEntries[n - 1])

        # fill x and y lists for graphs
        gPrim.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), parPrim[n - 1])
        gLame.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), parLam[n - 1])
        gSig.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), parSig[n - 1])
        gMat.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), parMat[n - 1])
        gFake.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), parFake[n - 1])
        gChi.SetPoint(n - 1, DCAxy_data.GetXaxis().GetBinCenter(n), ftot.GetChisquare() / ftot.GetNDF())

    # im not sure whats happening here but i might know when im finished
    c2 = ROOT.TCanvas("c2", "c2", 0, 0, 650, 550)
    c2.cd()
    gChi.SetLineColor(1)
    gChi.Draw()

    # calculate final pT weighted result
    m2tot = m2prim = m2lam = m2sig = m2mat = m2fake = 0
    for n in range(Num_pT_bins - 1):
        try:
            m2tot += dataEntries[n]
            m2prim += dataEntries[n]*parPrim[n]
            m2lam += dataEntries[n]*parLam[n]
            m2sig += dataEntries[n]*parSig[n]
            m2mat += dataEntries[n]*parMat[n]
            m2fake += dataEntries[n]*parFake[n]
        except:
            pass

    try:
        data_bins = DCAxy_data.GetXaxis().GetNbins()
        BinRanges_pT = [0]*(data_bins + 1)
        BinRanges_pT[data_bins] = DCAxy_data.GetXaxis().GetBinUpEdge(data_bins)
        pT_Weights = ROOT.TH1F("pT_Weights", "pT_Weights", Num_pT_bins, BinRanges_pT)
        for n in range(Num_pT_bins - 1):
            pT_Weigths.SetBinContent(n + 1, dataEntries[n] / m2tot)
    except:
        pass

    # what even is this
    d12 = ROOT.TCanvas("d12", "d12", 0, 0, 650, 550)
    d12.cd()
    gPrim.SetLineWidth(2)
    gLam.SetLineWidth(2)
    gSig.SetLineWidth(2)
    gMat.SetLineWidth(2)
    gFake.SetLineWidth(2)
    gChi.SetLineWidth(2)
    gPrim.GetXaxis().SetLabelSize(0.05)
    gPrim.GetYaxis().SetLabelSize(0.05)
    gPrim.GetXaxis().SetTitleSize(0.05)
    gPrim.GetXaxis().SetTitle("p (GeV)")
    gChi.GetXaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetTitle("chi2/ndf")
    gChi.GetXaxis().SetTitle("p (GeV)")
    gChi.GetXaxis().SetTitleSize(0.05)
    gChi.GetYaxis().SetTitleSize(0.05)

    gPrim.SetLineColor(2)
    gPrim.GetYaxis().SetRangeUser(0.,1.)
    gPrim.Draw("")
    gLam.SetLineColor(4)
    gLam.Draw("same")
    gSig.SetLineColor(6)
    gSig.Draw("same")
    gMat.SetLineColor(3)
    gMat.Draw("same")
    gFake.SetLineColor(95)
    gFake.Draw("same")

    try:
        primAvg = ROOT.TF1("primAvg", "[0]", 0, 5)
        primAvg.SetParameter(0, m2prim / m2tot)
        primAvg.SetLineColor(1)
        primAvg.SetLineStyle(2)
        primAvg.Draw("same")
    except:
        pass

    ofile.cd()
    gPrim.Write()
    gLam.Write()
    gSig.Write()
    gMat.Write()
    gFake.Write()
    gChi.Write()
    try:
        pT_Weights.Write()
    except:
        pass


