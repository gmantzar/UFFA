import ROOT

class ftotal:
    def __call__(self, arr, par):
        nbin = data.GetXaxis().FindBin(arr[0])
        prim_adj = par[0]*prim.GetBinContent(nbin)
        lam_adj = par[1]*prim.GetBinContent(nbin)
        sig_adj = par[2]*prim.GetBinContent(nbin)
        mat_adj = par[3]*prim.GetBinContent(nbin)
        fake_adj = par[4]*prim.GetBinContent(nbin)
        return prim_adj + lam_adj + sig_adj + mat_adj + fake_adj

def TemplateFit(filename):
    # TFile and TDirectory
    ifile = ROOT.TFile(filename, "read")
    std_dir = ifile.GetDirectory("_std")
    mc_dir = ifile.GetDirectory("_std/MC")

    # check if MC data exists in file
    if not mc_dir:
        print("Error: TDirectory \"MC\" does not exist!")
        exit()

    # Output file
    ofile = ROOT.TFile("TemplateFit_"+"filename", "recreate")

    # get data
    DCAxy_data = std_dir.Get("hDCAxy")

    # get MC templates
    DCAxy_prim  = mc_dir.Get("hDCAxy_Primary")
    DCAxy_lam   = mc_dir.Get("hDCAxy_DaughterLambda")
    DCAxy_sig   = mc_dir.Get("hDCAxy_DaughterSigmaplus")
    DCAxy_mat   = mc_dir.Get("hDCAxy_Material")
    DCAxy_fake  = mc_dir.Get("hDCAxy_Fake")

    # fitting function
    adj = ftotal()
    ftot = ROOT.TF1("ftot", adj, 0.9, 1, 4)

    # number of pt bins, might later be a function parameter
    Num_pT_bins = 8

    # x and y lists for graphs
    x_prim, y_prim = array('d'), array('d')
    x_lam, y_lam = array('d'), array('d')
    x_sig, y_sig = array('d'), array('d')
    x_mat, y_mat = array('d'), array('d')
    x_fake, y_fake = array('d'), array('d')
    x_chi, y_chi = array('d'), array('d')

    # dca cut
    DCA_cut = 0.1

    # ???
    dataEntries = []
    mcEntries = []
    parPrim = []
    parLam = []
    parSig = []
    parMat = []
    parFake = []

    for n in range(1, Num_pT_bins):
        # canvas for fitting
        c = ROOT.TCanvas()
        c.cd()
        gPad.SetLogy()

        # data
        data = DCAxy_data.ProjectionY("DCAxy_"+n, n, n)
        data.Scale(1./data.Integral())

        # MC templates
        prim   = DCAxy_prim.ProjectionY("prim_"+n, n, n)
        lam    = DCAxy_lam.ProjectionY("lam_"+n, n, n)
        sig    = DCAxy_sig.ProjectionY("sig_"+n, n, n)
        mat    = DCAxy_mat.ProjectionY("mat_"+n, n, n)
        fake   = DCAxy_fake.ProjectionY("fake_"+n, n, n)

        # normalize
        prim.Scale(1./prim.Integral())
        lam.Scale(1./lam.Integral())
        sig.Scale(1./sig.Integral())
        mat.Scale(1./mat.Integral())
        fake.Scale(1./fake.Integral())

        # fitting
        ftot.SetParameters(0.6, 0.15, 0.1, 0.05, 0.1)
        ftot.SetParLimits(0, 0., 1.)
        ftot.SetParLimits(1, 0., 1.)
        ftot.SetParLimits(2, 0., 1.)
        ftot.SetParLimits(3, 0., 1.)
        ftot.SetParLimits(4, 0., 1.)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", 0.9, 1)
        data.SetLineColor(kBlack)
        data.SetLineWidth(2)
        data.Draw("hist")

        prim.Scale(ftot.GetParameter(0))
        prim.SetLineColor(kRed + 2)
        prim.SetLineWidth(2)
        prim.Draw("same")

        lam.Scale(ftot.GetParameter(1))
        lam.SetLineColor(kBlue + 1)
        lam.SetLineWidth(2)
        lam.Draw("same")

        sig.Scale(ftot.GetParameter(2))
        sig.SetLineColor(kMagenta + 3)
        sig.SetLineWidth(2)
        sig.Draw("same")

        mat.Scale(ftot.GetParameter(3))
        mat.SetLineColor(kGreen + 3)
        mat.SetLineWidth(2)
        mat.Draw("same")

        fake.Scale(ftot.GetParameter(4))
        fake.SetLineColor(kOrange + 1)
        fake.SetLineWidth(2)
        fake.Draw("same")

        htot = prim.Clone("htot_"+n)
        htot.Add(lam)
        htot.Add(sig)
        htot.Add(mat)
        htot.Add(fake)
        htot.SetLineColor(kRed)
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
        dataEntries.append(data.Integral(data.FindBin(DCA_cut), data.FindBin(1)))
        mcEntries.append(htot.Integral(htot.FindBind(CPAcut), htot.FindBin(1)/mcEntries[n - 1]))
        parPrim.append(prim.Integral(prim.FindBind(CPAcut), prim.FindBin(1)/mcEntries[n - 1]))
        parLam.append(lam.Integral(lam.FindBind(CPAcut), lam.FindBin(1)/mcEntries[n - 1]))
        parSig.append(sig.Integral(sig.FindBind(CPAcut), sig.FindBin(1)/mcEntries[n - 1]))
        parMat.append(mat.Integral(mat.FindBind(CPAcut), mat.FindBin(1)/mcEntries[n - 1]))
        parFake.append(fake.Integral(fake.FindBind(CPAcut), fake.FindBin(1)/mcEntries[n - 1]))

        # fill x and y lists for graphs
        x_prim.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        x_lam.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        x_sig.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        x_mat.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        x_fake.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        x_chi.append(DCAxy_data.GetXaxis().GetBinCenter(n))
        y_prim.append(parPrim[n - 1])
        y_lam.append(parLam[n - 1])
        y_sig.append(parSig[n - 1])
        y_mat.append(parMat[n - 1])
        y_fake.append(parFake[n - 1])
        y_chi.append(ftot.GetChisquare() / ftot.GetNDF()))

    gPrim = TGraph(Num_pT_bin - 1, x_prim, y_prim)
    gLam = TGraph(Num_pT_bin - 1, x_lam, y_lam)
    gSig = TGraph(Num_pT_bin - 1, x_sig, y_sig)
    gMat = TGraph(Num_pT_bin - 1, x_mat, y_mat)
    gFake = TGraph(Num_pT_bin - 1, x_fake, y_fake)
    gChi = TGraph(Num_pT_bin - 1, x_chi, y_chi)

	gPrim->SetName("Primary");
	gLam->SetName("Lambda");
	gSig->SetName("Sigma+");
	gMat->SetName("Material");
	gFake->SetName("Fake");
	gChi->SetName("Chi");

    # im not sure whats happening here but i might know when im finished
    c2 = ROOT.TCanvas("c2", "c2", 0, 0, 650, 550)
    c2.cd()
    gChi.SetLineColor(1)
    gChi.Draw()

    # calculate final pT weighted result
    m2tot = m2prim = m2lam = m2sig = m2mat = m2fake = 0
    for n in range(Num_pT_bins - 1):
        m2tot += dataEntries[n]
        m2prim += dataEntries[n]*parPrim[n]
        m2lam += dataEntries[n]*parLam[n]
        m2sig += dataEntries[n]*parSig[n]
        m2mat += dataEntries[n]*parMat[n]
        m2fake += dataEntries[n]*parFake[n]

    data_bins = DCAxy_data.GetXaxis().GetNbins()
    BinRanges_pT = [0]*(data_bins + 1)
    BinRanges_pT[data_bins] = DCAxy_data.GetXaxis().GetBinUpEdge(data_bins)
    pT_Weights = ROOT.TH1F("pT_Weights", "pT_Weights", Num_pT_bins, BinRanges_pT)
    for n in range(Num_pT_bins - 1):
        pT_Weigths.SetBinContent(n + 1, dataEntries[n] / m2tot)

    # what even is this
    d12 = ROOT.TCanvas("d12", "d12", 0, 0, 650, 550)
    d12->cd()
    gPrim->SetLineWidth(2)
    gLam->SetLineWidth(2)
    gSig->SetLineWidth(2)
    gMat->SetLineWidth(2)
    gFake->SetLineWidth(2)
    gChi->SetLineWidth(2)
    gPrim->GetXaxis()->SetLabelSize(0.05)
    gPrim->GetYaxis()->SetLabelSize(0.05)
    gPrim->GetXaxis()->SetTitleSize(0.05)
    gPrim->GetXaxis()->SetTitle("p (GeV)")
    gChi->GetXaxis()->SetLabelSize(0.05)
    gChi->GetYaxis()->SetLabelSize(0.05)
    gChi->GetYaxis()->SetTitle("chi2/ndf")
    gChi->GetXaxis()->SetTitle("p (GeV)")
    gChi->GetXaxis()->SetTitleSize(0.05)
    gChi->GetYaxis()->SetTitleSize(0.05)

    gPrim->SetLineColor(kRed + 2)
    gPrim->GetYaxis()->SetRangeUser(0.,1.)
    gPrim->Draw("")
    gLam->SetLineColor(kBlue + 1)
    gLam->Draw("same")
    gSig->SetLineColor(kMagenta + 3)
    gSig->Draw("same")
    gMat->SetLineColor(kGreen + 3)
    gMat->Draw("same")
    gFake->SetLineColor(kOrange + 1)
    gFake->Draw("same")

    primAvg = TF1("primAvg", "[0]", 0, 5)
    primAvg->SetParameter(0, m2prim / m2tot)
    primAvg->SetLineColor(kBlack)
    primAvg->SetLineStyle(2)
    primAvg->Draw("same")

    ofile.cd()
    gPrim.Write()
    gLam.Write()
    gSig.Write()
    gMat.Write()
    gFake.Write()
    gChi.Write()
    pT_Weights.Write()


