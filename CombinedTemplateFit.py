import ROOT
import FemtoAnalysis as FA
import FileReader as FR
import array as arr
import numpy as np
import ctypes

class ftotal:
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
        self.axis = data.GetXaxis()
        self.ent = len(mcdata)

    def __call__(self, arr, par):
        total = 0
        nbin = self.axis.FindBin(arr[0])
        for n in range(len(self.histos)):
            total += par[n]*self.histos[n].GetBinContent(nbin)

        return total

class global_chi2:
    def __init__(self, fits, ent):
        self.fits = fits
        self.ent = ent

    def __call__(self, par):
        par_arr = np.frombuffer(par, dtype=np.float64, count=self.ent)
        par_arr_arr = [par_arr[range(self.ent)] for n in range(len(self.fits))]

        total = sum([self.fits[n](par_arr_arr[n]) for n in range(len(self.fits))])
        return total

def find_bin_reduce_on_lower_edge(axis, value):
    found_bin = axis.FindBin(value)
    if value == axis.GetBinLowEdge(found_bin):
        found_bin -= 1
    return found_bin

# merge templates from given list and return list where the merged template is in the first provided index
def merge_templates(temps, newname, index1, index2):
    merged = temps[index1].Clone(newname)
    merged.Add(temps[index2])
    newlist = temps[0:index1] + [merged] + temps[index1 + 1:index2] + temps[index2 + 1:]
    return newlist

# prepare pt bin ranges for splitting pt bins
def get_pt_bin_ranges(axis, pt_bins):
    # output: [(bin1, bin2), (bin2, bin3), ...]
    pt_count = 0
    pt_range = []
    if type(pt_bins) == int:            # int that describes how many bins to individually fit
        pt_count = pt_bins
        for n in range(pt_count):
            pt_range.append((n + 1, n + 1))
    elif type(pt_bins) == list:         # list of pt edges
        pt_count = len(pt_bins) - 1
        for n in range(pt_count):
            bin_value1 = axis.FindBin(pt_bins[n])
            bin_value2 = find_bin_reduce_on_lower_edge(axis, pt_bins[n + 1])
            pt_range.append((bin_value1, bin_value2))
    return pt_range, pt_count

# setup function for pt
def setup_pt(axis, pt_bins):
    pt_count, pt_range = None, None
    if type(pt_bins[0]) == tuple:
        pt_count = len(pt_bins)
        pt_range = pt_bins
    else:
        pt_range, pt_count = get_pt_bin_ranges(axis, pt_bins)
    return pt_range, pt_count

# setup of the fitting ranges for the given pt bins
def setup_fitrange(fitrange, pt_count):
    if type(fitrange) != list:
        fitrange = [(-fitrange, fitrange)]*pt_count
    elif type(fitrange) == list:
        if type(fitrange[0]) not in [list, tuple]:
            if len(fitrange) < pt_count:
                fitrange += [fitrange[-1]]*(pt_count - len(fitrange))
            fitrange = [(-fitrange[n], fitrange[n]) for n in range(len(fitrange))]
    return fitrange

# obtain fractions from purity th1 plot
def get_fractions_from_purity(pt_bins, hpurity):
    pt_range, pt_count = setup_pt(hpurity.GetXaxis(), pt_bins)
    purity = []
    for npt in range(pt_count):
        pt_int = hpurity.Integral(pt_range[npt][0], pt_range[npt][1])
        purity.append(1 - (pt_int / (pt_range[npt][1] + 1 - pt_range[npt][0])))
    return purity

# returns th2 (dca-pt) distributions from a list of thn's to be used as input for the combined fitter
def get_input_templates_from_thn(name_temp, tdir, namelist, pt_bins, fitrange):
    dim_counter = range(2)
    temp_counter = range(len(namelist))
    labels = ['xy', 'z']

    # get templates
    temps_file = FR.FileReader(name_temp, tdir)
    temps = [temps_file.GetHisto(f"hDCAxy_{name}", "Tracks_MC") for name in namelist]

    # setup for the bin ranges for the projections
    xAxis = temps[0].GetAxis(0)
    pt_range, pt_count = setup_pt(xAxis, pt_bins)

    # setup of the fitting ranges for the given pt bins
    fitrange = setup_fitrange(fitrange, pt_count)

    # prepare empty histograms
    temps_cut = [[] for dim in dim_counter]
    for dim in dim_counter:
        for ntemp in temp_counter:
            temps_cut[dim].append(temps[dim].Clone(f"empty_{dim}_{ntemp}").Projection(dim + 1, 0).Clone(f"hDCA{labels[dim]}_{namelist[ntemp]}"))
            temps_cut[dim][ntemp].Reset()
            temps_cut[dim][ntemp].SetDirectory(0)

    # generate input th2 templates
    for dim in dim_counter:
        for ntemp in temp_counter:
            proj_list = ROOT.TList()
            for npt in range(pt_count):
                thn = temps[ntemp].Clone()
                thn.GetAxis(0).SetRange(pt_range[npt][0], pt_range[npt][1])
                ax1 = thn.GetAxis(1)
                ax1.SetRange(ax1.FindBin(fitrange[npt][0]), find_bin_reduce_on_lower_edge(ax1, fitrange[npt][1]))
                ax2 = thn.GetAxis(2)
                ax2.SetRange(ax2.FindBin(fitrange[npt][0]), find_bin_reduce_on_lower_edge(ax2, fitrange[npt][1]))
                proj = thn.Projection(dim + 1, 0)
                proj.SetName(f"proj_{labels[dim]}_{ntemp}_{npt}")
                proj_list.Add(proj)
            temps_cut[dim][ntemp].Merge(proj_list)

    return temps_cut, fitrange

# returns absolute fractions of DCAxy and DCAz
def get_rel_yield(name_temps, tdir, namelist, target, pt_bins, fitrange):
    temps, fitrange = get_input_templates_from_thn(name_temps, tdir, namelist, pt_bins, fitrange)
    temps = temps[0]
    xAxis = temps[0].GetXaxis()
    pt_range, pt_count = setup_pt(xAxis, pt_bins)

    temps_total = temps[0].Clone()
    temps_total.Reset()

    temp_counter = range(len(temps))
    for ntemp in temp_counter:
        temps[ntemp] = temps[ntemp].ProjectionX()

    temp_list = ROOT.TList()
    for ntemp in temp_counter:
        temp_list.Add(temps[ntemp].Clone())
    temps_total.Merge(temp_list)

    index = 0
    if type(target) == int:
        index = target
    elif type(target) == str:
        if target in namelist:
            index = namelist.index(target)

    temps_target = temps[index].Clone()
    temps_target.Divide(temps_total.Clone())

    fractions = []
    for npt in range(pt_count):
        pt_int = temps_target.Integral(pt_range[npt][0], pt_range[npt][1]) / (pt_range[npt][1] - pt_range[npt][0] + 1)
        fractions.append(pt_int)

    return fractions

def CombinedFit(fname, dir_out, dca_data, dca_templates, dca_names, fit_range, pt_bins, pt_rebin, temp_init, temp_limits, temp_fraction):
    # constants
    fit_count = len(dca_templates)
    temp_count = len(dca_templates[0])
    fit_counter = range(fit_count)
    temp_counter = range(temp_count)
    names_for_saving = [dca_data[n].GetName() for n in fit_counter]
    xAxis = dca_data[0].GetXaxis()
    yAxis = dca_data[0].GetYaxis()

    # renaming stupid names
    for n, name in enumerate(dca_names):
        if name == "SecondaryDaughterLambda":
            dca_names[n] = "SecLambda"
        if name == "SecondaryDaughterSigmaplus":
            dca_names[n] = "SecSigmaPlus"

    # color palette
    color_data = ROOT.kBlack
    color_fit  = ROOT.kRed
    color_temp = {
        "Primary":          ROOT.kGreen,
        "Secondary":        ROOT.kYellow,
        "WrongCollision":   ROOT.kOrange,
        "Material":         ROOT.kBlue,
        "Fake":             ROOT.kGray,
        "SecLambda":        ROOT.kMagenta,
        "SecSigmaPlus":     ROOT.kCyan,
    }

    # Output file
    ofile = ROOT.TFile(dirOut + "CombinedFit_" + fname, "recreate") if type(fname) == str else fname
    odirs = [ofile.mkdir(f"{names_for_saving[nfit]}_fits") for nfit in fit_counter]

    # initialize list of values if single rebin value was provided
    if pt_rebin and (len(pt_rebin) < len(pt_bins)):
        pt_rebin = pt_rebin*len(pt_bins)

    # setup for the bin ranges for the projections
    # output: [(bin1, bin2), (bin2, bin3), ...]
    pt_range, pt_count = setup_pt(xAxis, pt_bins)

    # setup for the fraction fitting parameters
    if type(temp_fraction) == dict:
        temp_fraction = [temp_fraction]

    # graph titles
    pt_names = []
    for n in range(pt_count):
        name = "p_{T} range: " + f"[{pt_bins[n]}-{pt_bins[n + 1]}) GeV"
        if pt_rebin:
            name += " rebin: " + str(pt_rebin[n])
        pt_names.append(name)

    # setup of the fitting ranges for the given pt bins
    if type(fit_range) != list:
        fit_range = [(-fit_range, fit_range)]*pt_count
    elif type(fit_range) == list:
        if type(fit_range[0]) not in [list, tuple]:
            if len(fit_range) < pt_count:
                fit_range += [fit_range[-1]]*(pt_count - len(fit_range))
            fit_range = [(-fit_range[n], fit_range[n]) for n in range(len(fit_range))]

    # graph initialization for fractions
    temp_graph = [[ROOT.TGraph(pt_count - 1) for ntemp in temp_counter] for nfit in fit_counter]
    for nfit in fit_counter:
        for ntemp in temp_counter:
            temp_graph[nfit][ntemp].SetName(f"{names_for_saving[nfit]}_{dca_names[ntemp]}")
            temp_graph[nfit][ntemp].SetTitle(dca_names[ntemp])
            temp_graph[nfit][ntemp].SetLineColor(ntemp + 3)
            temp_graph[nfit][ntemp].SetLineWidth(2)
            temp_graph[nfit][ntemp].SetMarkerColor(ntemp + 3)
            temp_graph[nfit][ntemp].SetMarkerSize(1)
            temp_graph[nfit][ntemp].SetMarkerStyle(21)
            temp_graph[nfit][ntemp].GetXaxis().SetLabelSize(0.05)
            temp_graph[nfit][ntemp].GetXaxis().SetTitleSize(0.05)
            temp_graph[nfit][ntemp].GetXaxis().SetTitle("<p_{T}> (GeV)")
            temp_graph[nfit][ntemp].GetYaxis().SetLabelSize(0.05)
            if dca_names[ntemp] in color_temp:
                temp_graph[nfit][ntemp].SetLineColor(color_temp[dca_names[ntemp]])
                temp_graph[nfit][ntemp].SetMarkerColor(color_temp[dca_names[ntemp]])

    ### fitting loop ###
    chi_graph = [ROOT.TGraph(pt_count) for nfit in fit_counter]
    data_ent = [[] for nfit in fit_counter]
    temp_ent = [[] for nfit in fit_counter]
    par_ent = [[[] for ntemp in temp_counter] for nfit in fit_counter]
    m2ent = [[0]*temp_count for nfit in fit_counter]
    total_chi2 = ROOT.TGraph(pt_count)
    for n in range(pt_count):
        # canvas for fitting
        canvas_arr = [ROOT.TCanvas(f"{names_for_saving[nfit]}_canvas_{n + 1}", f"{names_for_saving[nfit]}_canvas_{n + 1}") for nfit in fit_counter]
        for canvas in canvas_arr:
            canvas.SetCanvasSize(1024, 768)
        ROOT.gPad.SetLogy()

        # setup data
        data = []
        for nfit in fit_counter:
            data.append(dca_data[nfit].ProjectionY(f"{names_for_saving[nfit]}_{n + 1}", pt_range[n][0], pt_range[n][1]))
            if pt_rebin and pt_rebin[n] != 1:
                data[nfit].Rebin(pt_rebin[n])
            data[nfit].SetAxisRange(fit_range[n][0], fit_range[n][1])
            data_ent[nfit].append(data[nfit].Integral(data[nfit].FindBin(fit_range[n][0]), data[nfit].FindBin(fit_range[n][1])))
            if data[nfit].Integral():
                data[nfit].Scale(1. / data[nfit].Integral())

        # setup templates
        temp_hist = [[] for nfit in fit_counter]
        for nfit in fit_counter:
            for ntemp in temp_counter:
                temp_hist[nfit].append(dca_templates[nfit][ntemp].ProjectionY(f"{names_for_saving[nfit]}_{dca_names[ntemp]}_{n + 1}", pt_range[n][0], pt_range[n][1]))
                if pt_rebin and pt_rebin[n] != 1:
                    temp_hist[nfit][ntemp].Rebin(pt_rebin[n])
                temp_hist[nfit][ntemp].SetAxisRange(fit_range[n][0], fit_range[n][1])
                temp_hist[nfit][ntemp].SetTitle(pt_names[n])
                if temp_hist[nfit][ntemp].Integral():
                    temp_hist[nfit][ntemp].Scale(1. / temp_hist[nfit][ntemp].Integral())

        # setup data for wrappers
        data_opt = ROOT.Fit.DataOptions()
        data_range = ROOT.Fit.DataRange()
        data_range.SetRange(fit_range[n][0], fit_range[n][1])
        data_wrapped = [ROOT.Fit.BinData(data_opt, data_range) for nfit in fit_counter]
        for nfit in fit_counter:
            ROOT.Fit.FillData(data_wrapped[nfit], data[nfit])

        # fit functions
        fit_objs = [ftotal(data[nfit], temp_hist[nfit]) for nfit in fit_counter]
        fit_funcs = [ROOT.TF1(f"{names_for_saving[nfit]}_ftot", fit_objs[nfit], fit_range[n][0], fit_range[n][1], temp_count) for nfit in fit_counter]
        fit_wrapped = [ROOT.Math.WrappedMultiTF1(fit_funcs[nfit], temp_count) for nfit in fit_counter]
        fit_chi2 = [ROOT.Fit.Chi2Function(data_wrapped[nfit], fit_wrapped[nfit]) for nfit in fit_counter]
        chi2 = global_chi2(fit_chi2, temp_count)

        # fit parameter initialization
        fitter = ROOT.Fit.Fitter()
        fitter.Config().MinimizerOptions().SetPrintLevel(0)
        fitter.Config().SetMinimizer("Minuit2", "Migrad")
        fitter.Config().SetParamsSettings(temp_count, np.array([1/temp_count]*temp_count))                           # default init 1/temp_count
        if temp_init:
            fitter.Config().SetParamsSettings(temp_count, np.array(temp_init))                                 # initialize to temp_init
        for ntemp in temp_counter:
            fitter.Config().ParSettings(ntemp).SetLimits(0., 1.)                                                # default limits [0, 1]
            if temp_limits:
                fitter.Config().ParSettings(ntemp).SetLimits(temp_limits[ntemp][0], temp_limits[ntemp][1])
            if temp_fraction:
                for entry in temp_fraction:
                    if entry['temp_name'] == dca_names[ntemp]:
                        if entry['temp_init']:
                            fitter.Config().ParSettings(ntemp).SetValue(entry['temp_init'][n])
                            fitter.Config().ParSettings(ntemp).SetLimits(entry['temp_init'][n], entry['temp_init'][n])
                        if entry['temp_limits']:
                            if type(entry['temp_limits'][n]) == list:
                                fitter.Config().ParSettings(ntemp).SetLimits(entry['temp_limits'][n][0], entry['temp_limits'][n][1])
                            else:
                                if entry['temp_limits'][n] == 0 and not entry['temp_init']:
                                    if not entry['temp_init']:
                                        fitter.Config().ParSettings(ntemp).SetValue(0)
                                        fitter.Config().ParSettings(ntemp).Fix()

        # fitting
        chi2_functor = ROOT.Math.Functor(chi2, temp_count)
        fitter.FitFCN(chi2_functor, 0, sum([entry.GetNbinsX() for entry in data]), True)

        result = fitter.Result()
        result.Print(ROOT.std.cout)

        for nfit in fit_counter:
            fit_funcs[nfit].SetFitResult(result, np.array(temp_counter, dtype=np.int32))
            # rerun the fit function on the data with fixed parameters to calculate the chi2
            for ntemp in temp_counter:
                fit_funcs[nfit].FixParameter(ntemp, fit_funcs[nfit].GetParameter(ntemp))
            data[nfit].Fit(f"{names_for_saving[nfit]}_ftot", "S, N, R, M", "", fit_range[n][0], fit_range[n][1])

        # drawing
        for nfit in fit_counter:
            for ntemp in temp_counter:
                temp_hist[nfit][ntemp].Scale(fit_funcs[nfit].GetParameter(ntemp))
                temp_hist[nfit][ntemp].SetLineColor(ntemp + 3)
                temp_hist[nfit][ntemp].SetMarkerColor(ntemp + 3)
                if dca_names[ntemp] in color_temp:
                    temp_hist[nfit][ntemp].SetLineColor(color_temp[dca_names[ntemp]])
                    temp_hist[nfit][ntemp].SetMarkerColor(color_temp[dca_names[ntemp]])
                temp_hist[nfit][ntemp].SetLineWidth(2)

        htot = [temp_hist[nfit][0].Clone(f"{names_for_saving[nfit]}_htot_{n + 1}") for nfit in fit_counter]
        for nfit in fit_counter:
            for ntemp in range(1, temp_count):
                htot[nfit].Add(temp_hist[nfit][ntemp])
            htot[nfit].SetLineColor(2)
            htot[nfit].SetLineWidth(2)

        for nfit in fit_counter:
            canvas_arr[nfit].cd()
            ROOT.gPad.SetLogy()
            data[nfit].Draw("hist")
            htot[nfit].Draw("same")
            for ntemp in temp_counter:
                temp_hist[nfit][ntemp].Draw("same")

        # legend
        legend = ROOT.TLegend(0.15, 0.40, 0.35, 0.8)
        legend.SetTextSize(0.05)
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(1001)

        legend.AddEntry(data[0], "total", "l")
        legend.AddEntry(htot[0], "fit", "l")
        for ntemp in temp_counter:
            legend.AddEntry(temp_hist[0][ntemp], dca_names[ntemp], "l")
        for nfit in fit_counter:
            canvas_arr[nfit].cd()
            legend.Draw()

        # write to file
        ofile.cd()
        for nfit in fit_counter:
            newdir = odirs[nfit].mkdir(f"pt_{pt_bins[n]}-{pt_bins[n + 1]}")
            newdir.cd()
            for ntemp in temp_counter:
                temp_hist[nfit][ntemp].Write()
            canvas_arr[nfit].Write()
        ofile.cd()

        # fractions
        for nfit in fit_counter:
            temp_ent[nfit].append(htot[nfit].Integral(htot[nfit].FindBin(fit_range[n][0]), htot[nfit].FindBin(fit_range[n][1])))
            for ntemp in temp_counter:
                tmp = temp_hist[nfit][ntemp].Integral(temp_hist[nfit][ntemp].FindBin(fit_range[n][0]), temp_hist[nfit][ntemp].FindBin(fit_range[n][1]))
                if temp_ent[nfit][n]:
                    par_ent[nfit][ntemp].append(tmp / temp_ent[nfit][n])
                else:
                    par_ent[nfit][ntemp].append(tmp)

        # fill qa graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1])) / 2
        for nfit in fit_counter:
            chi_graph[nfit].SetPoint(n, pt_avg, fit_funcs[nfit].GetChisquare() / (fit_funcs[nfit].GetNDF()))
            for ntemp in temp_counter:
                temp_graph[nfit][ntemp].SetPoint(n, pt_avg, par_ent[nfit][ntemp][n])
            fit_funcs[nfit].SetFitResult(result, np.array(temp_counter, dtype=np.int32))
        total_chi2.SetPoint(n, pt_avg, fit_funcs[0].GetChisquare() / fit_funcs[0].GetNDF())

    # canvas for chi2
    chi_canvas = [ROOT.TCanvas(f"{names_for_saving[nfit]}_chi_canvas", f"{names_for_saving[nfit]}_chi_canvas", 1024, 768) for nfit in fit_counter]
    for nfit in fit_counter:
        chi_canvas[nfit].cd()
        chi_graph[nfit].SetLineColor(1)
        chi_graph[nfit].SetLineWidth(2)
        chi_graph[nfit].GetXaxis().SetLabelSize(0.05)
        chi_graph[nfit].GetYaxis().SetLabelSize(0.05)
        chi_graph[nfit].GetYaxis().SetTitle("chi2/ndf")
        chi_graph[nfit].GetXaxis().SetTitle("<p_{T}> (GeV)")
        chi_graph[nfit].GetXaxis().SetTitleSize(0.05)
        chi_graph[nfit].GetYaxis().SetTitleSize(0.05)

    # canvas for all fractions
    fractions_canvas = [ROOT.TCanvas(f"{names_for_saving[nfit]}_fractions", f"{names_for_saving[nfit]}_fractions", 1440, 1080) for nfit in fit_counter]
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    for nfit in fit_counter:
        fractions_canvas[nfit].SetFillColor(0)
        fractions_canvas[nfit].SetBorderMode(0)
        fractions_canvas[nfit].SetBorderSize(2)
        fractions_canvas[nfit].SetFrameBorderMode(0)
        fractions_canvas[nfit].SetMargin(0.15, 0.05, 0.2, 0.05)

    # Draw Graphs
    empty_graph = [temp_graph[nfit][0] for nfit in fit_counter]
    for nfit in fit_counter:
        fractions_canvas[nfit].cd()
        empty_graph[nfit].SetTitle("")
        empty_graph[nfit].SetFillStyle(1000)
        empty_graph[nfit].SetMinimum(0.)
        empty_graph[nfit].SetMaximum(1.0)
        gx = empty_graph[nfit].GetXaxis()
        gy = empty_graph[nfit].GetYaxis()
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
        empty_graph[nfit].Draw("alp")
        for ntemp in range(1, temp_count):
            temp_graph[nfit][ntemp].Draw("lp same")
        for ntemp in range(1):
            for n in range(pt_count):
                line = ROOT.TLine(xAxis.GetBinLowEdge(pt_range[n][0]), par_ent[nfit][ntemp][n], xAxis.GetBinLowEdge(pt_range[n][1]), par_ent[nfit][ntemp][n])
                line.DrawClone("same")

    # averages
    m2tot = [0]*fit_count
    for nfit in fit_counter:
        for n in range(pt_count):
            m2tot[nfit] += data_ent[nfit][n]
            for ntemp in temp_counter:
                m2ent[nfit][ntemp] += data_ent[nfit][n]*par_ent[nfit][ntemp][n]

    averages = [[ROOT.TF1(f"{dca_names[ntemp]}_avg", "[0]", 0, 5) for ntemp in temp_counter] for nfit in fit_counter]
    for nfit in fit_counter:
        for ntemp in temp_counter:
            averages[nfit][ntemp].SetParameter(0, m2ent[nfit][ntemp] / m2tot[nfit])
            averages[nfit][ntemp].SetLineColor(12)
            if dca_names[ntemp] in color_temp:
                averages[nfit][ntemp].SetLineColor(color_temp[dca_names[ntemp]])
            averages[nfit][ntemp].SetLineStyle(2)

    # print fit info
    for nfit in fit_counter:
        print(f"{names_for_saving[nfit]} fit:")
        for ntemp in temp_counter:
            print("  %-12s %.6f" % (dca_names[ntemp], m2ent[nfit][ntemp] / m2tot[nfit]))

    # write all plots
    ofile.cd()
    for nfit in fit_counter:
        odirs[nfit].cd()
        chi_graph[nfit].Write("chi2/ndf")
    ofile.cd()
    for ntemp in temp_counter:
        temp_graph[0][ntemp].Write()
    for ntemp in temp_counter:
        averages[0][ntemp].Write()
    fractions_canvas[nfit].Write()
    total_chi2.Write("total chi2/ndf")
    ofile.Close()

def TemplateFit(fname, dca_data, dca_templates, dcacpa, dca_names, fit_range, pt_bins, pt_rebin, dirOut, temp_init, temp_limits, temp_fraction):
    xAxis = dca_data.GetXaxis()
    yAxis = dca_data.GetYaxis()
    dca_ent = len(dca_templates)

    # Output file
    ofile = fname
    if type(fname) == str:
        ofile = ROOT.TFile(dirOut + "TemplateFit_" + fname, "recreate")

    # initialize list of values if single rebin value was provided
    if pt_rebin and (len(pt_rebin) < len(pt_bins)):
        pt_rebin = pt_rebin*len(pt_bins)

    # setup for the bin ranges for the projections
    # output: [(bin1, bin2), (bin2, bin3), ...]
    if type(pt_bins) == int:            # int that describes how many bins to individually fit
        pt_ent = pt_bins
        pt_range = []
        for n in range(pt_ent):
            pt_range.append((n + 1, n + 1))
    elif type(pt_bins) == list:         # list of pt edges
        pt_ent = len(pt_bins) - 1
        pt_range = []
        for i in range(pt_ent):
            bin_value1 = xAxis.FindBin(pt_bins[i])
            bin_value2 = find_bin_reduce_on_lower_edge(xAxis, pt_bins[i + 1])
            pt_range.append((bin_value1, bin_value2))
    else:
        print("TemplateFit.py: pt_bins not an int or list of ranges: ", pt_bins)
        exit()

    # setup for the fraction fitting parameters
    if type(temp_fraction) == dict:
        temp_fraction = [temp_fraction]

    # titles of the graphs
    pt_names = []
    for n in range(pt_ent):
        #name = "p_{T} range: [%.3f-%.3f) GeV" % (xAxis.GetBinLowEdge(pt_range[n][0]), xAxis.GetBinLowEdge(pt_range[n][1]))
        name = "p_{T} range: " + f"[{pt_bins[n]}-{pt_bins[n + 1]}) GeV"
        if pt_rebin:
            name += " rebin: " + str(pt_rebin[n])
        pt_names.append(name)

    # should be removed or expanded for actual cpa use
    if type(dcacpa) == str:
        if dcacpa.lower() == "dca":
            dcacpa = 'dca'
        elif dcacpa.lower() == "cpa":
            dcacpa = 'cpa'
        else:
            print("Error in plot type: \"dca\" or \"cpa\"")

    # dca cut
    CPAcut = 0.99
    fitmax = []
    fitmin = []
    if dcacpa == 'dca':
        if type(fit_range) == list:
            for value in fit_range:
                fitmax.append(value)
                fitmin.append(-value)
            if len(fitmax) < pt_ent:
                fitmax += [fit_range[-1]]*(pt_ent - len(fitmax))    # append the last fitrange to be used for the rest of the pt bins
                fitmin += [-fit_range[-1]]*(pt_ent - len(fitmax))
        else:
            fitmax = [fit_range]*pt_ent
            fitmin = [-fit_range]*pt_ent
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
    for n in range(dca_ent):
        parDCA_mc.append([])
        gDCA_mc.append(ROOT.TGraph(pt_ent - 1))
        gDCA_mc[n].SetName(dca_names[n])
        gDCA_mc[n].SetTitle(dca_names[n])
        gDCA_mc[n].SetLineWidth(2)
        gDCA_mc[n].SetLineColor(n + 3)
        gDCA_mc[n].SetMarkerStyle(21)
        gDCA_mc[n].SetMarkerSize(1)
        gDCA_mc[n].SetMarkerColor(n + 3)
        gDCA_mc[n].GetXaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetYaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitleSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitle("<p_{T}> (GeV)")

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
        data.SetAxisRange(fitmin[n], fitmax[n])
        if dcacpa == "cpa":
            data_int = data.Integral(data.FindBin(CPAcut), data.FindBin(1))
        else:
            data_int = data.Integral(data.FindBin(fitmin[n]), data.FindBin(fitmax[n]))
        dataEntries.append(data_int)
        if data.Integral():
            data.Scale(1. / data.Integral())

        # MC templates
        hDCA_mc = []
        for i in range(dca_ent):
            hDCA_mc.append(dca_templates[i].ProjectionY(dca_names[i] + '_' + str(n + 1), pt_range[n][0], pt_range[n][1] - 1))
            if pt_rebin and pt_rebin[n] != 1:
                hDCA_mc[i].Rebin(pt_rebin[n])
            hDCA_mc[i].SetAxisRange(fitmin[n], fitmax[n])
            hDCA_mc[i].SetTitle(pt_names[n])
            if hDCA_mc[i].Integral():
                hDCA_mc[i].Scale(1. / hDCA_mc[i].Integral())

        # fit function
        adj = ftotal(data, hDCA_mc)
        ftot = ROOT.TF1("ftot", adj, fitmin[n], fitmax[n], dca_ent)

        for i in range(dca_ent):
            ftot.SetParameter(i, 1. / dca_ent)
            ftot.SetParLimits(i, 0., 1.)

            if temp_init:
                ftot.SetParameter(i, temp_init[i])
            if temp_limits:
                ftot.SetParLimits(i, temp_limits[i][0], temp_limits[i][1])
            if temp_fraction:
                for entry in temp_fraction:
                    if entry['temp_name'] == dca_names[i]:
                        if entry['temp_init']:
                            ftot.SetParameter(i, entry['temp_init'][n])
                            ftot.SetParLimits(i, entry['temp_init'][n], entry['temp_init'][n])
                        if entry['temp_limits']:
                            if type(entry['temp_limits'][n]) == list:
                                ftot.SetParLimits(i, entry['temp_limits'][n][0], entry['temp_limits'][n][1])
                            else:
                                if entry['temp_limits'][n] == 0 and not entry['temp_init']:
                                    if not entry['temp_init']:
                                        ftot.FixParameter(i, 0)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", fitmin[n], fitmax[n])
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
            legend.AddEntry(hDCA_mc[i], dca_names[i], "l")
        legend.Draw()

        # write to file
        ofile.cd()
        newdir = ofile.mkdir(f"pt_{pt_bins[n]}-{pt_bins[n + 1]}")
        newdir.cd()
        for i in range(dca_ent):
            hDCA_mc[i].Write()
        canvas.Write()
        ofile.cd()

        # fractions
        if dcacpa == "cpa":
            htot_int = htot.Integral(htot.FindBin(CPAcut), htot.FindBin(1))
        else:
            htot_int = htot.Integral(htot.FindBin(fitmin[n]), htot.FindBin(fitmax[n]))
        mcEntries.append(htot_int)
        if dcacpa == "cpa":
            for i in range(dca_ent):
                parDCA_mc[i].append(hDCA_mc[i].Integral(hDCA_mc[i].FindBin(CPAcut), hDCA_mc[i].FindBin(1)) / mcEntries[n])
        else:
            for i in range(dca_ent):
                tmp = hDCA_mc[i].Integral(hDCA_mc[i].FindBin(fitmin[n]), hDCA_mc[i].FindBin(fitmax[n]))
                if mcEntries[n]:
                    parDCA_mc[i].append(tmp / mcEntries[n])
                else:
                    parDCA_mc[i].append(tmp)

        # fill x and y lists for graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1])) / 2
        for i in range(dca_ent):
            gDCA_mc[i].SetPoint(n, pt_avg, parDCA_mc[i][n])
        gChi.SetPoint(n, pt_avg, ftot.GetChisquare() / (ftot.GetNDF()))

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
        print("%-12s %.6f" % (dca_names[i], m2ent[i] / m2tot))

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

