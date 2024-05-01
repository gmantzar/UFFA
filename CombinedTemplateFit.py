import ROOT
import FemtoAnalysis as FA
import FileReader as FR
import array as arr
import numpy as np
import ctypes

# fitting object to be passed to tf1
class ftotal:
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
        self.axis = data.GetXaxis()
        self.ent = len(mcdata)

    def __call__(self, arr, par):
        nbin = self.axis.FindBin(arr[0])
        total = 0
        for ntemp in range(len(self.histos)):
            total += par[ntemp]*self.histos[ntemp].GetBinContent(nbin)

        return total

# fitting object to be passed to tf1
class ftotal_chi2:
    def __init__(self, data, mcdata):
        self.data = data
        self.histos = mcdata
        self.axis = data.GetXaxis()
        self.ent = len(mcdata)

    def __call__(self, arr, par):
        nbin = self.axis.FindBin(arr[0])
        total = self.data.GetbinContent(nbin)
        error = self.data.GetBinError(nbin)**2
        for ntemp in range(len(self.histos)):
            total -= par[ntemp]*self.histos[ntemp].GetBinContent(nbin)
            error += self.histos[ntemp].GetBinError(nbin)**2

        return total

# global chi2 fitting object for combined fit
class global_chi2:
    def __init__(self, fits, ent):
        self.fits = fits
        self.ent = ent

    def __call__(self, par):
        par_arr = np.frombuffer(par, dtype=np.float64, count=self.ent)
        par_arr_arr = [par_arr[range(self.ent)] for n in range(len(self.fits))]

        total = sum([self.fits[n](par_arr_arr[n]) for n in range(len(self.fits))])
        return total

# finds the given bin on the axis and if its on the edge of the next bin, reduces the bin number
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

# normalize, reweight and merge histograms
def merge_templates_weighted(newname, temp1, temp2, weight1, weight2):
    hist1 = temp1.Clone(newname)
    hist2 = temp2.Clone()
    hist1.Scale(1/hist1.Integral())
    hist2.Scale(1/hist2.Integral())
    if weight1:
        hist1.Scale(weight1)
    if weight2:
        hist2.Scale(weight2)
    hist1.Add(hist2)
    return hist1

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

# colors for fit contributions
def setup_colors():
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
        "Flat":             ROOT.kBlue,
    }
    return color_data, color_fit, color_temp

# setup rebin factors
def setup_rebin(pt_rebin, pt_count):
    if type(pt_rebin) in [int, float]:
        pt_rebin = [pt_rebin]*pt_count
    if type(pt_rebin) == list:
        if len(pt_rebin) < pt_count:
            pt_rebin += [pt_rebin[-1]] * (pt_count - len(pt_rebin))

# setup function for pt ranges
def setup_pt(axis, pt_bins):
    pt_count, pt_range = None, None
    if type(pt_bins[0]) == tuple:
        pt_count = len(pt_bins)
        pt_range = pt_bins
    else:
        pt_range, pt_count = get_pt_bin_ranges(axis, pt_bins)
    return pt_range, pt_count

# setup of the fitting ranges for the given pt bins
def setup_fit_range(fitrange, pt_count):
    if type(fitrange) != list:
        fitrange = [(-fitrange, fitrange)]*pt_count
    elif type(fitrange) == list:
        if type(fitrange[0]) not in [list, tuple]:
            if len(fitrange) < pt_count:
                fitrange += [fitrange[-1]]*(pt_count - len(fitrange))
            fitrange = [(-fitrange[n], fitrange[n]) for n in range(len(fitrange))]
    return fitrange

# setup fit range for dca xy and z templates
def setup_fit_range_xyz(fitrange, pt_count):
    if type(fitrange) in [int, float]:
        fitrange = setup_fit_range(fitrange, pt_count)
        fitrange = [fitrange, fitrange]
    if type(fitrange) == list and len(fitrange) == 2:
        fitrange1 = setup_fit_range(fitrange[0], pt_count)
        fitrange2 = setup_fit_range(fitrange[1], pt_count)
        fitrange = [fitrange1, fitrange2]
    return fitrange

# setup signal range
def setup_signal_range(signal_range):
    if type(signal_range) in [int, float]:
        signal_range = [-signal_range, signal_range]
    return signal_range

# setup signal range for xy and z
def setup_signal_range_xyz(signal_range):
    if type(signal_range) in [int, float]:
        signal_range = [-signal_range, signal_range]
    if type(signal_range) == list and len(signal_range) == 2:
        if type(signal_range[0]) in [int, float]:
            signal_range = [signal_range, signal_range]
        elif type(signal_range[0]) == list and len(signal_range) < 2:
            signal_range = [signal_range, signal_range]
    return signal_range

# setup signal range as bin number
def setup_signal_bin_range(axis, signal_range):
    signal_range = setup_signal_range(signal_range)
    bin_low = axis.FindBin(signal_range[0])
    bin_up = find_bin_reduce_on_lower_edge(axis, signal_range[1])
    return [bin_low, bin_up]

# setup fraction dictionaries for individual template
def setup_temp_dict(temp_dict, pt_count):
    if type(temp_dict) == dict:
        temp_dict = [temp_dict]
    for ntemp in range(len(temp_dict)):
        init = temp_dict[ntemp]['temp_init']
        limits = temp_dict[ntemp]['temp_limits']

        if type(init) in [int, float]:
            init = [init] * pt_count
        elif type(init) == list:
            if len(init) < pt_count:
                init += [init[-1]] * (pt_count - len(init))

        if type(limits) in [int, float]:
            if limits:
                limits = [(0, limits)] * pt_count
        elif type(limits) == list:
            for n in range(len(limits)):
                if type(limits[n]) in [int, float]:
                    limits[n] = (0, limits[n])
            if len(limits) < pt_count:
                limits += [limits[-1]] * (pt_count - len(limits))

        temp_dict[ntemp]['temp_init'] = init
        temp_dict[ntemp]['temp_limits'] = limits

    return temp_dict

# obtain fractions from purity th1 plot
def get_fractions_from_purity(hpurity, pt_bins):
    pt_range, pt_count = setup_pt(hpurity.GetXaxis(), pt_bins)
    purity = []
    for npt in range(pt_count):
        pt_int = hpurity.Integral(pt_range[npt][0], pt_range[npt][1])
        purity.append(1 - (pt_int / (pt_range[npt][1] + 1 - pt_range[npt][0])))
    return purity

# project thnSparse from 8 dims to 4 dims, same structure like templates
def project_th8_to_th4(data):
    th8 = data.Clone()
    axis_pt = 2; axis_xy = 5; axis_z = 6; axis_mult = 0
    ndims = np.array([axis_pt, axis_xy, axis_z, axis_mult], dtype = ctypes.c_int)
    th4 = th8.Projection(4, ndims, "")
    return th4.Clone("hDCA_th4")

# apply xy, z cuts to pt for thn8 histos
def apply_dca_cuts_thn(data, dim, fit_range, pt_bins):
    labels = ['', 'xy', 'z']        # dim 1 = xy, dim 2 = z
    axis_pt = data.GetAxis(0)

    pt_range, pt_count = setup_pt(axis_pt, pt_bins)
    fit_range = setup_fit_range_xyz(fit_range, pt_count)

    # prepare empty histogram
    data_proj = data.Clone(f"hDCA{labels[dim]}").Projection(dim, 0).Clone(f"hDCA{labels[dim]}")
    data_proj.Reset()
    data_proj.SetDirectory(0)

    # generate input th2
    proj_list = ROOT.TList()
    for npt in range(pt_count):
        thn = data.Clone()
        axis_pt = thn.GetAxis(0)
        axis_xy = thn.GetAxis(1)
        axis_z  = thn.GetAxis(2)
        axis_pt.SetRange(pt_range[npt][0], pt_range[npt][1])
        axis_xy.SetRange(axis_xy.FindBin(fit_range[0][npt][0]), find_bin_reduce_on_lower_edge(axis_xy, fit_range[0][npt][1]))
        axis_z.SetRange(axis_z.FindBin(fit_range[1][npt][0]), find_bin_reduce_on_lower_edge(axis_z, fit_range[1][npt][1]))
        proj = thn.Projection(dim, 0)
        proj.SetName(f"proj_{labels[dim]}_{npt}")
        proj_list.Add(proj)
    data_proj.Merge(proj_list)

    return data_proj

# convert fraction to template weight
def convert_fraction_to_weight(data, temp, fraction, fit_range, signal_range):
    bin_signal_low, bin_signal_up = setup_signal_bin_range(data.GetXaxis(), signal_range)
    bin_fit_low, bin_fit_up = setup_signal_bin_range(data.GetXaxis(), fit_range)
    data_fit = data.Integral(bin_fit_low, bin_fit_up)
    temp_fit = temp.Integral(bin_fit_low, bin_fit_up)
    data_signal = data.Integral(bin_signal_low, bin_signal_up)
    temp_signal = temp.Integral(bin_signal_low, bin_signal_up)

    weight = fraction * ((data_signal * temp_fit) / (data_fit * temp_signal))
    return weight

# convert list of fractions to the template weights
def convert_fractions_to_weights(data, temp, fractions, fit_range, signal_range, pt_bins):
    pt_range, pt_count = setup_pt(data.GetXaxis(), pt_bins)
    data = data.ProjectionY()
    temp_slices = []
    for npt in range(pt_count):
        temp_slices.append(temp.ProjectionY(f"fraction_proj_{npt}", pt_range[npt][0], pt_range[npt][1]))

    weights = []
    for proj, fraction in zip(temp_slices, fractions):
        weights.append(convert_fraction_to_weight(data, proj, fraction, fit_range, signal_range))

    return weights

# returns th2 (dca-pt) distributions from a list of thn's to be used as input for the combined fitter
def get_input_templates_from_thn(name_temp, tdir, namelist, fitrange, pt_bins):
    dim_counter = range(2)
    temp_counter = range(len(namelist))
    labels = ['xy', 'z']

    # get templates
    temps_file = FR.FileReader(name_temp, tdir)
    temps = [temps_file.GetHisto(f"hDCAxy_{name}", "Tracks_MC") for name in namelist]

    # setup ranges
    xAxis = temps[0].GetAxis(0)
    pt_range, pt_count = setup_pt(xAxis, pt_bins)
    fitrange = setup_fit_range_xyz(fitrange, pt_count)

    # prepare empty histograms
    temps_cut = [[] for dim in dim_counter]
    for dim in dim_counter:
        for ntemp in temp_counter:
            temps_cut[dim].append(temps[0].Clone(f"empty_{dim}_{ntemp}").Projection(dim + 1, 0).Clone(f"hDCA{labels[dim]}_{namelist[ntemp]}"))
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
                ax1.SetRange(ax1.FindBin(fitrange[0][npt][0]), find_bin_reduce_on_lower_edge(ax1, fitrange[0][npt][1]))
                ax2 = thn.GetAxis(2)
                ax2.SetRange(ax2.FindBin(fitrange[1][npt][0]), find_bin_reduce_on_lower_edge(ax2, fitrange[1][npt][1]))
                proj = thn.Projection(dim + 1, 0)
                proj.SetName(f"proj_{labels[dim]}_{ntemp}_{npt}")
                proj_list.Add(proj)
            temps_cut[dim][ntemp].Merge(proj_list)

    return temps_cut

# returns absolute fractions of DCAxy and DCAz
def get_rel_yield(name_temps, tdir, namelist, target, pt_bins, fitrange):
    temps = get_input_templates_from_thn(name_temps, tdir, namelist, fitrange, pt_bins)
    xAxis = temps[0][0].GetXaxis()
    temps = temps[0]
    pt_range, pt_count = setup_pt(xAxis, pt_bins)

    temp_counter = range(len(temps))
    for ntemp in temp_counter:
        temps[ntemp] = temps[ntemp].ProjectionX().Clone(f"hDCA_proj_{namelist[ntemp]}")

    temps_total = temps[0].Clone()
    temps_total.Reset()

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

# internal function
def _setup_tgraphs(dca_names, temp_colors, pt_count):
    empty_graphs = [ROOT.TGraph(pt_count - 1) for ntemp in range(len(titles))]
    for ntemp in range(len(titles)):
        empty_graphs[n].SetName(dca_names[n])
        empty_graphs[n].SetTitle(dca_names[n])
        empty_graphs[n].SetLineWidth(2)
        empty_graphs[n].SetLineColor(n + 3)
        empty_graphs[n].SetMarkerStyle(21)
        empty_graphs[n].SetMarkerSize(2)
        empty_graphs[n].SetMarkerColor(n + 3)
        empty_graphs[n].GetXaxis().SetLabelSize(0.05)
        empty_graphs[n].GetYaxis().SetLabelSize(0.05)
        empty_graphs[n].GetXaxis().SetTitleSize(0.05)
        empty_graphs[n].GetXaxis().SetTitle("<p_{T}> (GeV)")
        if dca_names[n] in color_temp:
            empty_graphs[n].SetLineColor(color_temp[dca_names[n]])
            empty_graphs[n].SetMarkerColor(color_temp[dca_names[n]])

def CombinedFit(fname, dir_out, dca_data, dca_templates, dca_names, fit_range, signal_range, pt_bins, pt_rebin, temp_init, temp_limits, temp_fraction, print_canvas):
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
    ROOT.gStyle.SetOptStat(0)
    color_data, color_fit, color_temp = setup_colors()

    # Output file
    ofile = ROOT.TFile(dirOut + "CombinedFit_" + fname, "recreate") if type(fname) == str else fname
    odirs = [ofile.mkdir(f"{names_for_saving[nfit]}_fits") for nfit in fit_counter]

    # initialize list of values if single rebin value was provided
    if pt_rebin and (len(pt_rebin) < len(pt_bins)):
        pt_rebin = pt_rebin*len(pt_bins)

    # setup for the bin ranges for the projections
    # output: [(bin1, bin2), (bin2, bin3), ...]
    pt_range, pt_count = setup_pt(xAxis, pt_bins)

    # setup signal range
    signal_range = setup_signal_range(signal_range)

    # setup for the fraction fitting parameters
    if type(temp_fraction) == dict:
        temp_fraction = [temp_fraction]

    # graph titles
    pt_names = []
    for n in range(pt_count):
        name = "p_{T} #in " + f"[{pt_bins[n]} - {pt_bins[n + 1]}) GeV"
        if pt_rebin:
            name += f" - rebin {pt_rebin[n]}"
        pt_names.append(name)

    # setup of the fitting ranges for the given pt bins
    fit_range = setup_fit_range_xyz(fit_range, pt_count)

    # graph initialization
    total_chi2 = ROOT.TGraph(pt_count)
    total_chi2.SetLineWidth(2)
    total_chi2.GetXaxis().SetLabelSize(0.05)
    total_chi2.GetYaxis().SetLabelSize(0.05)
    total_chi2.GetYaxis().SetTitle("chi2/ndf")
    total_chi2.GetXaxis().SetTitle("<p_{T}> (GeV)")
    total_chi2.GetXaxis().SetTitleSize(0.05)
    total_chi2.GetYaxis().SetTitleSize(0.05)

    temp_graph = [[ROOT.TGraph(pt_count - 1) for ntemp in temp_counter] for nfit in fit_counter]
    for nfit in fit_counter:
        for ntemp in temp_counter:
            temp_graph[nfit][ntemp].SetName(f"{names_for_saving[nfit]}_{dca_names[ntemp]}")
            temp_graph[nfit][ntemp].SetTitle(dca_names[ntemp])
            temp_graph[nfit][ntemp].SetLineColor(ntemp + 3)
            temp_graph[nfit][ntemp].SetLineWidth(2)
            temp_graph[nfit][ntemp].SetMarkerColor(ntemp + 3)
            temp_graph[nfit][ntemp].SetMarkerSize(2)
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
    for n in range(pt_count):
        # canvas for fitting
        canvas_arr = [ROOT.TCanvas(f"{names_for_saving[nfit]}_canvas_{n + 1}", f"{names_for_saving[nfit]}_canvas_{n + 1}", 1440, 1080) for nfit in fit_counter]
        ROOT.gPad.SetLogy()

        # setup data
        data = []
        for nfit in fit_counter:
            data.append(dca_data[nfit].ProjectionY(f"{names_for_saving[nfit]}_{n + 1}", pt_range[n][0], pt_range[n][1]))
            if pt_rebin and pt_rebin[n] != 1:
                data[nfit].Rebin(pt_rebin[n])
            data[nfit].SetAxisRange(fit_range[nfit][n][0], fit_range[nfit][n][1])
            bin_low = data[nfit].FindBin(fit_range[nfit][n][0])
            bin_up = find_bin_reduce_on_lower_edge(data[nfit].GetXaxis(), fit_range[nfit][n][1])
            data_ent[nfit].append(data[nfit].Integral(bin_low, bin_up))
            if data[nfit].Integral():
                data[nfit].Scale(1. / data[nfit].Integral())
            data[nfit].SetLineColor(color_data)
            data[nfit].SetLineWidth(2)

        # setup templates
        temp_hist = [[] for nfit in fit_counter]
        for nfit in fit_counter:
            for ntemp in temp_counter:
                temp_hist[nfit].append(dca_templates[nfit][ntemp].ProjectionY(f"{names_for_saving[nfit]}_{dca_names[ntemp]}_{n + 1}", pt_range[n][0], pt_range[n][1]))
                if pt_rebin and pt_rebin[n] != 1:
                    temp_hist[nfit][ntemp].Rebin(pt_rebin[n])
                temp_hist[nfit][ntemp].SetAxisRange(fit_range[nfit][n][0], fit_range[nfit][n][1])
                temp_hist[nfit][ntemp].SetTitle(pt_names[n])
                if temp_hist[nfit][ntemp].Integral():
                    temp_hist[nfit][ntemp].Scale(1. / temp_hist[nfit][ntemp].Integral())

        # setup data for wrappers
        data_opt = ROOT.Fit.DataOptions()
        data_range = [ROOT.Fit.DataRange(), ROOT.Fit.DataRange()]
        for nfit in fit_counter:
            data_range[nfit].SetRange(fit_range[nfit][n][0], fit_range[nfit][n][1])

        data_wrapped = [ROOT.Fit.BinData(data_opt, data_range[nfit]) for nfit in fit_counter]
        for nfit in fit_counter:
            ROOT.Fit.FillData(data_wrapped[nfit], data[nfit])

        # fit functions
        fit_objs = [ftotal(data[nfit], temp_hist[nfit]) for nfit in fit_counter]
        fit_funcs = [ROOT.TF1(f"{names_for_saving[nfit]}_ftot", fit_objs[nfit], fit_range[nfit][n][0], fit_range[nfit][n][1], temp_count) for nfit in fit_counter]
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
            data[nfit].Fit(f"{names_for_saving[nfit]}_ftot", "S, N, R, M", "", fit_range[nfit][n][0], fit_range[nfit][n][1])

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
            htot[nfit].SetLineColor(color_fit)
            htot[nfit].SetLineWidth(2)

        for nfit in fit_counter:
            canvas_arr[nfit].cd()
            ROOT.gPad.SetLogy()
            data[nfit].SetTitle(pt_names[n])
            data[nfit].Draw("h")
            htot[nfit].Draw("h same")
            for ntemp in temp_counter:
                temp_hist[nfit][ntemp].Draw("h same")

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
            if print_canvas:
                canvas_arr[nfit].SaveAs(f"{canvas_arr[nfit].GetName()}.png")
        ofile.cd()

        # fractions
        for nfit in fit_counter:
            bin_low = htot[nfit].FindBin(signal_range[0])
            bin_up  = find_bin_reduce_on_lower_edge(htot[nfit].GetXaxis(), signal_range[1])
            temp_ent[nfit].append(htot[nfit].Integral(bin_low, bin_up))
            for ntemp in temp_counter:
                bin_low = temp_hist[nfit][ntemp].FindBin(signal_range[0])
                bin_up  = find_bin_reduce_on_lower_edge(temp_hist[nfit][ntemp].GetXaxis(), signal_range[1])
                tmp = temp_hist[nfit][ntemp].Integral(bin_low, bin_up)
                if temp_ent[nfit][n]:
                    par_ent[nfit][ntemp].append(tmp / temp_ent[nfit][n])
                else:
                    par_ent[nfit][ntemp].append(tmp)

        # fill qa graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1] + 1)) / 2
        for nfit in fit_counter:
            chi_graph[nfit].SetPoint(n, pt_avg, fit_funcs[nfit].GetChisquare() / (fit_funcs[nfit].GetNDF()))
            for ntemp in temp_counter:
                temp_graph[nfit][ntemp].SetPoint(n, pt_avg, par_ent[nfit][ntemp][n])
            fit_funcs[nfit].SetFitResult(result, np.array(temp_counter, dtype=np.int32))
        total_chi2.SetPoint(n, pt_avg, fit_funcs[0].GetChisquare() / fit_funcs[0].GetNDF())

    # canvas for chi2
    chi_canvas = [ROOT.TCanvas(f"{names_for_saving[nfit]}_chi_canvas", f"{names_for_saving[nfit]}_chi_canvas", 1440, 1080) for nfit in fit_counter]
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
    chi2_canvas = ROOT.TCanvas("total chi2/ndf", "total chi2/ndf", 1440, 1080)
    chi2_canvas.SetMargin(0.15, 0.05, 0.2, 0.05)
    chi2_canvas.cd()
    total_chi2.SetTitle()
    total_chi2.Draw()
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
        temp_graph[0][ntemp].Write(f"hDCA_{dca_names[ntemp]}")
    for ntemp in temp_counter:
        averages[0][ntemp].Write()
    total_chi2.Write("total chi2/ndf")
    fractions_canvas[nfit].Write("hDCA_fractions")
    if print_canvas:
        chi2_canvas.SaveAs("total_chi2-ndf.png")
        fractions_canvas[nfit].SaveAs("hDCA_fractions.png")
    ofile.Close()

def TemplateFit(fname, dca_data, dca_templates, dcacpa, dca_names, fit_range, signal_range, pt_bins, pt_rebin, dirOut, temp_init, temp_limits, temp_fraction, print_canvas):
    xAxis = dca_data.GetXaxis()
    yAxis = dca_data.GetYaxis()
    dca_ent = len(dca_templates)

    # renaming stupid names
    for n, name in enumerate(dca_names):
        if name == "SecondaryDaughterLambda":
            dca_names[n] = "SecLambda"
        if name == "SecondaryDaughterSigmaplus":
            dca_names[n] = "SecSigmaPlus"

    # Output file
    ofile = ROOT.TFile(dirOut + "TemplateFit_" + fname, "recreate") if type(fname) == str else fname

    # color palette
    ROOT.gStyle.SetOptStat(0)
    color_data, color_fit, color_temp = setup_colors()

    # initialize list of values if single rebin value was provided
    if pt_rebin and (len(pt_rebin) < len(pt_bins)):
        pt_rebin = pt_rebin*len(pt_bins)

    # parameter setup
    pt_range, pt_ent = setup_pt(xAxis, pt_bins)
    fitrange = setup_fit_range(fit_range, pt_ent)
    signal_range = setup_signal_range(signal_range)
    temp_fraction = setup_temp_dict(temp_fraction, pt_ent)

    # titles of the graphs
    pt_names = []
    for n in range(pt_ent):
        name = "p_{T} #in " + f"[{pt_bins[n]} - {pt_bins[n + 1]}) GeV"
        if pt_rebin:
            name += f" - rebin {pt_rebin[n]}"
        pt_names.append(name)

    # should be removed or expanded for actual cpa use
    if type(dcacpa) == str:
        if dcacpa.lower() == "dca":
            dcacpa = 'dca'
        elif dcacpa.lower() == "cpa":
            dcacpa = 'cpa'
        else:
            print("Error in plot type: \"dca\" or \"cpa\"")

    # graph initialization
    gChi = ROOT.TGraph(pt_ent)
    gChi_signal = ROOT.TGraph(pt_ent)
    m2ent = [0 for i in range(dca_ent)]
    dataEntries = []
    mcEntries = []
    parDCA_mc = [[] for n in range(dca_ent)]
    gDCA_mc = [ROOT.TGraph(pt_ent - 1) for ntemp in range(dca_ent)]
    for n in range(dca_ent):
        gDCA_mc[n].SetName(dca_names[n])
        gDCA_mc[n].SetTitle(dca_names[n])
        gDCA_mc[n].SetLineWidth(2)
        gDCA_mc[n].SetLineColor(n + 3)
        gDCA_mc[n].SetMarkerStyle(21)
        gDCA_mc[n].SetMarkerSize(2)
        gDCA_mc[n].SetMarkerColor(n + 3)
        gDCA_mc[n].GetXaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetYaxis().SetLabelSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitleSize(0.05)
        gDCA_mc[n].GetXaxis().SetTitle("<p_{T}> (GeV)")
        if dca_names[n] in color_temp:
            gDCA_mc[n].SetLineColor(color_temp[dca_names[n]])
            gDCA_mc[n].SetMarkerColor(color_temp[dca_names[n]])

    ### main loop for fitting ###
    for n in range(pt_ent):
        # canvas for fitting
        canvas = ROOT.TCanvas("canvas_" + str(n + 1), "canvas_" + str(n + 1))
        canvas.SetCanvasSize(1440, 1080)
        canvas.cd()
        ROOT.gPad.SetLogy()

        # data
        data = dca_data.ProjectionY(f"hDCAxy_{n + 1}", pt_range[n][0], pt_range[n][1])
        if pt_rebin and pt_rebin[n] != 1:
            data.Rebin(pt_rebin[n])
        data.SetAxisRange(fitrange[n][0], fitrange[n][1])
        bin_low = data.FindBin(fitrange[n][0])
        bin_up = find_bin_reduce_on_lower_edge(data.GetXaxis(), fitrange[n][1])
        data_int = data.Integral(bin_low, bin_up)
        dataEntries.append(data_int)
        if data.Integral():
            data.Scale(1. / data.Integral())

        # MC templates
        hDCA_mc = []
        for i in range(dca_ent):
            hDCA_mc.append(dca_templates[i].ProjectionY(f"{dca_names[i]}_{n + 1}", pt_range[n][0], pt_range[n][1]))
            if pt_rebin and pt_rebin[n] != 1:
                hDCA_mc[i].Rebin(pt_rebin[n])
            hDCA_mc[i].SetAxisRange(fitrange[n][0], fitrange[n][1])
            hDCA_mc[i].SetTitle(pt_names[n])
            if hDCA_mc[i].Integral():
                hDCA_mc[i].Scale(1. / hDCA_mc[i].Integral())

        # fit function
        adj = ftotal(data, hDCA_mc)
        ftot = ROOT.TF1("ftot", adj, fitrange[n][0], fitrange[n][1], dca_ent)

        # initialize fitting parameters
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
                            if entry['temp_limits'] == 0:
                                ftot.FixParameter(i, entry['temp_init'][n])
                            if type(entry['temp_limits']) == list and entry['temp_limits'][n][1] == 0:
                                ftot.FixParameter(i, entry['temp_init'][n])
                        else:
                            if entry['temp_limits'] == 0:
                                ftot.FixParameter(i, 0)
                        if entry['temp_limits']:
                            if type(entry['temp_limits']) == list and entry['temp_limits'][n][1] == 0:
                                ftot.FixParameter(i, 0)
                            else:
                                ftot.SetParLimits(i, entry['temp_limits'][n][0], entry['temp_limits'][n][1])

        # set data error to sqrt of sum of weights of data and temps
        for npt in range(data.GetNbinsX()):
            error = data.GetBinError(npt)**2
            for ntemp in range(dca_ent):
                error += hDCA_mc[i].GetBinError(npt)**2
            error = error**0.5
            data.SetBinError(npt, error)

        # draw histograms
        data.Fit("ftot", "S, N, R, M", "", fitrange[n][0], fitrange[n][1])
        data.SetLineColor(color_data)
        data.SetLineWidth(2)
        data.SetTitle(pt_names[n])
        data.Draw("h")

        for i in range(dca_ent):
            hDCA_mc[i].Scale(ftot.GetParameter(i))
            hDCA_mc[i].SetLineColor(i + 3)
            hDCA_mc[i].SetLineWidth(2)
            if dca_names[i] in color_temp:
                hDCA_mc[i].SetLineColor(color_temp[dca_names[i]])
                hDCA_mc[i].SetMarkerColor(color_temp[dca_names[i]])
            hDCA_mc[i].Draw("h same")

        htot = hDCA_mc[0].Clone("htot_" + str(n + 1))
        for i in range(1, dca_ent):
            htot.Add(hDCA_mc[i])
        htot.SetLineColor(color_fit)
        htot.SetLineWidth(2)
        htot.Draw("h same")

        # legend
        legend = ROOT.TLegend(0.15, 0.40, 0.35, 0.8)
        legend.SetTextSize(0.05)
        legend.SetBorderSize(0)
        legend.SetLineColor(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)

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
        if print_canvas:
            canvas.SaveAs(f"{canvas.GetName()}.png")
        ofile.cd()

        # fractions
        htot_int = htot.Integral(htot.FindBin(signal_range[0]), find_bin_reduce_on_lower_edge(htot.GetXaxis(), signal_range[1]))
        mcEntries.append(htot_int)
        for i in range(dca_ent):
            bin_low = hDCA_mc[i].FindBin(signal_range[0])
            bin_up = find_bin_reduce_on_lower_edge(hDCA_mc[i].GetXaxis(), signal_range[1])
            tmp = hDCA_mc[i].Integral(bin_low, bin_up)
            if mcEntries[n]:
                parDCA_mc[i].append(tmp / mcEntries[n])
            else:
                parDCA_mc[i].append(tmp)

        # fill x and y lists for graphs
        pt_avg = (xAxis.GetBinLowEdge(pt_range[n][0]) + xAxis.GetBinLowEdge(pt_range[n][1] + 1)) / 2
        for i in range(dca_ent):
            gDCA_mc[i].SetPoint(n, pt_avg, parDCA_mc[i][n])
        if ftot.GetNDF():
            gChi.SetPoint(n, pt_avg, ftot.GetChisquare() / (ftot.GetNDF()))

        chi2 = 0
        bin_low, bin_up = setup_signal_bin_range(data.GetXaxis(), signal_range)
        for dca_bin in range(data.GetNbinsX()):
            if dca_bin >= bin_low and dca_bin <= bin_up:
                data_value = data.GetBinContent(dca_bin)
                data_error = data.GetBinError(dca_bin)
                fit_value = ftot.Eval(data.GetBinCenter(dca_bin))
                if not fit_value:
                    continue
                chi2 += ((data_value - fit_value) / data_error)**2
                #chi2 += (data_value - fit_value)**2 / data_error
                #chi2 += ((data_value - fit_value) / fit_value)**2
                #chi2 += (data_value - fit_value)**2 / fit_value
        chi2ndf = chi2 / (bin_up - bin_low - dca_ent)
        gChi_signal.SetPoint(n, pt_avg, chi2ndf)

    # calculate final pT result
    m2tot = 0
    for i in range(pt_ent):
        m2tot += dataEntries[i]
        for j in range(dca_ent):
            m2ent[j] += dataEntries[i]*parDCA_mc[j][i]

    # chi2 graph
    gChi.SetLineColor(1)
    gChi.SetLineWidth(2)
    gChi.GetXaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetLabelSize(0.05)
    gChi.GetYaxis().SetTitle("chi2/ndf")
    gChi.GetXaxis().SetTitle("<p_{T}> (GeV)")
    gChi.GetXaxis().SetTitleSize(0.05)
    gChi.GetYaxis().SetTitleSize(0.05)

    gChi_signal.SetLineColor(1)
    gChi_signal.SetLineWidth(2)
    gChi_signal.GetXaxis().SetLabelSize(0.05)
    gChi_signal.GetYaxis().SetLabelSize(0.05)
    gChi_signal.GetYaxis().SetTitle("chi2/ndf")
    gChi_signal.GetXaxis().SetTitle("<p_{T}> (GeV)")
    gChi_signal.GetXaxis().SetTitleSize(0.05)
    gChi_signal.GetYaxis().SetTitleSize(0.05)

    chi2_canvas = ROOT.TCanvas("total chi2/ndf", "total chi2/ndf", 1440, 1080)
    chi2_canvas.SetMargin(0.15, 0.05, 0.2, 0.05)
    chi2_canvas.cd()
    gChi.SetTitle()
    gChi.Draw()

    # canvas for all fractions
    fractions = ROOT.TCanvas("fractions", "fractions", 1440, 1080)
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    fractions.cd()
    fractions.SetFillColor(0)
    fractions.SetBorderMode(0)
    fractions.SetBorderSize(2)
    fractions.SetFrameBorderMode(0)
    fractions.SetMargin(0.15, 0.05, 0.2, 0.05)

    # fraction graphs
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
        line = ROOT.TLine(xAxis.GetBinLowEdge(pt_range[i][0]), parDCA_mc[0][i], xAxis.GetBinLowEdge(pt_range[i][1] + 1), parDCA_mc[0][i])
        line.DrawClone("same")

    # average
    primAvg = ROOT.TF1("primAvg", "[0]", 0, 5)
    primAvg.SetParameter(0, m2ent[0] / (m2tot + 1))
    primAvg.SetLineColor(12)
    primAvg.SetLineStyle(2)
    primAvg.Draw("same")

    # dump info
    for i in range(dca_ent):
        print("%-12s %.6f" % (dca_names[i], m2ent[i] / (m2tot + 1)))

    # save output
    ofile.cd()
    for i in range(dca_ent):
        gDCA_mc[i].Write()
    gChi.Write("chi2/ndf")
    gChi_signal.Write("chi2/ndf signal")
    fractions.Write("hDCA_fractions")
    if print_canvas:
        chi2_canvas.SaveAs("total_chi2-ndf.png")
        fractions.SaveAs("hDCA_fractions.png")
    primAvg.Write()


