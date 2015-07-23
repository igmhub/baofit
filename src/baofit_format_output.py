"""
Transforms baofit output into readable outputs
In order to get the script working from any directory paste the following into your .bash_profile or .profile  
and reset the shell
    alias baofit_format_output='python $(PATH_TO_BAOFIT)/src/baofit_format_output.py'
Note that you should change $(PATH_TO_BAOFIT) for the actual path to baofit
USAGE:
    baofit_format_output BOSSDR11QSOLyaF.ini
    baofit_format_output -r BOSSDR11QSOLyaF.ini
    baofit_format_output --vmin="-18" --vmax=20 -r BOSSDR11QSOLyaF.ini
    baofit_format_output --vmin="-22" --vmax=25 --vmax-res=6 -r BOSSDR11QSOLyaF_k_no_rad_all_fixed.ini
    
OUTPUTS:
    *_results.tex - tex file containing the baofit results
    *_results.pdf - compiled pdf of *_results.tex

created 22 Jul 2015 by Ignasi Perez-Rafols <iprafols at icc.ub.edu>
"""
import argparse
import sys
import os

import numpy as np

import matplotlib.pyplot as plt



#------------------------------------------------------------------------
def main():
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="transforms baofit output into readable outputs")
    parser.add_argument('input', type=str, default=None, help='input file, required')
    parser.add_argument('-r', '--r-squared-weighting', action="store_true", help='plots r^2-weighted cross-correlation')
    parser.add_argument('--vmin', type=float, default=None, help='minimum value in contour plots')
    parser.add_argument('--vmax', type=float, default=None, help='maximum value in contour plots')
    parser.add_argument('--vmin-res', type=float, default=None, help='minimum value in residuals contour plots')
    parser.add_argument('--vmax-res', type=float, default=None, help='maximum value in residuals contour plots')
    args = parser.parse_args()


    # locate baofit output prefix
    output_prefix, pi_bins, sigma_bins = readIniFile(args.input)

    # initialize result file
    fit_results = output_prefix + "results.tex"
    results = open(fit_results, "w")
    initializeTexFile(results)

    # read results
    fit_file = output_prefix + "fit.config"
    results_dict = readFitResults(fit_file)
    
    # tabulate results
    tabulateResults(results, results_dict)

    # make contour plots
    residuals_file = output_prefix + "residuals.dat"
    fig_name = output_prefix + "residuals.eps"
    makeContourPlots(residuals_file, fig_name, pi_bins, sigma_bins, args)

    # add plots to tex file
    includeFigureToTex(results, fig_name)

    # close result file & compile
    results.write("\\end{document}\n")
    results.close()
    os.system("pdflatex " + output_prefix + "results.tex")

#------------------------------------------------------------------------
def includeFigureToTex(results, fig_name):
    try:
        assert type(results) == file
        assert type(fig_name) == str
    except AssertionError:
        print "errror including plot to tex file"
        return

    results.write("\n")
    results.write("\\begin{figure}\n")
    results.write("    \\centering\n")
    results.write("    \\includegraphics[width=\\textwidth]{" + fig_name +  "}\n")
    results.write("    \\caption{From left to right, 2D cross-correlation for the data measurement, the best fit model and the residuals. }\n")
    results.write("\\end{figure}\n")
    results.write("\n")
    
    return


#------------------------------------------------------------------------
def initializeTexFile(results):
    try:
        assert type(results) == file
    except AssertionError:
        print "errror initializing file"
        return
    
    results.write("\\documentclass{article}\n")
    results.write("\\usepackage{geometry}\n")
    results.write("\\geometry{a4paper}\n")
    results.write("\\usepackage{graphicx}\n")
    results.write("\\usepackage{epstopdf}\n")
    results.write("\n\\begin{document}\n")

    return

#------------------------------------------------------------------------
def makeContourPlots(residuals_file, fig_name, pi_bins, sigma_bins, args):
    try:
        assert type(residuals_file) == str
        assert type(fig_name) == str
        assert isinstance(pi_bins, np.ndarray)
        assert isinstance(sigma_bins, np.ndarray)
    except AssertionError:
        print "could not make contour plots"
        return

    # read data
    residuals = np.genfromtxt(residuals_file, names=["bin_index", "axis0_bin_center", "axis1_bin_center", "axis2_bin_center", "r", "mu", "z", "best_fit_model", "data", "diagonal_error"])

    # fill the voids in the data
    aux = 0
    data_list = []
    model_list = []
    r_list = []
    for (i, r, d, m) in zip(residuals["bin_index"], residuals["r"], residuals["data"], residuals["best_fit_model"]):
        while (aux < i):
            data_list.append(np.nan)
            model_list.append(np.nan)
            r_list.append(np.nan)
            aux += 1
        data_list.append(d)
        model_list.append(m)
        r_list.append(r)
        aux += 1
    while aux < sigma_bins.size*pi_bins.size:
        data_list.append(np.nan)
        model_list.append(np.nan)
        r_list.append(np.nan)
        aux += 1

    # reshape data to plot
    shape = (pi_bins.size, sigma_bins.size)
    data = np.reshape(np.array(data_list), shape)
    model = np.reshape(np.array(model_list), shape)
    r = np.reshape(np.array(r_list), shape)

    # contour settings
    cmap = plt.cm.get_cmap('RdYlBu')
    num_colors = 40.0
    if args.vmin == None:
        if args.r_squared_weighting:
            vmin = min(np.amin(residuals["best_fit_model"]*residuals["r"]*residuals["r"]), np.amin(residuals["data"]*residuals["r"]*residuals["r"]))
        else:
            vmin = min(np.amin(residuals["best_fit_model"]), np.amin(residuals["data"]))
    else:
        vmin = args.vmin
    if args.vmax == None:
        if args.r_squared_weighting:
            vmax = max(np.amax(residuals["best_fit_model"]*residuals["r"]*residuals["r"]), np.amax(residuals["data"]*residuals["r"]*residuals["r"]))
        else:
            vmax = max(np.amax(residuals["best_fit_model"]),  np.amax(residuals["data"]))
    else:
        vmax = args.vmax
    step = (vmax - vmin)/num_colors
    levels = np.arange(vmin, vmax + step, step)
    if args.vmin_res == None:
        if args.r_squared_weighting:
            vmin_res = np.amin((residuals["data"]-residuals["best_fit_model"])*residuals["r"]*residuals["r"])
        else:
            vmin_res = np.amin(residuals["data"]-residuals["best_fit_model"])
    else:
        vmin_res = args.vmin_res
    if args.vmax_res == None:
        if args.r_squared_weighting:
            vmax_res = np.amax((residuals["data"]-residuals["best_fit_model"])*residuals["r"]*residuals["r"])
        else:
            vmax_res = np.amax(residuals["data"]-residuals["best_fit_model"])
    else:
        vmax_res = args.vmax_res
    step_res = (vmax_res - vmin_res)/num_colors
    levels_res = np.arange(vmin_res, vmax_res + step_res, step_res)



    fig = plt.figure(figsize=(21, 7))
    # plot data
    ax = fig.add_subplot(1, 3, 1)

    ax.set_xlabel(r"$\sigma {\rm \left(h^{-1}Mpc\right)}$", fontsize=20)
    ax.set_ylabel(r"$\pi {\rm \left(h^{-1}Mpc\right)}$", fontsize=20)
    ax.tick_params(axis='both',which='major',labelsize=20,length=6,width=2)
    ax.tick_params(axis='both',which='minor',labelsize=0,length=4,width=1)
    if args.r_squared_weighting:
        ax.set_title(r"$r^{2}\xi_{\rm data}{\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
        cs = ax.contourf(sigma_bins, pi_bins, data*r*r, levels, cmap = cmap, vmin = vmin, vmax = vmax, fontsize=20)
    else:
        ax.set_title(r"$\xi_{\rm data}{\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
        cs = ax.contourf(sigma_bins, pi_bins, data, levels, cmap = cmap, vmin = vmin, vmax = vmax, fontsize=20)
    cbar = fig.colorbar(cs, ax = ax, shrink = 0.9, format="%.1e")

    # plot model
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.set_xlabel(r"$\sigma {\rm \left(h^{-1}Mpc\right)}$", fontsize=20)
    ax2.tick_params(axis='both',which='major',labelsize=20,length=6,width=2)
    ax2.tick_params(axis='both',which='minor',labelsize=0,length=4,width=1)
    if args.r_squared_weighting:
        cs2 = ax2.contourf(sigma_bins, pi_bins, model*r*r, levels, cmap = cmap, vmin = vmin, vmax = vmax, fontsize=20)
        ax2.set_title(r"$r^{2}\xi_{\rm model}{\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
    else:
        cs2 = ax2.contourf(sigma_bins, pi_bins, model, levels, cmap = cmap, vmin = vmin, vmax = vmax, fontsize=20)
        ax2.set_title(r"$\xi_{\rm model}{\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
    cbar2 = fig.colorbar(cs2, ax = ax2, shrink = 0.9, format="%.1e")

    # plot residuals
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.set_title(r"$r^{2}\left(\xi_{\rm data}-\xi_{\rm model}\right){\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
    ax3.set_xlabel(r"$\sigma {\rm \left(h^{-1}Mpc\right)}$", fontsize=20)
    ax3.tick_params(axis='both',which='major',labelsize=20,length=6,width=2)
    ax3.tick_params(axis='both',which='minor',labelsize=0,length=4,width=1)
    if args.r_squared_weighting:
        ax3.set_title(r"$r^{2}\left(\xi_{\rm data}-\xi_{\rm model}\right){\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
        cs3 = ax3.contourf(sigma_bins, pi_bins, (data-model)*r*r, levels_res, cmap = cmap, vmin = vmin_res, vmax = vmax_res, fontsize=20)
    else:
        ax3.set_title(r"$\xi_{\rm data}-\xi_{\rm model}{\rm \left(h^{-2}Mpc^{2}\right)}$", fontsize=20)
        cs3 = ax3.contourf(sigma_bins, pi_bins, data-model, levels_res, cmap = cmap, vmin = vmin_res, vmax = vmax_res, fontsize=20)
    cbar3 = fig.colorbar(cs3, ax = ax3, shrink = 0.9, format="%.1e")


    fig.subplots_adjust(wspace=0.3)
    fig.tight_layout()
    fig.savefig(fig_name)

    return


#------------------------------------------------------------------------
def readFitResults(fit_file):
    results_dict = {}
    try:
        assert type(fit_file) == str
    except AssertionError:
        print "could not read results"
        return results_dict

    try:
        for line in open(fit_file).readlines():
            cols = line.split(";")
            # check if the parameter was fixed
            if "fix" in line:
                fixed = True
            else:
                fixed = False
            name = cols[0].split("]")[0].split("[")[1]
            value = float(cols[0].split("=")[1])
            error = float(cols[1].split("=")[1])
            results_dict[name] = (value, error, fixed)
        return results_dict
    except IOError:
        print "could not read results"
        return results_dict

#------------------------------------------------------------------------
def readIniFile(input_file):
    try:
        assert type(input_file) == str
        assert input_file[-4:] == ".ini"
    except AssertionError:
        print "could not read ini file"
        return "", np.array(0), np.array(0)

    try:
        output_prefix = ""
        pi = np.array(0)
        sigma = np.array(0)
        
        for line in open(input_file).readlines():
            # read output_prefix
            if line.startswith("output-prefix"):
                cols = line.split()
                output_prefix = cols[-1]
            # read parallel separation binning
            if line.startswith("axis1-bins"):
                axis1 = line.split("=")[-1]
                try:
                    pi = np.array(axis1.split("{")[1].split("}")[0].split(","), dtype=float)
                except IndexError:
                    min_pi = float(axis1.split("[")[1].split(":")[0])
                    max_pi = float(axis1.split("[")[1].split(":")[1].split("]")[0])
                    num_pi_bins = int(axis1.split("[")[1].split(":")[1].split("*")[1])
                    pi = np.arange(min_pi, max_pi, (max_pi - min_pi)/num_pi_bins)
            # read perpendicular separation binning
            if line.startswith("axis2-bins"):
                axis2 = line.split("=")[-1]
                try:
                    sigma = np.array(axis2.split("{")[1].split("}")[0].split(","), dtype=float)
                except IndexError:
                    min_sigma = float(axis2.split("[")[1].split(":")[0])
                    max_sigma = float(axis2.split("[")[1].split(":")[1].split("]")[0])
                    num_sigma_bins = int(axis2.split("[")[1].split(":")[1].split("*")[1])
                    sigma = np.arange(min_sigma, max_sigma, (max_sigma - min_sigma)/num_sigma_bins)

        return output_prefix, pi, sigma
    except IOError:
        print "could not read ini file"
        return "", np.array(0), np.array(0)


#------------------------------------------------------------------------
def tabulateResults(results, results_dict):
    try:
        assert type(results) == file
        assert type(results_dict) == dict
    except AssertionError:
        print "could not tabulate results"
        return

    results.write("\n")
    results.write("\\begin{table}\n")
    results.write("    \\centering\n")
    results.write("    \\begin{tabular}{c|cccccc}\n")
    results.write("        & $\\beta_{F}$ & ${\\rm b_{q}}$ & ${\\rm a}$ & ${\\rm t_{q}}$ & ${\\rm b}_{\\Gamma}{\\rm r_{p}^{2}}$ & $\\lambda$ \\\\ \\hline \n")
    results.write("        fit ")

    for param in ["beta", "bias2", "Rad anisotropy", "Rad quasar lifetime", "Rad strength", "Rad mean free path"]:
        value = results_dict.get(param, ("-"))
        try:
            if value[2]:
                results.write("& $" + "{:.2F}".format(value[0]) + "$ ")
            else:
                results.write("& $" + "{:.2F}".format(value[0]) + " \\pm " + "{:.2F}".format(value[1]) + "$ ")
        except IndexError:
            results.write("& " + value[0]+" ")

    results.write("\\\\ \n")
    results.write("    \\end{tabular}\n")
    results.write("    \\caption{Fit values. Values without error correspond to those parameters that were fixed during the fitting process.}\n")
    results.write("\\end{table}\n")
    results.write("\n")

    return


__name__ = "__main__"
if __name__ == "__main__":
    main()