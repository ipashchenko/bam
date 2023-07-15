import os
import sys
import matplotlib
import scienceplots
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines
from matplotlib.patches import Circle
import matplotlib.ticker as ticker
from cycler import cycler
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from spydiff import export_difmap_model_from_tuples, join_difmap_models, create_difmap_file_from_single_component

matplotlib.use("TkAgg")
# For tics and line widths. Re-write colors and figure size later.
plt.style.use('science')
# Default color scheme
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
matplotlib.rcParams['axes.prop_cycle'] = cycler('color', colors)
# Default figure size
matplotlib.rcParams['figure.figsize'] = (6.4, 4.8)

label_size = 18
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def convert_posterior_file_to_pandas_df(posterior_file):
    from io import StringIO
    with open(posterior_file, "r") as fo:
        line = fo.readline()
    line = line.strip()
    line = line.lstrip("#")
    names = line.split('\t')
    names = [name.strip() for name in names]
    new_header = " ".join(names)
    with open(posterior_file, "r") as fo:
        lines = fo.readlines()
    new_file = new_header + "\n"
    for line in lines:
        if line.startswith("#"):
            continue
        else:
            new_file += line

    df = pd.read_csv(StringIO(new_file), delim_whitespace=True)
    for old in ("dx", "dy", "logsize", "lognu_max", "logS_max", "alpha_thick", "alpha_thin"):
        df = df.rename(columns={old: old + ".0"})

    return df


def count_jet_components(df, n_bands, jitter):
    if jitter:
        n_jitters = n_bands
    else:
        n_jitters = 0
    return int((len(df.columns) - n_jitters - 2*n_bands - 6)/7)


def plot_posterior_samples_on_map(posterior_file, n_bands, freqs_ghz, ra_lims=(-5, 5), dec_lims=(-5, 5),
                                  alpha_jet=0.01, alpha_core=0.2, jitter=True, each=1):
    freqs_ghz = np.sort(freqs_ghz)
    df = convert_posterior_file_to_pandas_df(posterior_file)
    n_jc = count_jet_components(df, n_bands, jitter)
    print("Found {} jet components".format(n_jc))
    fig, axes = plt.subplots(1, 1)
    axes.set_xlim(ra_lims)
    axes.set_ylim(dec_lims)
    if jitter:
        n_jitters = n_bands
    else:
        n_jitters = 0
    idx_0_jc = int(n_jitters + 2*n_bands + 6)
    idx_0_cc = int(n_jitters + 2*n_bands)
    axes.scatter(0, 0, marker="+", color="black", s=20)
    for i, freq_ghz in enumerate(freqs_ghz):
        axes.scatter([], [], alpha=1.0, color=colors[i], label="{} GHz core".format(freq_ghz))
    for _, row in df.iterrows():
        if _ % each != 0:
            continue
        a = row[idx_0_cc + 0]
        PA = row[idx_0_cc + 1]
        k_r = row[idx_0_cc + 3]
        ra_core_shifts = list()
        dec_core_shifts = list()
        for i, freq_ghz in enumerate(freqs_ghz):
            distance = a*freq_ghz**(-1./k_r)
            ra_core_shift = distance*np.sin(PA)
            dec_core_shift = distance*np.cos(PA)
            ra_core_shifts.append(ra_core_shift)
            dec_core_shifts.append(dec_core_shift)
            axes.scatter(ra_core_shift, dec_core_shift, s=10, alpha=alpha_core, color=colors[i])
        axes.plot(ra_core_shifts, dec_core_shifts, lw=0.5, color="gray", alpha=0.5)
        for i_jc in range(n_jc):
            ra = row[idx_0_jc + i_jc*7]
            dec = row[idx_0_jc + 1 + i_jc*7]
            size = np.exp(row[idx_0_jc + 2 + i_jc*7])
            # print("Component #{}. ra = {:.2f}, dec = {:.2f}, size = {:.2f}".format(i_jc, ra, dec, size))
            circle = Circle(xy=(ra, dec), radius=size/2, facecolor='gray', alpha=alpha_jet)
            # axes.scatter(ra, dec, s=5, alpha=alpha_jet, color="chocolate")
            axes.add_patch(circle)
    axes.set_xlabel("RA, mas")
    axes.set_ylabel("DEC, mas")
    axes.invert_xaxis()
    axes.set_aspect("equal")
    plt.legend()
    plt.show()
    return fig


def optically_thin_spectr(nu, logS_max, lognu_max, alpha_thin, alpha_thick):
    return logS_max + alpha_thick*np.log(nu/np.exp(lognu_max)) - np.log(1 - np.exp(-1)) + np.log(1 - np.exp(-(nu/np.exp(lognu_max))**(alpha_thin - alpha_thick)))


def plot_components_info(df, n_bands, jitter, freqs_ghz, save_dir=None, opacity_each_line=0.05):
    freqs_ghz = sorted(freqs_ghz)
    n_jc = count_jet_components(df, n_bands, jitter)
    nu_grid = np.logspace(np.log10(np.min(freqs_ghz))-0.25, np.log10(np.max(freqs_ghz))+0.25, 100, base=10)
    for i in range(n_jc):
        RA, DEC = np.median(df[f"dx.{i}"]), np.median(df[f"dy.{i}"])
        logS_max = df[f"logS_max.{i}"].values
        lognu_max = df[f"lognu_max.{i}"].values
        alpha_thin = df[f"alpha_thin.{i}"].values
        alpha_thick = df[f"alpha_thick.{i}"].values
        print("For component at RA = {:.2f}, DEC = {:.2f} the fluxes are :".format(RA, DEC))
        for nu in freqs_ghz:
            print("For nu = {} GHz: {:.3f} Jy".format(nu, np.median(np.exp(optically_thin_spectr(nu, logS_max, lognu_max, alpha_thin, alpha_thick)))))
        print("Spectral indexes are :")
        print("alpha_thin = {:.2f}, alpha_thick = {:.2f}".format(np.median(df[f"alpha_thin.{i}"]), np.median(df[f"alpha_thick.{i}"])))
        print("nu_max = {:.2f} GHz".format(np.exp(np.median(df[f"lognu_max.{i}"]))))
        print("---------------------------------------------------------------")
        fig, axes = plt.subplots(1, 1)
        for j in range(len(df)):
            logflux = optically_thin_spectr(nu_grid, logS_max[j], lognu_max[j], alpha_thin[j], alpha_thick[j])
            axes.plot(nu_grid, np.exp(logflux), alpha=opacity_each_line, color="r", lw=2)
        axes.set_xlabel(r"frequency, GHz")
        axes.set_ylabel(r"Flux, Jy")
        axes.set_xscale("log", base=10)
        axes.set_yscale("log", base=10)
        axes.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
        axes.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y), 0)))).format(y)))
        axes.set_title("RA: {:.2f}, DEC: {:.2f}, size: {:.2f} mas".format(RA, DEC, np.median(np.exp(df[f"logsize.{i}"]))))

        # remove the minor ticks
        axes.yaxis.set_minor_formatter(ticker.NullFormatter())
        for icolor, nu in enumerate(freqs_ghz):
            axes.axvline(nu, lw=1, color="k", ls="--")
        if save_dir is not None:
            plt.savefig(os.path.join(save_dir, "spectra_component_RA_{:.2f}_DEC_{:.2f}.png".format(RA, DEC)), bbox_inches="tight", dpi=300)
        plt.show()


if __name__ == "__main__":
    jitter = True
    # freqs_ghz = (15.4, 23.8, 43.2)
    freqs_ghz = (4.6, 5.0, 8.1, 8.43, 15.4, 23.8, 43.2)
    n_bands = len(freqs_ghz)
    posterior_file = "/home/ilya/github/bam/mf/all_bands/posterior_sample.txt"
    # posterior_file = "/home/ilya/github/bam/mf/posterior_sample_6jc.txt"
    save_dir = "/home/ilya/github/bam/mf_noj/6jc"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    df = convert_posterior_file_to_pandas_df(posterior_file)
    n_jc = count_jet_components(df, n_bands, jitter)
    plot_components_info(df, n_bands, jitter, freqs_ghz, save_dir=None, opacity_each_line=0.5)
    fig = plot_posterior_samples_on_map(posterior_file, n_bands, freqs_ghz, ra_lims=(-10, 10), dec_lims=(-10, 10),
                                        alpha_jet=0.3, alpha_core=0.5, jitter=jitter, each=1)
