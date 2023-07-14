import sys
import matplotlib
import scienceplots
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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


if __name__ == "__main__":
    n_bands = 3
    jitter = False
    freqs_ghz = (15.4, 23.8, 43.2)
    posterior_file = "/home/ilya/github/bam/mf_noj/posterior_sample_5jc_noj.txt"
    # posterior_file = "/home/ilya/github/bam/mf/posterior_sample_6jc.txt"
    df = convert_posterior_file_to_pandas_df(posterior_file)
    n_jc = count_jet_components(df, n_bands, jitter)
    fig = plot_posterior_samples_on_map(posterior_file, n_bands, freqs_ghz, ra_lims=(-10, 10), dec_lims=(-10, 10),
                                        alpha_jet=0.002, alpha_core=0.1, jitter=jitter, each=1)
