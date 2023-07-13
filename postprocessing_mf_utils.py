import sys
import matplotlib
import scienceplots
import pandas as pd
import numpy as np
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib.ticker as ticker
from cycler import cycler
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from spydiff import export_difmap_model_from_tuples, join_difmap_models, create_difmap_file_from_single_component

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


def count_jet_components(df, n_bands):
    return int((len(df.columns) - n_bands - 2*n_bands - 6)/7)


def plot_posterior_samples_on_map(posterior_file, n_bands, ra_lims=(-5, 5), dec_lims=(-5, 5)):
    df = convert_posterior_file_to_pandas_df(posterior_file)
    n_jc = count_jet_components(df, n_bands)
    print("Found {} jet components".format(n_jc))
    fig, axes = plt.subplots(1, 1)
    axes.set_xlim(ra_lims)
    axes.set_ylim(dec_lims)
    idx_0_jc = int(n_bands + 2*n_bands + 6)
    print("Jet components parameters start at idx = {}".format(idx_0_jc))
    for _, row in df.iterrows():
        for i_jc in range(n_jc):
            ra = row[idx_0_jc + i_jc*7]
            dec = row[idx_0_jc + 1 + i_jc*7]
            size = np.exp(row[idx_0_jc + 2 + i_jc*7])
            # print("Component #{}. ra = {:.2f}, dec = {:.2f}, size = {:.2f}".format(i_jc, ra, dec, size))
            circle = Circle(xy=(ra, dec), radius=size, facecolor='green', alpha=0.01)
            axes.scatter(ra, dec, s=5, alpha=0.2, color=colors[i_jc])
            # axes.add_patch(circle)
    axes.set_xlabel("RA, mas")
    axes.set_ylabel("DEC, mas")
    axes.invert_xaxis()
    axes.set_aspect("equal")
    plt.show()
    return fig


if __name__ == "__main__":
    n_bands = 3
    posterior_file = "/home/ilya/github/bam/mf/posterior_sample.txt"
    df = convert_posterior_file_to_pandas_df(posterior_file)
    n_jc = count_jet_components(df, n_bands)
    fig = plot_posterior_samples_on_map(posterior_file, n_bands, ra_lims=(-5, 5), dec_lims=(-5, 5))
