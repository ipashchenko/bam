import os
import pathlib
import glob
from collections import OrderedDict
import numpy as np
import pandas as pd
import pickle
from itertools import cycle
import seaborn as sns
from corner import corner
from astropy import units as u
from astropy import constants as const
from sklearn.cluster import DBSCAN, HDBSCAN, OPTICS, cluster_optics_dbscan, SpectralClustering
import ehtim as eh
import matplotlib
# import scienceplots
matplotlib.use("Agg")
# import scienceplots
import matplotlib.pyplot as plt
from cycler import cycler
# For tics and line widths.
# plt.style.use('science')
# Default color scheme
matplotlib.rcParams['axes.prop_cycle'] = cycler('color', ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
# Default figure size
matplotlib.rcParams['figure.figsize'] = (6.4, 4.8)
label_size = 14
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from matplotlib.patches import Ellipse, Circle

from data_utils import radplot, gaussian_circ_ft, gaussian_ell_ft, optically_thin_sphere_ft, point_ft
import sys
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from from_fits import create_clean_image_from_fits_file
from image import plot as iplot
from spydiff import find_bbox, find_image_std, import_difmap_model


degree_to_rad = u.deg.to(u.rad)
mas_to_rad = u.mas.to(u.rad)
# Speed of light [cm / s]
c = const.c.cgs.value
k = const.k_B.cgs.value


def get_r(sample, comp_length, n_comp, jitter_first=True, n_jitters=1):
    j = 0
    if jitter_first:
        j += n_jitters
    return [np.hypot(sample[i*comp_length+j], sample[i*comp_length+j+1]) for i in
            range(n_comp)]


def sort_sample_by_r(sample, n_comp, comp_length=4, jitter_first=True,
                     n_jitters=1):
    r = get_r(sample, comp_length, n_comp, jitter_first, n_jitters)
    indices = np.argsort(r)
    # Construct re-labelled sample
    j = 0
    if jitter_first:
        j += n_jitters
    else:
        n_jitters = 0
    result = np.hstack([sample[j+i*comp_length: j+(i+1)*comp_length] for i in
                        indices])
    return np.hstack((sample[: n_jitters], result))


def sort_samples_by_r(samples, n_comp, comp_length=4, jitter_first=True,
                      n_jitters=1):
    new_samples = list()
    for sample in samples:
        sorted_sample = sort_sample_by_r(sample, n_comp, comp_length, jitter_first,
                                         n_jitters)
        new_samples.append(sorted_sample)
    return np.atleast_2d(new_samples)


def rj_plot_ncomponents_distribution(posterior_file="posterior_sample.txt",
                                     picture_fn=None, jitter_first=True,
                                     n_jitters=1, skip_hyperparameters=False,
                                     type="cg", normed=False, show=True):
    samples = np.atleast_2d(np.loadtxt(posterior_file))
    if type == "cg":
        nn = 0
    elif type == "eg":
        nn = 2
    if jitter_first:
        ncomp_index = 6 + nn + n_jitters
    else:
        ncomp_index = 6 + nn
    if skip_hyperparameters:
        ncomp_index -= (4 + nn)
    values, counts = np.unique(samples[:, ncomp_index], return_counts=True)
    fig, axes = plt.subplots(1, 1)
    if normed:
        counts = counts/sum(counts)
    axes.vlines(values, 0, counts, color='C0', lw=8)
    axes.set_ylim(0, max(counts) * 1.06)
    axes.set_xticks(np.arange(min(values), max(values)+1))
    axes.set_xlabel("# components")
    if normed:
        axes.set_ylabel("P")
    else:
        axes.set_ylabel("N")
    if picture_fn is not None:
        fig.savefig(picture_fn, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def get_samples_for_each_n(samples, jitter_first=True, n_jitters=1, n_max=30,
                           skip_hyperparameters=False,
                           type="cg"):
    samples = np.atleast_2d(samples)
    if type == "cg":
        nn = 0
        comp_length = 4
    elif type == "eg":
        nn = 2
        comp_length = 6
    j = 0
    if jitter_first:
        j += n_jitters
    # dim, max num components, 4 hyperparameters + num components
    j += (7 + nn)
    if skip_hyperparameters:
        j -= (4 + nn)
    n_components = samples[:, j-1]
    out_samples = dict()
    for n in np.array(np.unique(n_components), dtype=int):
        samples_with_n_components = list()
        for sample in samples[n_components == n]:
            one_post_point = list()
            if jitter_first:
                one_post_point.extend(sample[:n_jitters])
            for k in range(n):
                one_post_point.extend([sample[j+k+i*n_max] for i in range(comp_length)])
            samples_with_n_components.append(one_post_point)
        out_samples.update({n: np.atleast_2d(samples_with_n_components)})
    return out_samples


def plot_size_distance_posterior(samples, savefn=None, s=0.6,
                                 sorted_components=False, type="cg"):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = cycle(colors)

    if type == "cg":
        comp_length = 4
    elif type == "eg":
        comp_length = 6

    fig, axes = plt.subplots(1, 1)
    sizes = dict()
    rs = dict()
    n_comps = int(len(samples[0])/comp_length)

    for i_comp in range(n_comps):
        rs[i_comp] = np.hypot(samples[:, 0+i_comp*comp_length], samples[:, 1+i_comp*comp_length])
        sizes[i_comp] = np.exp(samples[:, 3+i_comp*comp_length])

    for i_comp, color in zip(range(n_comps), colors):
        if not sorted_components:
            color = "gray"
        axes.scatter(rs[i_comp], sizes[i_comp], s=s, color=color)

    axes.set_xlabel("r [mas]")
    axes.set_ylabel("FWHM [mas]")
    axes.set_xscale('log')
    axes.set_yscale('log')

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def tb_comp(flux, bmaj, freq, z=0., bmin=None, D=1.):
    """
    :param flux:
        Flux in Jy.
    :param freq:
        Frequency in GHz.
    """
    bmaj *= mas_to_rad
    if bmin is None:
        bmin = bmaj
    else:
        bmin *= mas_to_rad
    freq *= 10**9
    flux *= 10**(-23)
    return 2.*np.log(2)*(1.+z)*flux*c**2/(freq**2*np.pi*k*bmaj*bmin*D)


def plot_tb_distance_posterior(samples, freq_ghz, z=0.0, type="cg", savefn=None, s=0.6,
                               sorted_compnents=False):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = cycle(colors)
    fig, axes = plt.subplots(1, 1)
    tbs = dict()
    rs = dict()
    bmins = dict()
    if type == "cg":
        comp_length = 4
    elif type == "eg":
        comp_length = 6
    else:
        raise Exception("keyword parameter ``type`` should be eg or cg!")
    n_comps = int(len(samples[0])/comp_length)

    for i_comp in range(n_comps):
        rs[i_comp] = np.hypot(samples[:, 0+i_comp*comp_length], samples[:, 1+i_comp*comp_length])
        if type == "eg":
            bmins[i_comp] = samples[:, 4+i_comp*comp_length]
        else:
            bmins[i_comp] = None
        tbs[i_comp] = tb_comp(np.exp(samples[:, 2+i_comp*comp_length]),
                              np.exp(samples[:, 3+i_comp*comp_length]),
                              freq_ghz, z=z,
                              bmin=bmins[i_comp])

    for i_comp, color in zip(range(n_comps), colors):
        if not sorted_compnents:
            color = "gray"
        axes.scatter(rs[i_comp], tbs[i_comp], s=s, color=color)

    # Need manually set ylim because of matplotlib bug
    lg_tb_min = np.floor((np.log10(np.min([tbs[i] for i in range(n_comps)]))))
    lg_tb_max = np.ceil(np.log10(np.max([tbs[i] for i in range(n_comps)])))
    axes.set_ylim([10**lg_tb_min, 10**lg_tb_max])
    axes.set_xlabel("r [mas]")
    axes.set_ylabel("Tb [K]")
    axes.set_xscale('log')
    axes.set_yscale('log')

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


# FIXME: Implement elliptical gaussians
def plot_model_predictions(posterior_file, data_file, rj=True, n_jitters=0, style="reim", savefname=None, show=True,
                           n_samples_to_plot=100, component_type="cg", alpha_model=0.03, skip_hyperparameters=False):

    allowed_component_types = ("cg", "eg", "sphere", "delta")
    if component_type not in allowed_component_types:
        raise Exception(f"Allowed component types are: {allowed_component_types}")
    if style not in ("reim", "ap"):
        raise Exception("Only reim or ap style plotting is supported!")


    function_dict = {"cg": gaussian_circ_ft,
                     "eg": gaussian_ell_ft,
                     "sphere": optically_thin_sphere_ft,
                     "delta": point_ft}

    data_df = pd.read_csv(data_file)
    uv = data_df[["u", "v"]].values
    r = np.hypot(uv[:, 0], uv[:, 1])/10**6
    # dr = np.max(r) - np.min(r)
    # rr = np.linspace(0, np.max(r)+0.1*dr, 1000)
    df = data_df.copy()

    # Zero visibilities
    df["vis_re"] = 0
    df["vis_im"] = 0

    post_samples = np.atleast_2d(np.loadtxt(posterior_file))

    n_post, len_of_sample = post_samples.shape
    print("Posterior shape : ", n_post, len_of_sample)

    samples = post_samples[np.random.randint(0, n_post, n_samples_to_plot)]

    samples_predictions_re = list()
    samples_predictions_im = list()

    if n_jitters > 0:
        jitter_first = True
    else:
        jitter_first = False


    if rj:

        n_max = int(post_samples[0, n_jitters + 1])
        comp_length = int(post_samples[0, n_jitters])

        print(n_max, comp_length)

        if component_type in ("cg", "shpere"):
            assert comp_length == 4
        elif component_type == "eg":
            assert comp_length == 6
        else:
            assert comp_length == 3

        samples_for_each_n = get_samples_for_each_n(samples, jitter_first,
                                                    n_jitters=n_jitters, n_max=n_max,
                                                    skip_hyperparameters=skip_hyperparameters,
                                                    type=component_type)

        for n_local, samples_local in samples_for_each_n.items():
            for sample in samples_local:
                sample_prediction_re = np.zeros(len(df))
                sample_prediction_im = np.zeros(len(df))
                for i_comp in range(n_local):
                    params = sample[n_jitters+i_comp*comp_length: n_jitters + (i_comp + 1)*comp_length]
                    # log of size (except the point)
                    try:
                        params[3] = np.exp(params[3])
                    except IndexError:
                        pass
                    # Using bmin n sampling
                    if component_type == "eg":
                        params[3] = params[3]/params[4]

                    re, im = function_dict[component_type](uv, *params)
                    sample_prediction_re += re
                    sample_prediction_im += im
                samples_predictions_re.append(sample_prediction_re)
                samples_predictions_im.append(sample_prediction_im)


    else:

        # TODO
        raise NotImplementedError


    # Original data plot
    fig = radplot(data_df, style=style, show=False, color="black")
    axes = fig.get_axes()
    for sample_prediction_re, sample_prediction_im in zip(samples_predictions_re, samples_predictions_im):
        if style == "reim":
            axes[0].plot(r, sample_prediction_re, "_", alpha=alpha_model, color="red", ms=5, mew=1.)
            axes[1].plot(r, sample_prediction_im, "_", alpha=alpha_model, color="red", ms=5, mew=1.)
        else:
            axes[0].plot(r, np.hypot(sample_prediction_re, sample_prediction_im), "_", alpha=alpha_model, color="red", ms=5, mew=1.)
            axes[1].plot(r, np.arctan2(sample_prediction_im, sample_prediction_re), "_", alpha=alpha_model, color="red", ms=5, mew=1.)
    if savefname is not None:
        fig.savefig(savefname, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    return fig


def export_difmap_model(comps, out_fname, freq_hz):
    """
    :param comps:
        Iterable of tuples with (flux, x, y, bmaj).
    :param out_fname:
        Path for saving file.
    """
    with open(out_fname, "w") as fo:
        fo.write("! Flux (Jy) Radius (mas)  Theta (deg)  Major (mas)  Axial ratio   Phi (deg) T\n\
! Freq (Hz)     SpecIndex\n")
        for comp in comps:
            if len(comp) == 4:
                # Jy, mas, mas, mas
                flux, x, y, bmaj = comp
                e = "1.00000"
                bpa = "000.000"
                type = "1"
                bmaj = "{:.7f}v".format(bmaj)
            elif len(comp) == 6:
                # Jy, mas, mas, mas, -, deg
                flux, x, y, bmaj, e, bpa = comp
                e = "{}v".format(e)
                bpa = "{}v".format((bpa-np.pi/2)/degree_to_rad)
                bmaj = "{}v".format(bmaj)
                type = "1"
            elif len(comp) == 3:
                flux, x, y = comp
                e = "1.00000"
                bmaj = "0.0000"
                bpa = "000.000"
                type = "0"
            else:
                raise Exception
            # mas
            r = np.hypot(x, y)
            # rad
            theta = np.arctan2(x, y)
            theta /= degree_to_rad
            fo.write("{:>11.7f}v {:>13.7f}v {:>13.5f}v {:>13} {:>13} {:>13} {:>3} {:>12.5e} {:>12d}\n".format(flux, r, theta,
                                                                                                              bmaj, e, bpa, type,
                                                                                                              freq_hz, 0))


def convert_sample_to_difmap_model(sample, out_fname, freq_ghz, type="cg"):
    if type == "cg":
        comp_length = 4
    elif type == "eg":
        comp_length = 6
    else:
        raise Exception("Converting sample to difmap model works for Gaussians only")
    n_comps = int(len(sample)/comp_length)
    components = list()
    for i in range(n_comps):
        subsample = sample[comp_length*i:comp_length*i+comp_length]
        print(subsample)
        flux = subsample[2]
        bmaj = np.exp(subsample[3])
        x = subsample[0]
        y = subsample[1]
        if type == "eg":
            e = subsample[4]
            bpa = subsample[5]
            comp = (flux, x, y, bmaj, e, bpa)
        elif type == "cg":
            comp = (flux, x, y, bmaj)
        components.append(comp)
    components = sorted(components, key=lambda comp: comp[0], reverse=True)
    export_difmap_model(components, out_fname, 1e9*freq_ghz)
    return components


# FIXME: Implement elliptical gaussians
def plot_position_posterior(samples, savefn=None, ra_lim=(-10, 10),
                            dec_lim=(-10, 10), difmap_model_fn=None,
                            n_relative_posterior=None, s=0.6,
                            n_relative_difmap=None, type="cg", figsize=None,
                            sorted_componets=False,
                            fig=None,
                            inverse_xaxis=True,
                            alpha_opacity=0.01):

    """

    :param samples:
        Already sorted posterior samples. For nonRJ gain samples you should
        do first ``postprocess_labebels_gains.sort_samples_by_r`` than feed
        only non-jitter part of the sorted sampler here.
    :param savefn:
    :param ra_lim:
    :param dec_lim:
    :param difmap_model_fn:
    :param n_relative_posterior:
        Number of component in already sorted posterior that should be in phase
        center.
    :param n_relative_difmap:
        Number of component in difmap model that should be in phase center.
    :return:
    """
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = cycle(colors)
    if fig is None:
        fig, axes = plt.subplots(1, 1, figsize=figsize)
    else:
        axes = fig.gca()
    xs = dict()
    ys = dict()
    fluxes = dict()


    if type == "cg":
        comp_length = 4
    elif type == "eg":
        comp_length = 6
    n_comps = int(len(samples[0])/comp_length)
    print("# comp = ", n_comps)


    if n_relative_posterior is not None:
        shift_x = samples[:, 0+n_relative_posterior*comp_length]
        shift_y = samples[:, 1+n_relative_posterior*comp_length]
    else:
        shift_x = 0
        shift_y = 0

    for i_comp in range(n_comps):
        xs[i_comp] = samples[:, 0+i_comp*comp_length] - shift_x
        ys[i_comp] = samples[:, 1+i_comp*comp_length] - shift_y
        fluxes[i_comp] = samples[:, 2+i_comp*comp_length]

    for i_comp, color in zip(range(n_comps), colors):
        if not sorted_componets:
            color = "red"
        axes.scatter(xs[i_comp], ys[i_comp], s=s, color=color, edgecolors='none', alpha=alpha_opacity)

    if difmap_model_fn is not None:
        comps = import_difmap_model(difmap_model_fn)
        print("Read {} components from {}".format(len(comps), difmap_model_fn))
        if n_relative_difmap is not None:
            c7comp = comps[n_relative_difmap]
            shift_x = c7comp.p[1]
            shift_y = c7comp.p[2]
        else:
            shift_x = 0
            shift_y = 0

        for comp in comps:
            if comp.size == 3:
                axes.scatter(-(comp.p[1]-shift_x), -(comp.p[2]-shift_y), s=80, color="black", alpha=1, marker="x")
            elif comp.size == 4:
                # FIXME: Here C7 is putted in the phase center
                e = Circle((-(comp.p[1]-shift_x), -(comp.p[2]-shift_y)), comp.p[3],
                           edgecolor="black", facecolor="red",
                           alpha=0.05)
                axes.add_patch(e)
            elif comp.size == 6:
                raise Exception("Not implemented for elliptical component!")

    axes.set_xlim(ra_lim)
    axes.set_ylim(dec_lim)
    axes.set_xlabel("RA, mas")
    axes.set_ylabel("DEC, mas")
    if inverse_xaxis:
        axes.invert_xaxis()
    axes.set_aspect("equal")

    if savefn is not None:
        fig.savefig(savefn, dpi=600, bbox_inches="tight")
    return fig


def plot_position_posterior_clustered(samples, savefn=None, ra_lim=(-10, 10),
                            dec_lim=(-10, 10), difmap_model_fn=None,
                            n_relative_posterior=None, s=0.6,
                            n_relative_difmap=None, type="cg", figsize=None,
                            sorted_componets=False,
                            fig=None,
                            inverse_xaxis=True,
                            alpha_opacity=0.01,
                            cluster_membership=None):

    """

    :param samples:
        Already sorted posterior samples. For nonRJ gain samples you should
        do first ``postprocess_labebels_gains.sort_samples_by_r`` than feed
        only non-jitter part of the sorted sampler here.
    :param savefn:
    :param ra_lim:
    :param dec_lim:
    :param difmap_model_fn:
    :param n_relative_posterior:
        Number of component in already sorted posterior that should be in phase
        center.
    :param n_relative_difmap:
        Number of component in difmap model that should be in phase center.
    :return:
    """
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = cycle(colors)
    if fig is None:
        fig, axes = plt.subplots(1, 1, figsize=figsize)
    else:
        axes = fig.gca()
    xs = dict()
    ys = dict()
    fluxes = dict()

    if type == "cg":
        comp_length = 4
    elif type == "eg":
        comp_length = 6
    n_comps = int(len(samples[0])/comp_length)
    print("# comp = ", n_comps)

    components = list()
    components_cluster_dict = dict()
    for i in np.unique(cluster_membership):
        components_cluster_dict[i] = list()
    for sample in samples:
        # print(f"Splitting sample {sample} into components")
        splitted = np.split(sample, n_comps)
        for comp in splitted:
            # print(f"Appending component {comp}")
            components.append(comp[:2])


    for comp, cluster_id in zip(components, cluster_membership):
        components_cluster_dict[cluster_id].append(comp)

    for i, color in zip(sorted(np.unique(cluster_membership)), colors):
        xs = [comp[0] for comp in components_cluster_dict[i]]
        ys = [comp[1] for comp in components_cluster_dict[i]]
        axes.scatter(xs, ys, s=s, color=color, edgecolors='none', alpha=alpha_opacity)

    axes.set_xlim(ra_lim)
    axes.set_ylim(dec_lim)
    axes.set_xlabel("RA, mas")
    axes.set_ylabel("DEC, mas")
    if inverse_xaxis:
        axes.invert_xaxis()
    axes.set_aspect("equal")

    if savefn is not None:
        fig.savefig(savefn, dpi=600, bbox_inches="tight")
    return fig


def plot_per_antenna_jitters_and_offsets(samples, uvfits=None, n_antennas=10, save_dir=None, save_basename=None,
                                         plot_jitters=True, plot_offsets=True):
    if uvfits is not None:
        obs = eh.obsdata.load_uvfits(uvfits)
        print(obs.tkey)
        tkey = {i: j for (j, i) in obs.tkey.items()}
    else:
        tkey = {i: str(i) for i in range(n_antennas)}

    if save_basename is None:
        save_basename = ""
    else:
        save_basename = f"{save_basename}_"

    if save_dir is None:
        save_dir = os.getcwd()

    samples = np.atleast_2d(samples)

    # Jitters
    if plot_jitters:
        data = [samples[:, i] for i in range(n_antennas)]
        labels = [tkey[i] for i in range(n_antennas)]
        df = pd.DataFrame.from_dict(OrderedDict(zip(labels, data)))
        fig, axes = plt.subplots(1, 1)
        axes = sns.boxplot(data=df, orient='h', ax=axes)
        axes.set_xlabel(r"$\log{\sigma_{\rm ant}}$")
        axes.set_ylabel("Antenna")
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f"{save_basename}jitters.png"), bbox_inches="tight")
        plt.show()

    # Offsets
    if plot_offsets:
        data = [samples[:, i+n_antennas] for i in range(n_antennas)]
        labels = [tkey[i] for i in range(n_antennas)]
        df = pd.DataFrame.from_dict(OrderedDict(zip(labels, data)))
        fig, axes = plt.subplots(1, 1)
        axes = sns.boxplot(data=df, orient='h', ax=axes)
        axes.set_xlabel(r"Offset")
        axes.set_ylabel("Antenna")
        axes.axvline(1., lw=1, color="k")
        plt.tight_layout()
        fig.savefig(os.path.join(save_dir, f"{save_basename}offsets.png"), bbox_inches="tight")
        plt.show()

    return axes


def get_points_for_clustering_from_samples(samples, n_comp, standardize=False):
    from sklearn.preprocessing import StandardScaler
    components = list()
    for sample in samples:
        # print(f"Splitting sample {sample} into components")
        splitted = np.split(sample, n_comp)
        for comp in splitted:
            # print(f"Appending component {comp}")
            components.append(comp[:2])
    components = np.atleast_2d(components)
    print("========================================================================")
    print(f"Number of components in current model: {n_comp}, # of posterior samples : {len(samples)}")
    if standardize:
        components = StandardScaler().fit_transform(components)
    return components


def clusterization(posterior_file, jitter_first=True, n_jitters=1, n_max=30,
                   skip_hyperparameters=False, component_type="cg",
                   algorithm="hdbscan",
                   dbscan_min_core_samples_frac_of_posterior_size=0.5,
                   dbscan_eps=0.1,
                   hdbscan_min_cluster_frac_of_posterior_size=0.5):

    # if type == "cg":
    #     comp_length = 4
    # elif type == "eg":
    #     comp_length = 6
    # else:
    #     raise Exception

    posterior_samples = np.loadtxt(posterior_file)
    samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
                                                n_jitters=n_jitters, n_max=n_max,
                                                skip_hyperparameters=skip_hyperparameters,
                                                type=component_type)

    labels_dict = dict()
    for n_comp, samples in samples_for_each_n.items():
        if jitter_first and n_jitters > 0:
            samples = samples[:, n_jitters:]

        X = get_points_for_clustering_from_samples(samples, n_comp)

        print("========================================================================")
        print(f"Number of components in current model: {n_comp}, # of posterior samples : {len(samples)}")
        # The number of samples (or total weight) in a neighborhood for a point
        # to be considered as a core point. It is set as a fraction of posterior samples per component. 
        print("number of points to cluster : ", len(X))
        min_core_samples = int(dbscan_min_core_samples_frac_of_posterior_size*len(X)/n_comp)
        if min_core_samples == 0:
            min_core_samples = 1
        print("min core # of points : ", min_core_samples)
        if algorithm == "dbscan":
            db = DBSCAN(eps=dbscan_eps, min_samples=min_core_samples).fit(X)
        elif algorithm == "hdbscan":
            min_cluster_size = int(hdbscan_min_cluster_frac_of_posterior_size*len(samples))
            if min_cluster_size < 2:
                min_cluster_size = 2
            db = HDBSCAN(min_cluster_size=min_cluster_size).fit(X)
        else:
            raise Exception

        labels = db.labels_

        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)

        print("Estimated number of clusters: {}".format(n_clusters_))
        print("Estimated percent of noise points: {:.2f}".format(100*n_noise_/len(X)))

        labels_dict[n_comp] = db.labels_

    return labels_dict


def cluster_optics(samples):
    clust = OPTICS(min_samples=50, xi=0.05, min_cluster_size=0.05)
    X = samples.copy()
    clust.fit(X)
    space = np.arange(len(X))
    reachability = clust.reachability_[clust.ordering_]
    labels = clust.labels_[clust.ordering_]
    fig, ax1 = plt.subplots(1, 1)
    # Reachability plot
    colors = ["g.", "r.", "b.", "y.", "c."]
    for klass, color in zip(range(0, 5), colors):
        Xk = space[labels == klass]
        Rk = reachability[labels == klass]
        ax1.plot(Xk, Rk, color, alpha=0.3)
    ax1.plot(space[labels == -1], reachability[labels == -1], "k.", alpha=0.3)
    ax1.plot(space, np.full_like(space, 2.0, dtype=float), "k-", alpha=0.5)
    ax1.plot(space, np.full_like(space, 0.5, dtype=float), "k-.", alpha=0.5)
    ax1.set_ylabel("Reachability (epsilon distance)")
    ax1.set_title("Reachability Plot")
    plt.show()


def plot(X, labels, probabilities=None, parameters=None, ground_truth=False, ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 4))
    labels = labels if labels is not None else np.ones(X.shape[0])
    probabilities = probabilities if probabilities is not None else np.ones(X.shape[0])
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    # The probability of a point belonging to its labeled cluster determines
    # the size of its marker
    proba_map = {idx: probabilities[idx] for idx in range(len(labels))}
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_index = np.where(labels == k)[0]
        for ci in class_index:
            ax.plot(
                X[ci, 0],
                X[ci, 1],
                "x" if k == -1 else "o",
                markerfacecolor=tuple(col),
                markeredgecolor="k",
                markersize=4 if k == -1 else 1 + 5 * proba_map[ci],
            )
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    preamble = "True" if ground_truth else "Estimated"
    title = f"{preamble} number of clusters: {n_clusters_}"
    if parameters is not None:
        parameters_str = ", ".join(f"{k}={v}" for k, v in parameters.items())
        title += f" | {parameters_str}"
    ax.set_title(title)
    plt.tight_layout()


def cluster_hdbscan(samples, n_comp):
    X = get_points_for_clustering_from_samples(samples, n_comp, standardize=True)
    # First, tune ``min_cluster_size``, then - ``min_samples``
    hdb = HDBSCAN(min_cluster_size=int(0.5*len(samples))).fit(X)
    plot(X, hdb.labels_, hdb.probabilities_)


def plot_corner(samples, n_comps, cluster_membership=None, savefig=None, component_type="cg"):
    """
    :param samples:
        Samples w/o jitter of type returned by ``get_samples_for_each_n`` function.
    :param n_comps:
        Number of components in a model.
    :param cluster_membership: (optional)
        labels_dict[n_comp]. If ``None`` then assume that samples are sorted somehow.
        (default: ``None``)
    """
    # Construct labels
    labels = list()
    if component_type == "cg" or component_type == "sphere":
        comp_length = 4
        for i in range(n_comps):
            labels.extend([r"RA$_{}$".format(i+1),
                           r"DEC$_{}$".format(i+1),
                           r"$S_{}$".format(i+1),
                           r"$\theta_{}$".format(i+1)])
    elif component_type == "eg":
        comp_length = 6
        for i in range(n_comps):
            labels.extend([r"RA$_{}$".format(i+1),
                           r"DEC$_{}$".format(i+1),
                           r"$S_{}$".format(i+1),
                           r"$\bmaj_{}$".format(i+1),
                           r"$\e_{}$".format(i+1),
                           r"$\bpa_{}$".format(i+1)])
    else:
        raise Exception("component_type must be cg, eg or sphere!")

    # Cluster components
    if cluster_membership is not None:
        components = list()
        components_cluster_dict = dict()
        unique_cluster_ids = list(np.unique(cluster_membership))
        for cluster_id in unique_cluster_ids:
            components_cluster_dict[cluster_id] = list()
        for sample in samples:
            # print(f"Splitting sample {sample} into components")
            splitted = np.split(sample, n_comps)
            for comp in splitted:
                # print(f"Appending component {comp}")
                components.append(comp)

        for comp, cluster_id in zip(components, cluster_membership):
            components_cluster_dict[cluster_id].append(comp)

        # Drop non-cluster points
        try:
            unique_cluster_ids.remove(-1)
            components_cluster_dict.pop(-1)
        # If all samples lie in clusters
        except ValueError:
            pass

        for cluster_id in unique_cluster_ids:
            components_cluster_dict[cluster_id] = np.atleast_2d(components_cluster_dict[cluster_id])

        # Convert size from log scale
        for cluster_id in unique_cluster_ids:
            components_cluster_dict[cluster_id][:, 3] = np.exp(components_cluster_dict[cluster_id][:, 3]) 

        # Some components are in -1 class (no cluster). Trim other clusters to the minimal common size.
        minimal_common_len = min([components_cluster_dict[cluster_id].shape[0] for cluster_id in unique_cluster_ids])
        for cluster_id in unique_cluster_ids:
            components_cluster_dict[cluster_id] = components_cluster_dict[cluster_id][:minimal_common_len, :]

        median_fluxes_dict = dict()
        for cluster_id in unique_cluster_ids:
            median_fluxes_dict[cluster_id] = np.median(components_cluster_dict[cluster_id][:, 2])
        # Sort cluster IDs in decreasing of the median flux
        clusters_id_flux_sorted = sorted(median_fluxes_dict, key=lambda x: median_fluxes_dict[x], reverse=True)
        # Concatenate in the same order in a single 2D array
        new_samples = np.concatenate([components_cluster_dict[i] for i in clusters_id_flux_sorted], axis=1)

    else:
        # Just sorted samples
        new_samples = samples
        # Convert sizes from log-scale
        for i in range(n_comps):
            new_samples[:, 3 + i*comp_length] = np.exp(new_samples[:, 3 + i*comp_length])


    fig = corner(new_samples, show_titles=True, title_fmt=".3f", quantiles=[0.16, 0.50, 0.84], labels=labels)
    if savefig is not None:
        fig.savefig(savefig, bbox_inches="tight", dpi=100)
    plt.close()


def postprocess_run(save_basename, posterior_file, data_file, original_ccfits, save_dir,
                    n_max, has_jitter, component_type, skip_hp=True, pixsize_mas=None,
                    plot_type="reim", corner_with_clustered=True):
    save_rj_ncomp_distribution_file = os.path.join(save_dir, f"{save_basename}_ncomponents_distribution.png")
    n_jitters = 1
    # Plot all samples - for easy handling component cluster membership
    n_max_samples_to_plot = 1000
    jitter_first = has_jitter
    if component_type == "cg":
        comp_length = 4
    elif component_type == "eg":
        comp_length = 6
    else:
        raise Exception
    posterior_samples = np.loadtxt(posterior_file)
    import matplotlib
    matplotlib.use("Agg")


    fig = rj_plot_ncomponents_distribution(posterior_file, picture_fn=save_rj_ncomp_distribution_file,
                                           jitter_first=jitter_first, n_jitters=n_jitters, type=component_type,
                                           normed=True, show=False, skip_hyperparameters=skip_hp)

    labels_dict = clusterization(posterior_file, jitter_first=jitter_first, n_jitters=n_jitters, n_max=n_max,
                                 skip_hyperparameters=skip_hp, component_type=component_type,
                                 algorithm="hdbscan",
                                 dbscan_min_core_samples_frac_of_posterior_size=0.3,
                                 dbscan_eps=0.2,
                                 hdbscan_min_cluster_frac_of_posterior_size=0.5)

    plot_model_predictions(posterior_file, data_file, rj=True, n_jitters=n_jitters, component_type=component_type,
                           style=plot_type, n_samples_to_plot=1000, alpha_model=0.01, skip_hyperparameters=skip_hp,
                           savefname=os.path.join(save_dir, f"{save_basename}_radplot.png"))


    samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
                                                n_jitters=n_jitters, n_max=n_max,
                                                skip_hyperparameters=skip_hp,
                                                type=component_type)
    n_components_spread = samples_for_each_n.keys()

    ccimage = create_clean_image_from_fits_file(original_ccfits)
    beam = ccimage.beam
    if pixsize_mas is not None:
        npixels_beam = int(np.pi*beam[0]*beam[1]/(4.*np.log(2)*pixsize_mas**2))
        print(f"Estimated beam area (pixels) = {npixels_beam}")
    else:
        npixels_beam = 30

    std = find_image_std(ccimage.image, npixels_beam, min_num_pixels_used_to_estimate_std=100,
                         blc=None, trc=None)
    blc, trc = find_bbox(ccimage.image, level=4*std, min_maxintensity_jyperbeam=10*std,
                         min_area_pix=3*npixels_beam, delta=5)
    fig = iplot(ccimage.image, x=ccimage.x, y=ccimage.y,
                min_abs_level=3*std, beam=(beam[0], beam[1], np.rad2deg(beam[2])), show_beam=True, blc=blc, trc=trc,
                components=None, close=False, plot_colorbar=False, show=False,
                contour_linewidth=0.25, contour_color='k')
    fig.savefig(os.path.join(save_dir, f"{save_basename}_CLEAN_image.png"), dpi=600)
    for n_component in n_components_spread:
        # Remove jitters & offsets
        samples_to_plot = samples_for_each_n[n_component][:, n_jitters:]
        # Sort samples by r
        sorted_samples_to_plot = sort_samples_by_r(samples_to_plot, n_component, comp_length=comp_length, jitter_first=False)
        n_samples = len(samples_to_plot)
        if n_samples > n_max_samples_to_plot:
            n_samples = n_max_samples_to_plot
        fig_p = pickle.loads(pickle.dumps(fig))
        fig_p1 = pickle.loads(pickle.dumps(fig))
        fig_p2 = pickle.loads(pickle.dumps(fig))

        # Vanilla
        fig_out = plot_position_posterior(samples_to_plot[:n_max_samples_to_plot, :],
                                          savefn=None, ra_lim=None, dec_lim=None,
                                          difmap_model_fn=None, type=component_type, s=3.0, figsize=None,
                                          sorted_componets=False, fig=fig_p,
                                          inverse_xaxis=False,
                                          alpha_opacity=0.3)
        fig_out.savefig(os.path.join(save_dir, f"{save_basename}_CLEAN_image_ncomp_{n_component}.png"), dpi=600)
        plt.close(fig_out)

        # Sorted samples
        fig_out = plot_position_posterior(sorted_samples_to_plot[:n_max_samples_to_plot, :],
                                          savefn=None, ra_lim=None, dec_lim=None,
                                          difmap_model_fn=None, type=component_type, s=3.0, figsize=None,
                                          sorted_componets=True, fig=fig_p1,
                                          inverse_xaxis=False,
                                          alpha_opacity=0.3)
        fig_out.savefig(os.path.join(save_dir, f"{save_basename}_CLEAN_image_ncomp_{n_component}_sorted.png"), dpi=600)
        plt.close(fig_out)

        fig_out = plot_position_posterior_clustered(samples_to_plot,
                                                    savefn=None, ra_lim=None, dec_lim=None,
                                                    difmap_model_fn=None, type=component_type, s=1.0, figsize=None,
                                                    sorted_componets=False, fig=fig_p2,
                                                    inverse_xaxis=False,
                                                    alpha_opacity=0.3,
                                                    cluster_membership=labels_dict[n_component])
        fig_out.savefig(os.path.join(save_dir, f"{save_basename}_CLEAN_image_ncomp_{n_component}_clusters.png"), dpi=600)
        plt.close(fig_out)

        if corner_with_clustered:
            plot_corner(samples_to_plot, n_component, cluster_membership=labels_dict[n_component],
                        savefig=os.path.join(save_dir, f"{save_basename}_corner.png"))
        else:
            plot_corner(sorted_samples_to_plot, n_component, cluster_membership=None,
                        savefig=os.path.join(save_dir, f"{save_basename}_corner.png"))

        # f = plot_size_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :],
        #                                  savefn=os.path.join(save_dir, f"r_R_ncomp_{n_component}.png"),
        #                                  type=component_type)
        # plt.close(f)
        #
        # f = plot_tb_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :], freq_ghz, type=component_type,
        #                                savefn=os.path.join(save_dir, f"r_Tb_ncomp_{n_component}.png"))
        # plt.close(f)


if __name__ == "__main__":
    pixsize_band_dict = {"Q": 0.03, "U": 0.1, "X": 0.2, "C": 0.3, "S": 0.5}
    data_files_dir = "/home/ilya/data/VLBI_Gaia/3comp"
    data_files = glob.glob(os.path.join(data_files_dir, "*.csv"))
    print("Found data files: ", data_files)
    for data_file in data_files:
        if data_file != "/home/ilya/data/VLBI_Gaia/3comp/J0000-3221_S_2017_01_16_pet_vis.csv":
            continue
        # data_file = "/home/ilya/data/VLBI_Gaia/2comp/J0823-0939_X_2017_02_24_pus_vis.csv"
        data_fn = os.path.split(data_file)[-1]
        source, band, year, month, day, author, product = data_fn.split("_")
        product = product.split(".")[0]
        assert product == "vis"
        save_basename = f"{source}_{band}_{year}_{month}_{day}"
        results_dir = os.path.join(data_files_dir, save_basename)
        posterior_file = os.path.join(results_dir, f"posterior_sample_{save_basename}.txt")
        if not os.path.exists(posterior_file):
            raise Exception(f"No posterior file : {posterior_file}")
        original_ccfits = os.path.join(data_files_dir, f"{save_basename}_{author}_map.fits")
        if not os.path.exists(original_ccfits):
            raise Exception(f"No CCFITS image : {original_ccfits}")


        save_dir = results_dir
        # pathlib.Path(save_dir).mkdir(parents=True, exist_ok=True)
            
        n_max = 3 
        has_jitter = True
        component_type = "cg"
        pixsize_mas = pixsize_band_dict[band]
        plot_type = "reim"

        postprocess_run(save_basename, posterior_file, data_file, original_ccfits, save_dir,
                        n_max, has_jitter, component_type, skip_hp=True, pixsize_mas=pixsize_mas,
                        plot_type=plot_type, corner_with_clustered=False)




    sys.exit(0)
    data_file = "/home/ilya/data/VLBI_Gaia/2comp/J0823-0939_X_2017_02_24_pus_vis.csv"
    df = pd.read_csv(data_file)
    posterior_file = "/home/ilya/data/VLBI_Gaia/2comp/posterior_sample_J0823-0939_X_2017_02_24.txt"
    # save_dir = "/home/ilya/data/rjbam/0212+735/2019_08_15/jitters_offsets/circular_2Dprior/"
    # save_dir = "/home/ilya/data/rjbam/1502+106/4comp_jitters_2D"
    save_dir = "/home/ilya/data/VLBI_Gaia/2comp"
    save_rj_ncomp_distribution_file = os.path.join(save_dir, "ncomponents_distribution.png")
    original_ccfits = "/home/ilya/data/VLBI_Gaia/2comp/J0823-0939_X_2017_02_24_pus_map.fits"
    n_max = 10 
    n_jitters = 1
    # Plot all samples - for easy handling component cluster membership
    n_max_samples_to_plot = 1000
    jitter_first = True
    skip_hp = True
    component_type = "cg"
    if component_type == "cg":
        comp_length = 4
    elif component_type == "eg":
        comp_length = 6
    else:
        raise Exception
    pixsize_mas = 0.3
    freq_ghz = 4.6
    posterior_samples = np.loadtxt(posterior_file)
    import matplotlib
    matplotlib.use("TkAgg")
    save_basename = os.path.split(data_file)[-1].split(".csv")[0]


    # plot_per_antenna_jitters_and_offsets(posterior_samples, uvfits=uvfits, save_dir=save_dir,
    #                                      save_basename=save_basename, plot_offsets=False,
    #                                      n_antennas=n_antennas)

    fig = rj_plot_ncomponents_distribution(posterior_file, picture_fn=save_rj_ncomp_distribution_file,
                                           jitter_first=jitter_first, n_jitters=n_jitters, type=component_type,
                                           normed=True, show=False, skip_hyperparameters=skip_hp)
    # samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
    #                                             n_jitters=n_jitters, n_max=n_max,
    #                                             skip_hyperparameters=skip_hp,
    #                                             type=component_type)

    # sys.exit(0)

    labels_dict = clusterization(posterior_file, jitter_first=jitter_first, n_jitters=n_jitters, n_max=n_max,
                                 skip_hyperparameters=skip_hp, component_type=component_type,
                                 algorithm="hdbscan",
                                 dbscan_min_core_samples_frac_of_posterior_size=0.3,
                                 dbscan_eps=0.2,
                                 hdbscan_min_cluster_frac_of_posterior_size=0.5)

    # sys.exit(0)

    plot_model_predictions(posterior_file, data_file, rj=True, n_jitters=n_jitters, component_type=component_type,
                           style="ap", n_samples_to_plot=1000, alpha_model=0.01, skip_hyperparameters=skip_hp,
                           savefname=os.path.join(save_dir, "radplot.png"))


    # sys.exit(0)

    samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
                                                n_jitters=n_jitters, n_max=n_max,
                                                skip_hyperparameters=skip_hp,
                                                type=component_type)

    n_components_spread = samples_for_each_n.keys()


    ccimage = create_clean_image_from_fits_file(original_ccfits)
    beam = ccimage.beam
    # Number of pixels in beam
    npixels_beam = 100

    std = find_image_std(ccimage.image, npixels_beam, min_num_pixels_used_to_estimate_std=100,
                         blc=None, trc=None)
    blc, trc = find_bbox(ccimage.image, level=4*std, min_maxintensity_jyperbeam=10*std,
                         min_area_pix=3*npixels_beam, delta=30)
    fig = iplot(ccimage.image, x=ccimage.x, y=ccimage.y,
                min_abs_level=3*std, beam=(beam[0], beam[1], np.rad2deg(beam[2])), show_beam=True, blc=blc, trc=trc,
                components=None, close=False, plot_colorbar=False, show=False,
                contour_linewidth=0.25, contour_color='k')
    fig.savefig(os.path.join(save_dir, "CLEAN_image.png"), dpi=600)
    for n_component in n_components_spread:
        # Remove jitters & offsets
        samples_to_plot = samples_for_each_n[n_component][:, n_jitters:]
        # Sort samples by r
        sorted_samples_to_plot = sort_samples_by_r(samples_to_plot, n_component, comp_length=comp_length, jitter_first=False)
        n_samples = len(samples_to_plot)
        if n_samples > n_max_samples_to_plot:
            n_samples = n_max_samples_to_plot
        fig_p = pickle.loads(pickle.dumps(fig))
        fig_p1 = pickle.loads(pickle.dumps(fig))
        fig_p2 = pickle.loads(pickle.dumps(fig))

        # Vanilla
        fig_out = plot_position_posterior(samples_to_plot[:n_max_samples_to_plot, :],
                                          savefn=None, ra_lim=None, dec_lim=None,
                                          difmap_model_fn=None, type=component_type, s=3.0, figsize=None,
                                          sorted_componets=False, fig=fig_p,
                                          inverse_xaxis=False,
                                          alpha_opacity=0.3)
        fig_out.savefig(os.path.join(save_dir, f"CLEAN_image_ncomp_{n_component}.png"), dpi=600)
        plt.close(fig_out)

        # Sorted samples
        fig_out = plot_position_posterior(sorted_samples_to_plot[:n_max_samples_to_plot, :],
                                          savefn=None, ra_lim=None, dec_lim=None,
                                          difmap_model_fn=None, type=component_type, s=3.0, figsize=None,
                                          sorted_componets=True, fig=fig_p1,
                                          inverse_xaxis=False,
                                          alpha_opacity=0.3)
        fig_out.savefig(os.path.join(save_dir, f"CLEAN_image_ncomp_{n_component}_sorted.png"), dpi=600)
        plt.close(fig_out)

        fig_out = plot_position_posterior_clustered(samples_to_plot,
                                                    savefn=None, ra_lim=None, dec_lim=None,
                                                    difmap_model_fn=None, type=component_type, s=1.0, figsize=None,
                                                    sorted_componets=False, fig=fig_p2,
                                                    inverse_xaxis=False,
                                                    alpha_opacity=0.3,
                                                    cluster_membership=labels_dict[n_component])
        fig_out.savefig(os.path.join(save_dir, f"CLEAN_image_ncomp_{n_component}_clusters.png"), dpi=600)
        plt.close(fig_out)

        # f = plot_size_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :],
        #                                  savefn=os.path.join(save_dir, f"r_R_ncomp_{n_component}.png"),
        #                                  type=component_type)
        # plt.close(f)
        #
        # f = plot_tb_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :], freq_ghz, type=component_type,
        #                                savefn=os.path.join(save_dir, f"r_Tb_ncomp_{n_component}.png"))
        # plt.close(f)
