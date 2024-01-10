import os
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from matplotlib.patches import Ellipse, Circle
from itertools import cycle
import seaborn as sns
from astropy import units as u
from astropy import constants as const
import ehtim as eh
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
                           n_samples_to_plot=100, component_type="cg", alpha_model=0.03):

    allowed_component_types = ("cg", "eg", "sphere", "delta")
    if component_type not in allowed_component_types:
        raise Exception(f"Allowed component types are: {allowed_component_types}")
    if style not in ("reim", "ap"):
        raise Exception("Only reim or ap style plotting is supported!")


    function_dict = {"cg": gaussian_circ_ft,
                     "eg": gaussian_ell_ft,
                     "sphere": optically_thin_sphere_ft,
                     "delta": point_ft}

    data_df = pd.read_csv(data_file, names=["t1", "t2", "u", "v", "vis_re", "vis_im", "error"], delim_whitespace=True)
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
        if component_type in ("cg", "shpere"):
            assert comp_length == 4
        elif component_type == "eg":
            assert comp_length == 6
        else:
            assert comp_length == 3

        samples_for_each_n = get_samples_for_each_n(samples, jitter_first,
                                                    n_jitters=n_jitters, n_max=n_max,
                                                    skip_hyperparameters=False,
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


def plot_per_antenna_jitters(samples, uvfits=None, n_antennas=10):
    if uvfits is not None:
        obs = eh.obsdata.load_uvfits(uvfits)
        tkey = {i: j for (j, i) in obs.tkey.items()}
    else:
        tkey = {i: str(i) for i in range(n_antennas)}
    data = [samples[:, i] for i in range(n_antennas)]
    labels = [tkey[i] for i in range(n_antennas)]
    df = pd.DataFrame.from_dict(OrderedDict(zip(labels, data)))
    axes = sns.boxplot(data=df, orient='h')
    axes.set_xlabel(r"$\log{\sigma_{\rm ant}}$")
    axes.set_ylabel("Antenna")
    plt.tight_layout()
    plt.show()
    return axes


def plot_per_antenna_jitters_and_offsets(samples, uvfits=None, n_antennas=10, save_dir=None, save_basename=None):
    if uvfits is not None:
        obs = eh.obsdata.load_uvfits(uvfits)
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


if __name__ == "__main__":

    # uvfits = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2023_05_03.uvf"
    uvfits = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2012_11_11.uvf"
    # uvfits = "/home/ilya/data/rjbam/0212+735/2019_08_15/0212+735.u.2019_08_15.uvf"
    # data_file = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2023_07_01_60sec.txt"
    # df = pd.read_csv(data_file, names=["u", "v", "vis_re", "vis_im", "error"], delim_whitespace=True)
    # data_file = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2023_05_03_60sec_antennas.txt"
    data_file = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2012_11_11_60sec_antennas.txt"
    # data_file = "/home/ilya/data/rjbam/0212+735/2019_08_15/0212+735.u.2019_08_15_60sec_antennas.txt"
    df = pd.read_csv(data_file, names=["t1", "t2", "u", "v", "vis_re", "vis_im", "error"], delim_whitespace=True)
    posterior_file = "/home/ilya/github/bam/posterior_sample.txt"
    # posterior_file = "/home/ilya/github/bam/Release/posterior_sample.txt"
    # old
    # posterior_file = "/home/ilya/github/bam/posterior_sample_rjell.txt"
    # save_dir = "/home/ilya/data/rjbam/0851+202/2023_05_03/jitters_offsets"
    save_dir = "/home/ilya/data/rjbam/0851+202/2012_11_11/jitters_offsets/circular"
    # save_dir = "/home/ilya/data/rjbam/0212+735/2019_08_15/jitters_offsets"
    # save_dir = "/home/ilya/data/rjbam/0212+735/2019_08_15/old"
    save_rj_ncomp_distribution_file = os.path.join(save_dir, "ncomponents_distribution.png")
    # original_ccfits = "/home/ilya/data/rjbam/0851+202/0851+202.u.2023_05_03.icn.fits"
    original_ccfits = "/home/ilya/data/rjbam/0851+202/0851+202.u.2012_11_11.icn.fits"
    # original_ccfits = "/home/ilya/data/rjbam/0212+735/2019_08_15/0212+735.u.2019_08_15.icn.fits"
    n_max = 20
    n_jitters = 20
    n_max_samples_to_plot = 500
    jitter_first = True
    component_type = "eg"
    pixsize_mas = 0.1
    freq_ghz = 15.4
    posterior_samples = np.loadtxt(posterior_file)
    import matplotlib
    matplotlib.use("TkAgg")
    save_basename = os.path.split(uvfits)[-1].split(".uvf")[0]

    plot_model_predictions(posterior_file, data_file, rj=True, n_jitters=n_jitters, component_type="eg",
                           style="reim", n_samples_to_plot=1000, alpha_model=0.01)


    plot_per_antenna_jitters_and_offsets(posterior_samples, uvfits=uvfits, save_dir=save_dir, save_basename=save_basename)

    fig = rj_plot_ncomponents_distribution(posterior_file, picture_fn=save_rj_ncomp_distribution_file,
                                           jitter_first=jitter_first, n_jitters=n_jitters, type=component_type,
                                           normed=True, show=False)
    samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
                                                n_jitters=n_jitters, n_max=n_max,
                                                skip_hyperparameters=False,
                                                type=component_type)

    n_components_spread = samples_for_each_n.keys()


    ccimage = create_clean_image_from_fits_file(original_ccfits)
    beam = ccimage.beam
    # Number of pixels in beam
    npixels_beam = np.pi*beam[0]*beam[1]/(4*np.log(2)*pixsize_mas**2)

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
        samples_to_plot = samples_for_each_n[n_component][:, n_jitters:]
        n_samples = len(samples_to_plot)
        if n_samples > n_max_samples_to_plot:
            n_samples = n_max_samples_to_plot
        fig_p = pickle.loads(pickle.dumps(fig))
        fig_out = plot_position_posterior(samples_to_plot[:n_max_samples_to_plot, :],
                                          savefn=None, ra_lim=None, dec_lim=None,
                                          difmap_model_fn=None, type=component_type, s=1.0, figsize=None,
                                          sorted_componets=False, fig=fig_p,
                                          inverse_xaxis=False,
                                          alpha_opacity=0.03)
        fig_out.savefig(os.path.join(save_dir, f"CLEAN_image_ncomp_{n_component}.png"), dpi=600)
        plt.close(fig_out)
        f = plot_size_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :],
                                         savefn=os.path.join(save_dir, f"r_R_ncomp_{n_component}.png"),
                                         type="eg")
        plt.close(f)
        f = plot_tb_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :], freq_ghz, type="eg",
                                       savefn=os.path.join(save_dir, f"r_Tb_ncomp_{n_component}.png"))
        plt.close(f)
        # f = plot_model_predictions(samples_to_plot[:n_max_samples_to_plot, :], df,
        #                            savefname=os.path.join(save_dir, f"radplot_{n_component}.png"), show=False)
        # plt.close(f)
