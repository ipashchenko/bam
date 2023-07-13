import numpy as np
from astropy import units as u
from astropy import constants as const
import corner
import matplotlib
label_size = 16
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt

mas_to_rad = u.mas.to(u.rad)
degree_to_rad = u.deg.to(u.rad)
# Speed of light [cm / s]
c = const.c.cgs.value
k = const.k_B.cgs.value


def get_r(sample, comp_length):
    """
    Get distances of components from the phase center.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        List of component distances.
    """
    length = len(sample)
    return [np.hypot(sample[i*comp_length+0], sample[i*comp_length+1]) for i in
            range(int(length/comp_length))]


def get_x(sample, comp_length):
    """
    Get x-coordinate of components from the phase center.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        List of component x-coordinates.
    """
    length = len(sample)
    return [sample[i*comp_length+0] for i in range(int(length/comp_length))]


def get_y(sample, comp_length):
    """
    Get y-coordinate of components from the phase center.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        List of component y-coordinates.
    """
    length = len(sample)
    return [sample[i*comp_length+1] for i in range(int(length/comp_length))]


def get_F(sample, comp_length):
    """
    Get fluxes of components.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        List of component fluxes.
    """
    length = len(sample)
    return [(sample[i*comp_length+2]) for i in range(int(length/comp_length))]


def get_size(sample, comp_length):
    """
    Get sizes of components.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        List of component sizes.
    """
    length = len(sample)
    return [(sample[i*comp_length+3]) for i in range(int(length/comp_length))]


def sort_sample_by_r(sample, comp_length=4):
    """
    Sort sample such that component with lowest distance from phase center goes
    first.


    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        Array with sorted sample.
    """
    r = get_r(sample, comp_length)
    indices = np.argsort(r)
    # Construct re-labelled sample
    return np.hstack([sample[i*comp_length: (i+1)*comp_length] for i in
                      indices])


def sort_sample_by_DEC(sample, comp_length=4, inverse=False):
    dec = get_x(sample, comp_length)
    indices = np.argsort(dec)
    if inverse:
        indices = indices[::-1]
    # Construct re-labelled sample
    return np.hstack([sample[i*comp_length: (i+1)*comp_length] for i in
                      indices])


def sort_sample_by_RA(sample, comp_length=4):
    ra = get_y(sample, comp_length)
    indices = np.argsort(ra)
    # Construct re-labelled sample
    return np.hstack([sample[i*comp_length: (i+1)*comp_length] for i in
                      indices])


def sort_sample_by_F(sample, comp_length=4):
    """
    Sort sample such that component with highest flux goes first.

    :param sample:
        Iterable of component properties, i.e. [x1, x2, F1, size1, x2, y2, F2,
        size2, ...].
    :param comp_length:
        Number of component parameters (4 for circular Gaussian).
    :return:
        Array with sorted sample.
    """
    F = get_F(sample, comp_length)
    indices = np.argsort(F)[::-1]
    # Construct re-labelled sample
    return np.hstack([sample[i*comp_length: (i+1)*comp_length] for i in
                      indices])


def sort_samples_by_r(samples, comp_length=4):
    """
    Sort each sample by distance from phase center..
    """
    new_samples = list()
    for sample in samples:
        new_samples.append(sort_sample_by_r(sample, comp_length))
    return np.atleast_2d(new_samples)


def sort_samples_by_dec(samples, comp_length=4, inverse=False):
    """
    Sort each sample by DEC.
    """
    new_samples = list()
    for sample in samples:
        new_samples.append(sort_sample_by_DEC(sample, comp_length, inverse=inverse))
    return np.atleast_2d(new_samples)


def sort_samples_by_ra(samples, comp_length=4):
    """
    Sort each sample by RA.
    """
    new_samples = list()
    for sample in samples:
        new_samples.append(sort_sample_by_RA(sample, comp_length))
    return np.atleast_2d(new_samples)


def sort_samples_by_F(samples, comp_length=4):
    """
    Sort each sample by flux.
    """
    new_samples = list()
    for sample in samples:
        new_samples.append(sort_sample_by_F(sample, comp_length))
    return np.atleast_2d(new_samples)


def tb_comp(flux, bmaj, freq, z=0, bmin=None, D=1):
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


def lg_size_for_given_flux_and_tb(flux_jy, lg_tb, freq_ghz=15.4, z=0.0, D=1.0):
    freq = freq_ghz * 10**9
    flux = flux_jy * 10**(-23)
    lg_bmaj = 0.5*(np.log10(2.*np.log(2)*(1.+z)*c**2/(np.pi*k*freq**2*D))
                   + np.log10(flux) - lg_tb)
    return lg_bmaj - np.log10(mas_to_rad)


def plot_position_posterior(samples, savefn=None, ra_lim=(-10, 10), dec_lim=(-10, 10), n_relative_posterior=None,
                            comp_length=4):
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
    fig, axes = plt.subplots(1, 1)
    xs = dict()
    ys = dict()
    fluxes = dict()
    n_comps = int(len(samples[0])/comp_length)

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
        axes.scatter(xs[i_comp], ys[i_comp], s=0.6, color=color)

    axes.set_xlim(ra_lim)
    axes.set_ylim(dec_lim)
    axes.set_xlabel("RA, mas")
    axes.set_ylabel("DEC, mas")
    axes.invert_xaxis()
    axes.set_aspect("equal")

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def plot_flux_size_posterior(samples, savefn=None):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1)
    fluxes = dict()
    sizes = dict()
    n_comps = int(len(samples[0])/4)

    for i_comp in range(n_comps):
        fluxes[i_comp] = np.exp(samples[:, 2+i_comp*4])
        sizes[i_comp] = np.exp(samples[:, 3+i_comp*4])

    for i_comp, color in zip(range(n_comps), colors):
        axes.scatter(fluxes[i_comp], sizes[i_comp], s=0.6, color=color)

    axes.set_xlabel("flux [Jy]")
    axes.set_ylabel("FWHM [mas]")
    axes.set_xscale('log')
    axes.set_yscale('log')

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def plot_flux_size_posterior_isoT(samples, freq_ghz=15.4, z=0, D=1, savefn=None):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1)
    log_fluxes = dict()
    log_sizes = dict()
    n_comps = int(len(samples[0])/4)

    for i_comp in range(n_comps):
        log_fluxes[i_comp] = samples[:, 2+i_comp*4]
        log_sizes[i_comp] = samples[:, 3+i_comp*4]

    # Find range of fluxes to plot iso-Tb lines
    lg_all_fluxes = np.log10(np.e)*np.concatenate(list(log_fluxes.values()))
    x = np.linspace(np.min(lg_all_fluxes)-0.1*np.ptp(lg_all_fluxes), np.max(lg_all_fluxes), 100)
    for lg_tb in (9, 10.5, 12, 13):
        y = lg_size_for_given_flux_and_tb(10**x, lg_tb, freq_ghz, z, D)
        axes.plot(x, y, color='black', linestyle='--', label=r"$10^{%s}$ K" % (str(lg_tb)))

    from labellines import labelLines
    lines = axes.get_lines()
    labelLines(lines, xvals=[x[int(len(x)/10)]]*len(lines))

    for i_comp, color in zip(range(n_comps), colors):
        axes.scatter(np.log10(np.e)*log_fluxes[i_comp], np.log10(np.e)*log_sizes[i_comp], s=0.6, color=color)

    axes.set_xlabel("lg(flux [Jy])")
    axes.set_ylabel("lg(size [mas])")

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def plot_tb_distance_posterior(samples, freq_ghz, z=0.0, savefn=None):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1)
    tbs = dict()
    rs = dict()
    n_comps = int(len(samples[0])/4)

    for i_comp in range(n_comps):
        rs[i_comp] = np.hypot(samples[:, 0+i_comp*4], samples[:, 1+i_comp*4])
        tbs[i_comp] = tb_comp(np.exp(samples[:, 2+i_comp*4]),
                              np.exp(samples[:, 3+i_comp*4]),
                              freq_ghz, z=z)

    for i_comp, color in zip(range(n_comps), colors):
        axes.scatter(rs[i_comp], tbs[i_comp], s=0.6, color=color)

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


def plot_size_distance_posterior(samples, savefn=None):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1)
    sizes = dict()
    rs = dict()
    n_comps = int(len(samples[0])/4)

    for i_comp in range(n_comps):
        rs[i_comp] = np.hypot(samples[:, 0+i_comp*4], samples[:, 1+i_comp*4])
        sizes[i_comp] = np.exp(samples[:, 3+i_comp*4])

    for i_comp, color in zip(range(n_comps), colors):
        axes.scatter(rs[i_comp], sizes[i_comp], s=0.6, color=color)

    axes.set_xlabel("r [mas]")
    axes.set_ylabel("FWHM [mas]")
    axes.set_xscale('log')
    axes.set_yscale('log')

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def plot_size_tb_posterior(samples, freq_ghz, z=0.0, savefn=None):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    fig, axes = plt.subplots(1, 1)
    tbs = dict()
    sizes = dict()
    n_comps = int(len(samples[0])/4)

    for i_comp in range(n_comps):
        sizes[i_comp] = np.exp(samples[:, 3+i_comp*4])
        tbs[i_comp] = tb_comp(np.exp(samples[:, 2+i_comp*4]),
                              np.exp(samples[:, 3+i_comp*4]),
                              freq_ghz, z=z)

    for i_comp, color in zip(range(n_comps), colors):
        axes.scatter(sizes[i_comp], tbs[i_comp], s=0.6, color=color)

    # Need manually set ylim because of matplotlib bug
    lg_tb_min = np.floor((np.log10(np.min([tbs[i] for i in range(n_comps)]))))
    lg_tb_max = np.ceil(np.log10(np.max([tbs[i] for i in range(n_comps)])))
    axes.set_ylim([10**lg_tb_min, 10**lg_tb_max])
    # axes.set_xlim([0.01, 1])
    axes.set_xlabel("size [mas]")
    axes.set_ylabel("Tb [K]")
    axes.set_xscale('log')
    axes.set_yscale('log')

    if savefn is not None:
        fig.savefig(savefn, dpi=300, bbox_inches="tight")
    return fig


def process_norj_samples(post_file, jitter_first=True,
                         ra_lim=(-10, 10), dec_lim=(-10, 10), freq_ghz=15.4,
                         z=0.0,
                         difmap_model_fn=None, sort_by="r",
                         savefn_position_post=None,
                         savefn_fluxsize_isot_post=None,
                         savefn_fluxsize_post=None,
                         savefn_rtb_post=None, savefn_sizer_post=None,
                         savefn_sizetb_post=None, savefn_corner=None, inverse=False,
                         comp_length=4, dfm_outname=None):
    if sort_by not in ("flux", "r", "ra", "dec"):
        raise Exception

    data = np.loadtxt(post_file)
    if jitter_first:
        data = data[:, 1:]

    if sort_by == "flux":
        data = sort_samples_by_F(data, comp_length=comp_length)
    elif sort_by == "r":
        data = sort_samples_by_r(data, comp_length=comp_length)
    elif sort_by == "ra":
        data = sort_samples_by_ra(data, comp_length=comp_length)
    elif sort_by == "dec":
        data = sort_samples_by_dec(data, inverse=inverse, comp_length=comp_length)

    fig1 = plot_position_posterior(data, savefn_position_post, ra_lim, dec_lim, difmap_model_fn, comp_length=comp_length)
    fig2 = plot_flux_size_posterior_isoT(data, freq_ghz, z, savefn=savefn_fluxsize_isot_post)
    fig3 = plot_flux_size_posterior(data, savefn=savefn_fluxsize_post)
    fig4 = plot_tb_distance_posterior(data, freq_ghz, z, savefn=savefn_rtb_post)
    fig5 = plot_size_distance_posterior(data, savefn=savefn_sizer_post)
    fig6 = plot_size_tb_posterior(data, freq_ghz, z, savefn=savefn_sizetb_post)
    fig7 = plot_corner_gen(data, savefn=savefn_corner, use_rtheta=True)
    if dfm_outname is not None:
        median_sample = np.median(data, axis=0)
        convert_sample_to_difmap_model(median_sample, dfm_outname, freq_ghz, comp_length)
    return fig1, fig2, fig3, fig4, fig5, fig6, fig7


def plot_corner_gen(samples, savefn=None, truths=None, range_frac=1.0,
                    jitter_first=False, n_jitters=1, plot_range=None,
                    comp_coordinates_to_skip=None, comp_type="circ",
                    plot_jitters=False, use_rtheta=False):

    samples_cp = samples.copy()

    if comp_type not in ("pt", "circ", "ell"):
        raise Exception("Component type must be: pt, circ or ell!")

    n_params = {"pt": 3, "circ": 4, "ell": 6}[comp_type]

    coordinates_sym_0 = r"r" if use_rtheta else "x"
    coordinates_sym_1 = r"\theta" if use_rtheta else "y"
    coordinates_sym_units_1 = "deg" if use_rtheta else "mas"

    columns = list()
    j = 0
    if jitter_first:
        j = n_jitters

    n_comps = int(len(samples_cp[0, j:]) / n_params)
    for i in range(1, n_comps+1):
        if comp_type == "pt":
            columns.append([r"${}_{}$[mas]".format(coordinates_sym_0, i), r"${}_{}$[{}]".format(coordinates_sym_1, i, coordinates_sym_units_1),
                            r"$flux_{%s}$[Jy]" % str(i)])
        elif comp_type == "circ":
            columns.append([r"${}_{}$[mas]".format(coordinates_sym_0, i), r"${}_{}$[{}]".format(coordinates_sym_1, i, coordinates_sym_units_1),
                            r"$flux_{%s}$[Jy]" % str(i), r"$bmaj_{%s}$[mas]" % str(i)])
        else:
            columns.append([r"${}_{}$[mas]".format(coordinates_sym_0, i), r"${}_{}$[{}]".format(coordinates_sym_1, i, coordinates_sym_units_1),
                            r"$flux_{%s}$[Jy]" % str(i), r"$bmaj_{%s}$[mas]" % str(i),
                            r"$e_{}$".format(i), r"$bpa_{}$[deg]".format(i)])

    for i in range(n_comps):
        if comp_type == "ell":
            # BPA
            samples_cp[:, j + 5 + i * n_params] = np.rad2deg(np.unwrap(sorted(2 * samples_cp[:, j + 5 + j * n_params])) / 2)
        # Size
        samples_cp[:, j + 3 + i * n_params] = np.exp(samples_cp[:, j + 3 + i * n_params])
        # Flux
        samples_cp[:, j + 2 + i * n_params] = np.exp(samples_cp[:, j + 2 + i * n_params])
        if use_rtheta:
            x = samples_cp[:, j + 0 + i * n_params]
            y = samples_cp[:, j + 1 + i * n_params]
            # mas
            r = np.hypot(x, y)
            # rad
            theta = np.arctan2(x, y)
            theta /= degree_to_rad
            samples_cp[:, j + 0 + i * n_params] = r
            samples_cp[:, j + 1 + i * n_params] = theta

    columns = [item for sublist in columns for item in sublist]
    if comp_coordinates_to_skip is not None:
        if not jitter_first:
            idx_to_delete = [0 + comp_coordinates_to_skip*n_params,
                             1 + comp_coordinates_to_skip*n_params]
        else:
            idx_to_delete = [n_jitters + 0 + comp_coordinates_to_skip*n_params,
                             n_jitters + 1 + comp_coordinates_to_skip*n_params]
        columns = columns[:0 + comp_coordinates_to_skip*n_params] + \
                  columns[2 + comp_coordinates_to_skip*n_params:]
        samples_cp = np.delete(samples_cp, idx_to_delete, axis=1)
    if jitter_first and plot_jitters:
        if n_jitters == 1:
            columns.insert(0, r"$\log{\sigma_{\rm jitter}}$")
        else:
            for i in range(1, n_jitters+1):
                columns.insert(0, r"$\log{\sigma_{%s}_{\rm jitter}}$" % str(i))
    if plot_range is None:
        plot_range = [range_frac] * len(columns)
    fig = corner.corner(samples_cp, labels=columns, truths=truths,
                        show_titles=True, quantiles=[0.16, 0.5, 0.84],
                        color="gray", truth_color="#1f77b4",
                        plot_contours=True, range=plot_range, title_fmt=".3f", title_kwargs={"size": 10},
                        plot_datapoints=False, fill_contours=True,
                        levels=(0.393, 0.865, 0.989),
                        hist2d_kwargs={"plot_datapoints": False,
                                       "plot_density": False,
                                       "plot_contours": True,
                                       "no_fill_contours": True},
                        hist_kwargs={'ls': 'solid',
                                     'density': True})
    if savefn is not None:
        fig.savefig(savefn, dpi=100, bbox_inches="tight")
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
                type_ = "1"
                bmaj = "{:.7f}v".format(bmaj)
            elif len(comp) == 6:
                # Jy, mas, mas, mas, -, deg
                flux, x, y, bmaj, e, bpa = comp
                e = "{}v".format(e)
                bpa = "{}v".format((bpa-np.pi/2)/degree_to_rad)
                bmaj = "{}v".format(bmaj)
                type_ = "1"
            elif len(comp) == 3:
                flux, x, y = comp
                e = "1.00000"
                bmaj = "0.0000"
                bpa = "000.000"
                type_ = "0"
            else:
                raise Exception
            # mas
            r = np.hypot(x, y)
            # rad
            theta = np.arctan2(x, y)
            theta /= degree_to_rad
            fo.write("{:>11.7f}v {:>13.7f}v {:>13.5f}v {:>13} {:>13} {:>13} {:>3} {:>12.5e} {:>12d}\n".format(flux, r, theta,
                                                              bmaj, e, bpa, type_,
                                                             freq_hz, 0))


def convert_sample_to_difmap_model(sample, out_fname, freq_ghz, comp_length=4):
    if comp_length != 4:
        raise NotImplementedError("Only Circular Gaussians are implemented")
    n_comps = int(len(sample)/comp_length)
    components = list()
    for i in range(n_comps):
        subsample = sample[comp_length*i:comp_length*i+comp_length]
        print(subsample)
        flux = np.exp(subsample[2])
        bmaj = np.exp(subsample[3])
        x = subsample[0]
        y = subsample[1]
        cg = (flux, x, y, bmaj)
        components.append(cg)
    components = sorted(components, key=lambda comp: comp[0], reverse=True)
    export_difmap_model(components, out_fname, 1e9*freq_ghz)