import os
import numpy as np
from astropy import units as u
import astropy.io.fits as pf
from astropy.stats import mad_std
import matplotlib.pyplot as plt
import pandas as pd
import ehtim as eh
import sys
from subprocess import Popen, PIPE

mas_to_rad = u.mas.to(u.rad)


def time_average(uvfits, outfname, time_sec=60, show_difmap_output=True,
                 reweight=True):
    if reweight:
        cmd = "observe " + uvfits + ", {}, true\n".format(time_sec)
    else:
        cmd = "observe " + uvfits + ", {}, false\n".format(time_sec)
    cmd += "wobs {}\n".format(outfname)
    cmd += "exit\n"

    with Popen('difmap', stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True) as difmap:
        outs, errs = difmap.communicate(input=cmd)
    if show_difmap_output:
        print(outs)
        print(errs)


def point_ft(uv, ra, dec, flux):
    """
    FT of delta function at ``uv`` points.

    :param flux:
        Full flux [Jy].
    :param ra:
        Distance from phase center [mas].
    :param dec:
        Distance from phase center [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [ra*mas_to_rad, dec*mas_to_rad]
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    ft = flux
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def gaussian_circ_ft(uv, ra, dec, flux, bmaj):
    """
    FT of circular gaussian at ``uv`` points.

    :param flux:
        Full flux of gaussian [Jy].
    :param ra:
        Distance from phase center [mas].
    :param dec:
        Distance from phase center [mas].
    :param bmaj:
        FWHM of a gaussian [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [ra*mas_to_rad, dec*mas_to_rad]
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = uv[:, 0]**2 + uv[:, 1]**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def gaussian_ell_ft(uv, ra, dec, flux, bmaj, e, bpa):
    """
    FT of elliptical gaussian at ``uv`` points.

    :param flux:
        Full flux of gaussian [Jy].
    :param ra:
        Distance from phase center [mas].
    :param dec:
        Distance from phase center [mas].
    :param bmaj:
        FWHM of a gaussian [mas].
    :param e:
        Eccentricity.
    :param bpa:
        Positional angle [rad]. As usual - from N to E (positive RA).
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [ra*mas_to_rad, dec*mas_to_rad]
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = e**2 * (uv[:, 0]*np.cos(bpa) - uv[:, 1]*np.sin(bpa))**2 + (uv[:, 0]*np.sin(bpa) + uv[:, 1]*np.cos(bpa))**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def optically_thin_sphere_ft(uv, ra, dec, flux, size):
    """
    FT of circular gaussian at ``uv`` points.

    :param flux:
        Full flux [Jy].
    :param ra:
        Distance from phase center [mas].
    :param dec:
        Distance from phase center [mas].
    :param size:
        Size of a sphere [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [ra*mas_to_rad, dec*mas_to_rad]
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    pi_D_rho = np.pi*size*mas_to_rad*np.sqrt(uv[:, 0]**2 + uv[:, 1]**2)
    ft = 3*flux*(np.sin(pi_D_rho) - pi_D_rho*np.cos(pi_D_rho))/pi_D_rho**3
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def get_rms(obs, t1, t2):
    rec = obs.unpack_bl(t1, t2, "vis")
    real = rec["vis"].real
    imag = rec["vis"].imag
    sigma_real = mad_std(np.diff(real, axis=0))/np.sqrt(2)
    sigma_imag = mad_std(np.diff(imag, axis=0))/np.sqrt(2)
    return 0.5*(sigma_real + sigma_imag)


# For now better average using difmap, that assures that weights are reasonable
# 08.05.2024: Now it uses successive differences approach and does not trust the weights
def get_data_file_from_ehtim(uvfits, outname, avg_time_sec=0, average_using="difmap"):

    assert average_using in ("difmap", "eht-imager")

    if avg_time_sec > 0:
        if average_using == "difmap":
            uvfits_dir, uvfits_fname = os.path.split(uvfits)
            uvfits_ta = os.path.join(uvfits_dir, f"ta{avg_time_sec}_{uvfits_fname}")
            # difmap uses vector weighted averaging
            time_average(uvfits, uvfits_ta, time_sec=avg_time_sec)
            uvfits = uvfits_ta
    # eht-imager estimates errors from weights!
    obs = eh.obsdata.load_uvfits(uvfits)

    if avg_time_sec > 0:
        if average_using == "eht":
            obs = obs.avg_coherent(avg_time_sec)

    rec = obs.unpack(["t1", "t2", "u", "v", "vis", "sigma"])
    df = pd.DataFrame.from_records(rec)
    df["vis_re"] = np.real(df["vis"])
    df["vis_im"] = np.imag(df["vis"])
    # From weights
    df["error"] = df["sigma"]
    # Successive differences approach
    # df["error"] = df.apply(lambda x: get_rms(obs, x.t1, x.t2), axis=1)
    df = df[["t1", "t2", "u", "v", "vis_re", "vis_im", "error"]]
    df = df.replace({"t1": obs.tkey})
    df = df.replace({"t2": obs.tkey})
    df.to_csv(outname, sep=",", header=True, index=False)
    return df


def get_data_file_from_ehtim_multiple_IFs(uvfits, outname, avg_time_sec=0, average_using="difmap"):

    # Find n_IF
    header = pf.getheader(uvfits)
    n_IF = None
    if header["CTYPE5"] == "IF":
        n_IF = int(header["NAXIS5"])
    else:
        raise Exception("Can't find n_IF from header")

    # Find frequencies of all IFs
    # Origin frequency of first IF
    ref_fr = header["CRVAL4"]
    # Width of the IF
    IF_bandwidth = header["CDELT4"]
    IF_center_freqs = np.ones(n_IF)*ref_fr + (0.5 + np.arange(n_IF))*IF_bandwidth


    assert average_using in ("difmap", "eht-imager")

    if avg_time_sec > 0:
        if average_using == "difmap":
            uvfits_dir, uvfits_fname = os.path.split(uvfits)
            uvfits_ta = os.path.join(uvfits_dir, f"ta{avg_time_sec}_{uvfits_fname}")
            # difmap uses vector weighted averaging
            time_average(uvfits, uvfits_ta, time_sec=avg_time_sec)
            uvfits = uvfits_ta

    dfs = list()
    for i_IF in range(n_IF):
        # eht-imager estimates errors from weights!
        obs = eh.obsdata.load_uvfits(uvfits, IF=[i_IF])
        obs.tlist()
        if avg_time_sec > 0:
            if average_using == "eht":
                obs = obs.avg_coherent(avg_time_sec)



        rec = obs.unpack(["t1", "t2", "u", "v", "vis", "sigma"])
        # Correct for wrong conversion of (u,v) from light seconds to lambdas
        rec["u"] = IF_center_freqs[i_IF]/ref_fr * rec["u"]
        rec["v"] = IF_center_freqs[i_IF]/ref_fr * rec["v"]
        df = pd.DataFrame.from_records(rec)

        # Queries
        # query("u == -@row.u and v == -@row.v")

        df["vis_re"] = np.real(df["vis"])
        df["vis_im"] = np.imag(df["vis"])
        df["error"] = df["sigma"]

        # df["error"] = df.apply(lambda x: get_rms(obs, x.t1, x.t2), axis=1)

        df = df[["t1", "t2", "u", "v", "vis_re", "vis_im", "error"]]
        df = df.replace({"t1": obs.tkey})
        df = df.replace({"t2": obs.tkey})
        dfs.append(df)

    df = pd.concat(dfs, axis=0, ignore_index=False)
    df.to_csv(outname, sep=",", header=True, index=False)
    return df


def add_noise(df, use_global_median_noise=True, global_noise_scale=None, logjitters_dict=None, offsets_dict=None):
    """
    Add noise as specified in ``error`` columns to ``vis_re`` and ``vis_im`` columns.

    :param df:
        DataFrame with columns ``error``, ``vis_re`` and ``vis_im``.
    :param use_global_median_noise: (optional)
        Noise is estimated using median over all visibilities. Useful if template uv-data
        has shitty baselines with high noise. (default: ``True``)
    :param logjitters_dict: (optional)
        Dictionary with keys - antenna numbers (as in ``df``) and values - logjitters for each
        antenna. If ``None`` then do not add per-antenna jitters. (default: ``None``)
    :param offsets_dict: (optional)
        Dictionary with keys - antenna numbers (as in ``df``) and values - systematical amplitude offsets
        for each antenna. If ``None`` then do not add offsets. (default: ``None``)
    :return:
        New DataFrame with noise added.
    """
    df_ = df.copy()

    if offsets_dict is not None:
        df_["offset_i"] = df_["t1"].map(lambda x: offsets_dict[x])
        df_["offset_j"] = df_["t2"].map(lambda x: offsets_dict[x])
        df_["vis_re"] *= (df_["offset_i"]*df_["offset_j"])
        df_["vis_im"] *= (df_["offset_i"]*df_["offset_j"])
        df_.drop(["offset_i", "offset_j"], axis=1, inplace=True)

    if logjitters_dict is not None:
        df_["logjitter_i"] = df_["t1"].map(lambda x: logjitters_dict[x])
        df_["logjitter_j"] = df_["t2"].map(lambda x: logjitters_dict[x])
        df_["vis_re"] += np.hypot(df_["logjitter_i"].map(lambda x: np.random.normal(0, np.exp(x), size=1)[0]),
                                  df_["logjitter_j"].map(lambda x: np.random.normal(0, np.exp(x), size=1)[0]))
        df_["vis_im"] += np.hypot(df_["logjitter_i"].map(lambda x: np.random.normal(0, np.exp(x), size=1)[0]),
                                  df_["logjitter_j"].map(lambda x: np.random.normal(0, np.exp(x), size=1)[0]))
        df_.drop(["logjitter_i", "logjitter_j"], axis=1, inplace=True)

    if use_global_median_noise:
        error = df_["error"].median()
        if global_noise_scale is not None:
            error = global_noise_scale*error
            df_["error"] *= global_noise_scale
        df_["vis_re"] += np.random.normal(0, error, size=df.shape[0])
        df_["vis_im"] += np.random.normal(0, error, size=df.shape[0])
    # Use individual visibilities noise estimated
    else:
        df_["vis_re"] += df_["error"].map(lambda x: np.random.normal(0, x, size=1)[0])
        df_["vis_im"] += df_["error"].map(lambda x: np.random.normal(0, x, size=1)[0])
    return df_


def radplot(df, fig=None, color=None, label=None, style="ap", savefname=None, show=True):
    uv = df[["u", "v"]].values
    r = np.hypot(uv[:, 0], uv[:, 1])/10**6

    if style == "ap":
        value1 = np.hypot(df["vis_re"].values, df["vis_im"].values)
        value2 = np.arctan2(df["vis_im"].values, df["vis_re"].values)
    elif style == "reim":
        value1 = df["vis_re"].values
        value2 = df["vis_im"].values
        error = df["error"].values
    else:
        raise Exception("style can be ap or reim")

    if fig is None:
        fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    else:
        axes = fig.axes

    if color is None:
        color = "#1f77b4"

    if style == "ap":
        axes[0].plot(r, value1, '.', color=color, label=label)
        axes[1].plot(r, value2, '.', color=color)
    else:
        axes[0].errorbar(r, value1, yerr=error, fmt=".", color=color, label=label, capsize=3)
        axes[1].errorbar(r, value2, yerr=error, fmt=".", color=color)
    if style == "ap":
        axes[0].set_ylabel("Amp, Jy")
        axes[1].set_ylabel("Phase, rad")
    else:
        axes[0].set_ylabel("Re, Jy")
        axes[1].set_ylabel("Im, Jy")
    axes[1].set_xlabel(r"$r_{\rm uv}$, $M\lambda$")
    if label is not None:
        axes[0].legend()
    plt.tight_layout()
    if savefname:
        fig.savefig(savefname, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    return fig


if __name__ == "__main__":
    import glob
    uvfits_dir = "/home/ilya/data/VLBI_Gaia/2comp"
    uvfits_files = glob.glob(os.path.join(uvfits_dir, "*_vis.fits"))
    for uvfits_file in uvfits_files:
        print("Processing UVFITS : ", uvfits_file)
        uvfits_fn = os.path.split(uvfits_file)[-1]
        source, band, year, month, day, author, product = uvfits_fn.split("_")
        product = product.split(".")[0]
        assert product == "vis"
        save_basename = f"{source}_{band}_{year}_{month}_{day}"
        out_fn = f"{save_basename}_{author}_{product}.csv"
        out_fname = os.path.join(uvfits_dir, out_fn)
        df = get_data_file_from_ehtim(uvfits_file, out_fname, avg_time_sec=60, average_using="difmap")

    sys.exit(0)


    # Artificial source
    # Zero observed data
    df["vis_re"] = 0
    df["vis_im"] = 0
    # Add model
    re, im = gaussian_ell_ft(flux=1.0, ra=0.0, dec=0.0, bmaj=0.2, e=0.2, bpa=np.pi/3., uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_ell_ft(flux=0.5, ra=1.0, dec=0.1, bmaj=0.3, e=0.5, bpa=np.pi/4., uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_ell_ft(flux=0.25, ra=2.5, dec=0.5, bmaj=0.5, e=0.7, bpa=np.pi/2., uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_ell_ft(flux=0.1, ra=10., dec=-1.0, bmaj=2.0, e=0.8, bpa=np.pi/3., uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im

    try:
        fig = radplot(df, label="Data", show=False)
    except AttributeError:
        fig = radplot(df, label="Data", style="reim")

    df_updated = add_noise(df, use_global_median_noise=True, global_noise_scale=1)
                           # logjitters_dict={0: -8, 1: -8, 2: -8, 3: -8, 4: -8, 5: -8, 6: -8, 7: -8, 8: -3, 9: -4},
                           # offsets_dict={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1.3, 9: 1})
    fig = radplot(df_updated, color="#ff7f0e", fig=fig, label="With noise")
    df_updated.to_csv(out_fname, sep=",", index=False, header=True)
