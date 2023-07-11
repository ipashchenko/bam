import os
import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
from astropy import units
import pandas as pd
from tqdm import tqdm
import sys
sys.path.insert(0, '/home/ilya/github/agn_abc')
from data import Data
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from spydiff import time_average

mas_to_rad = units.mas.to(units.rad)


def gaussian_circ_ft(flux, RA, DEC, bmaj, uv):
    """
    FT of circular gaussian at ``uv`` points.

    :param flux:
        Full flux of gaussian.
    :param RA:
        Distance from phase center [mas].
    :param DEC:
        Distance from phase center [mas].
    :param bmaj:
        FWHM of a gaussian [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [RA*mas_to_rad, DEC*mas_to_rad]
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = uv[:, 0]**2 + uv[:, 1]**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def create_data_file(uvfits, outfile=None):
    """
    :param uvfits:
        Path to UVFITS file.
    :param outfile:
        Path to output txt-file.
    """
    hdus = pf.open(uvfits)
    data = hdus[0].data
    header = hdus[0].header
    freq = header["CRVAL4"]

    df = pd.DataFrame(columns=["u", "v", "vis_re", "vis_im", "error"])

    for group in data:
        try:
            u = group["UU"]
            v = group["VV"]
        except KeyError:
            u = group["UU--"]
            v = group["VV--"]
        if abs(u) < 1.:
            u *= freq
            v *= freq
        # IF, STOKES, COMPLEX
        # data_ = group["DATA"].squeeze()
        data_ = group["DATA"][0, 0, :, 0, :, :]
        weights = data_[:, :, 2]
        mask = weights <= 0
        mask = np.repeat(mask[:, :, np.newaxis], 3, axis=2)
        masked_data = np.ma.array(data_, mask=mask)
        if np.alltrue(weights <= 0):
            continue
        weights = np.ma.array(weights, mask=mask[..., 0])
        weight = np.ma.sum(weights)
        error = 1/np.sqrt(weight)

        # FIXME: Comment out this block for some specific STOKES data sets (e.g. 3C84 RA)
        if data_.shape[1] > 1:
            # print("Stokes I")
            vis_re = 0.5*(masked_data[:, 0, 0] + masked_data[:, 1, 0])
            vis_im = 0.5*(masked_data[:, 0, 1] + masked_data[:, 1, 1])
        else:
            # print("Stokes RR or LL")
            vis_re = masked_data[:, 0, 0]
            vis_im = masked_data[:, 0, 1]

        # # For LL (3 84 RA)
        # vis_re = masked_data[:, 1, 0]
        # vis_im = masked_data[:, 1, 1]

        vis_re = np.ma.mean(vis_re)
        vis_im = np.ma.mean(vis_im)
        df_ = pd.Series({"u": u, "v": v, "vis_re": vis_re, "vis_im": vis_im, "error": error})
        df = df.append(df_, ignore_index=True)

    if outfile is not None:
        df.to_csv(outfile, sep=" ", index=False, header=True)
    return df


def create_data_file_v2(uvfits, outfile, time_average_sec=None, error_from_weight=False):
    """
    :param uvfits:
        Path to UVFITS file.
    :param outfile:
        Path to output txt-file.
    """
    if time_average_sec is not None:
        path, fname = os.path.split(uvfits)
        to_read_uvfits = os.path.join(path, "ta{}sec_{}".format(time_average_sec, fname))
        time_average(uvfits, to_read_uvfits, time_sec=time_average_sec)
    else:
        to_read_uvfits = uvfits
    all_data = Data(to_read_uvfits)
    if error_from_weight:
        error = all_data.error
    else:
        error = all_data.error_wt
    df = pd.DataFrame.from_dict({"u": all_data.uv[:, 0],
                                 "v": all_data.uv[:, 1],
                                 "vis_re": all_data.data.real,
                                 "vis_im": all_data.data.imag,
                                 "error": error})
    df = df.dropna()
    if outfile is not None:
        df.to_csv(outfile, sep=" ", index=False, header=False)
    return df


def add_noise(df, use_global_median_noise=True, global_noise_scale=None):
    """
    Add noise as specified in ``error`` columns to ``vis_re`` and ``vis_im`` columns.

    :param df:
        DataFrame with columns ``error``, ``vis_re`` and ``vis_im``.
    :param use_global_median_noise: (optional)
        Noise is estimated using median over all visibilities. Useful if template uv-data
        has shitty baselines with high noise. (default: ``True``)
    :return:
        New DataFrame with noise added.
    """
    df_ = df.copy()

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


def radplot(df, fig=None, color=None, label=None, style="ap"):
    u = df["u"].values
    v = df["v"].values
    uv = np.vstack((u, v)).T
    r = np.hypot(uv[:, 0], uv[:, 1])

    if style == "ap":
        value1 = np.hypot(df["vis_re"].values, df["vis_im"].values)
        value2 = np.arctan2(df["vis_im"].values, df["vis_re"].values)
    elif style == "reim":
        value1 = df["vis_re"].values
        value2 = df["vis_im"].values
    else:
        raise Exception("style can be ap or reim")

    import matplotlib.pyplot as plt
    if fig is None:
        fig, axes = plt.subplots(2, 1, sharex=True)
    else:
        axes = fig.axes

    if color is None:
        color = "#1f77b4"

    axes[0].plot(r, value1, '.', color=color, label=label)
    axes[1].plot(r, value2, '.', color=color)
    if style == "ap":
        axes[0].set_ylabel("Amp, Jy")
        axes[1].set_ylabel("Phase, rad")
    else:
        axes[0].set_ylabel("Re, Jy")
        axes[1].set_ylabel("Im, Jy")
    axes[1].set_xlabel(r"$r_{\rm uv}$, $\lambda$")
    if label is not None:
        axes[0].legend()
    plt.tight_layout()
    fig.show()
    return fig


if __name__ == "__main__":

    # TODO: Use this to create data files for several MOJAVE epochs
    # import glob
    # source = "2233-148"
    # uvfits_files = glob.glob(os.path.join("/home/ilya/github/bam/data/{}".format(source), "{}.u.*.uvf".format(source)))
    # for uvfits_file in uvfits_files:
    #     fname = os.path.split(uvfits_file)[-1]
    #     epoch = fname.split(".")[2]
    #     out_fname = "/home/ilya/github/bam/data/{}/{}_{}_120s.txt".format(source, source, epoch)
    #     df = create_data_file_v2(uvfits_file, out_fname, time_average_sec=120)

    out_fname = "/home/ilya/github/bam/data/1800+7828/1800+7828_1998_10_01_X_120s.txt"
    uvfits_file = "/home/ilya/github/bam/data/J1800+7828_X_1998_10_01_pus_vis.fits"
    df = create_data_file_v2(uvfits_file, out_fname, time_average_sec=120)
    sys.exit(0)


# # Zero observed data
    # df["vis_re"] = 0
    # df["vis_im"] = 0
    # # Add model
    # re, im = gaussian_circ_ft(flux=2.0, dx=0.0, dy=0.0, bmaj=0.1, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=1.0, dx=0.5, dy=0.0, bmaj=0.2, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.5, dx=1.5, dy=1.0, bmaj=0.3, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.25, dx=3.5, dy=2.0, bmaj=0.5, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.125, dx=5, dy=5.0, bmaj=0.5, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.075, dx=7.5, dy=8.0, bmaj=0.5, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.05, dx=8.0, dy=9.0, bmaj=0.75, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im
    # re, im = gaussian_circ_ft(flux=0.035, dx=9.0, dy=9.5, bmaj=0.75, uv=df[["u", "v"]].values)
    # df["vis_re"] += re
    # df["vis_im"] += im

    try:
        fig = radplot(df, label="Data")
    except AttributeError:
        fig = radplot(df, label="Data", style="reim")

    #
    # # # Add noise and plot
    # df_updated = add_noise(df, use_global_median_noise=True, global_noise_scale=10)
    # fig = radplot(df_updated, color="#ff7f0e", fig=fig, label="With noise")
    df_updated = df
    df_updated.to_csv(out_fname, sep=" ", index=False, header=False)