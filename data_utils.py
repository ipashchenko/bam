import os
import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
import pandas as pd
import ehtim as eh
from tqdm import tqdm
import sys
sys.path.insert(0, '/home/ilya/github/agn_abc')
from data import Data
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from spydiff import time_average

mas_to_rad = u.mas.to(u.rad)


def gaussian_circ_ft(flux, dx, dy, bmaj, uv):
    """
    FT of circular gaussian at ``uv`` points.

    :param flux:
        Full flux of gaussian.
    :param dx:
        Distance from phase center [mas].
    :param dy:
        Distance from phase center [mas].
    :param bmaj:
        FWHM of a gaussian [mas].
    :param uv:
        2D numpy array of (u,v)-coordinates (dimensionless).
    :return:
        Tuple of real and imaginary visibilities parts.
    """
    shift = [dx*mas_to_rad, dy*mas_to_rad]
    result = np.exp(-2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = uv[:, 0]**2 + uv[:, 1]**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


# For now better average using difmap, that assures that weights are reasonable
def get_data_file_from_ehtim(uvfits, outname, avg_time_sec=0, average_using="difmap"):
    assert average_using in ("difmap", "eht")
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
    df["error"] = df["sigma"]
    df = df[["t1", "t2", "u", "v", "vis_re", "vis_im", "error"]]
    df = df.replace({"t1": obs.tkey})
    df = df.replace({"t2": obs.tkey})
    df.to_csv(outname, sep=" ", header=False, index=False)
    return df


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
        axes[0].errorbar(r, value1, yerr=error, fmt=".", color=color, label=label)
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
    # uvfits_file = "/home/ilya/Downloads/0212+735/u/0212+735.u.2019_08_15.uvf"
    # out_fname = "/home/ilya/Downloads/0212+735/u/0212+735.u.2019_08_15.txt"
    uvfits_file = "/home/ilya/Downloads/mojave/0851+202/0851+202.u.2023_07_01.uvf"
    out_fname = "/home/ilya/Downloads/mojave/0851+202/3comp.txt"
    # df = create_data_file_v2(uvfits_file, out_fname, time_average_sec=60)

    # This uses antenna sites for per-antenna jitter.
    df = get_data_file_from_ehtim(uvfits_file, out_fname, avg_time_sec=60, average_using="difmap")


    sys.exit(0)

    # Artificial source
    # Zero observed data
    df["vis_re"] = 0
    df["vis_im"] = 0
    # Add model
    re, im = gaussian_circ_ft(flux=2.0, dx=0.0, dy=0.0, bmaj=0.1, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_circ_ft(flux=1.0, dx=1.0, dy=0.0, bmaj=0.5, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im
    re, im = gaussian_circ_ft(flux=0.5, dx=2.5, dy=2.0, bmaj=0.5, uv=df[["u", "v"]].values)
    df["vis_re"] += re
    df["vis_im"] += im

    df_updated = add_noise(df, use_global_median_noise=True, global_noise_scale=1)
                           # logjitters_dict={0: -8, 1: -8, 2: -8, 3: -8, 4: -8, 5: -8, 6: -8, 7: -8, 8: -3, 9: -4},
                           # offsets_dict={0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1.3, 9: 1})
    df_updated.to_csv(out_fname, sep=" ", index=False, header=False)

    sys.exit(0)



    import glob
    uvfits_files = glob.glob(os.path.join("/home/ilya/github/dterms/data/mojave/", "0506+056.u*.uvf"))
    for uvfits_file in uvfits_files:
        fname = os.path.split(uvfits_file)[-1]
        epoch = fname.split(".")[2]
        out_fname = "/home/ilya/github/bam/data/0506+056_{}_120s.txt".format(epoch)
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
