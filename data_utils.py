import os
import numpy as np
import astropy.io.fits as pf
from astropy.time import Time
from astropy import units
import pandas as pd
from tqdm import tqdm
import sys
import ehtim as eh
from subprocess import Popen, PIPE

mas_to_rad = units.mas.to(units.rad)


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


# For now better average using difmap, that assures that weights are reasonable
# 08.05.2024: Now it uses successive differences approach and does not trust the weights
def get_data_file_from_ehtim(uvfits, outname, avg_time_sec=0, average_using="difmap",
                             working_dir=None):

    assert average_using in ("difmap", "eht-imager")
    if working_dir is None:
        working_dir = os.getcwd()

    if avg_time_sec > 0:
        if average_using == "difmap":
            uvfits_dir, uvfits_fname = os.path.split(uvfits)
            uvfits_ta = os.path.join(working_dir, f"ta{avg_time_sec}_{uvfits_fname}")
            # difmap uses vector weighted averaging
            time_average(uvfits, uvfits_ta, time_sec=avg_time_sec)
            uvfits = uvfits_ta
    # eht-imager estimates errors from weights!
    obs = eh.obsdata.load_uvfits(uvfits)

    os.unlink(uvfits_ta)

    if avg_time_sec > 0:
        if average_using == "eht-imager":
            obs = obs.avg_coherent(avg_time_sec)

    rec = obs.unpack(["u", "v", "vis", "sigma"])
    df = pd.DataFrame.from_records(rec)
    df["vis_re"] = np.real(df["vis"])
    df["vis_im"] = np.imag(df["vis"])
    # From weights
    df["error"] = df["sigma"]
    # Successive differences approach
    # df["error"] = df.apply(lambda x: get_rms(obs, x.t1, x.t2, vis_type="vis"), axis=1)
    df = df[["u", "v", "vis_re", "vis_im", "error"]]
    df.to_csv(outname, sep=",", header=True, index=False)
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

