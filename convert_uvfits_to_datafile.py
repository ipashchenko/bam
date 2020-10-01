import os
import numpy as np
import pandas as pd
import astropy.io.fits as pf
from astropy import units as u
import datetime
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
    # Note "+" sign here. 
    result = np.exp(2.0*np.pi*1j*(uv @ shift))
    c = (np.pi*bmaj*mas_to_rad)**2/(4. * np.log(2.))
    b = uv[:, 0]**2 + uv[:, 1]**2
    ft = flux*np.exp(-c*b)
    ft = np.array(ft, dtype=complex)
    result *= ft
    return result.real, result.imag


def time_average(uvfits, outfname, time_sec=120, show_difmap_output=True,
                 reweight=True):
    stamp = datetime.datetime.now()
    command_file = "difmap_commands_{}".format(stamp.isoformat())

    difmapout = open(command_file, "w")
    if reweight:
        difmapout.write("observe " + uvfits + ", {}, true\n".format(time_sec))
    else:
        difmapout.write("observe " + uvfits + ", {}, false\n".format(time_sec))
    difmapout.write("wobs {}\n".format(outfname))
    difmapout.write("exit\n")
    difmapout.close()
    # TODO: Use subprocess for silent cleaning?
    shell_command = "difmap < " + command_file + " 2>&1"
    if not show_difmap_output:
        shell_command += " >/dev/null"
    os.system(shell_command)

    # Remove command file
    os.unlink(command_file)


def create_data_file(uvfits, outfile=None, keep_header=False, time_average_sec=None):
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

    hdus = pf.open(to_read_uvfits)
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

        # FIXME: Comment out this block for some specific STOKES data sets (e.g. 3C84 RA with only Stokes I present)
        if data_.shape[1] > 1:
            vis_re = 0.5*(masked_data[:, 0, 0] + masked_data[:, 1, 0])
            vis_im = 0.5*(masked_data[:, 0, 1] + masked_data[:, 1, 1])
        else:
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
        if keep_header:
            header = True
        else:
            header = False
        df.to_csv(outfile, sep=" ", index=False, header=header)
    return df


def plot_model_predictions(post_samples, data_df, jitter_first=True, style="reim", savefname=None, show=True,
                           n_samples_to_plot=24):
    if style not in ("reim", "ap"):
        raise Exception("Only reim or ap style plotting is supported!")

    uv = data_df[["u", "v"]].values
    r = np.hypot(uv[:, 0], uv[:, 1])/10**6
    # dr = np.max(r) - np.min(r)
    # rr = np.linspace(0, np.max(r)+0.1*dr, 1000)
    df = data_df.copy()

    # Zero visibilities
    df["vis_re"] = 0
    df["vis_im"] = 0

    # Add model
    n_post, n_components = post_samples.shape
    if jitter_first:
        n_components -= 1
    n_components = int(n_components/4)
    samples = post_samples[np.random.randint(0, n_post, n_samples_to_plot)]
    samples_predictions_re = list()
    samples_predictions_im = list()

    if jitter_first:
        jitter = 1
    else:
        jitter = 0

    for sample in samples:
        sample_prediction_re = np.zeros(len(df))
        sample_prediction_im = np.zeros(len(df))
        for n_comp in range(n_components):
            dx, dy, flux, bmaj = sample[jitter+n_comp*4:jitter+(n_comp+1)*4]
            flux = np.exp(flux)
            bmaj = np.exp(bmaj)
            re, im = gaussian_circ_ft(flux=flux, dx=dx, dy=dy, bmaj=bmaj, uv=uv)
            sample_prediction_re += re
            sample_prediction_im += im
        samples_predictions_re.append(sample_prediction_re)
        samples_predictions_im.append(sample_prediction_im)

    # Original data plot
    fig = radplot(data_df, style=style, show=False)
    axes = fig.get_axes()
    for sample_prediction_re, sample_prediction_im in zip(samples_predictions_re, samples_predictions_im):
        if style == "reim":
            axes[0].plot(r, sample_prediction_re, "_", alpha=0.25, color="C1", ms=5, mew=0.2)
            axes[1].plot(r, sample_prediction_im, "_", alpha=0.25, color="C1", ms=5, mew=0.2)
        else:
            axes[0].plot(r, np.hypot(sample_prediction_re, sample_prediction_im), "_", alpha=0.25, color="C1", ms=5, mew=0.2)
            axes[1].plot(r, np.arctan2(sample_prediction_im, sample_prediction_re), "_", alpha=0.25, color="C1", ms=5, mew=0.2)
    if savefname is not None:
        fig.savefig(savefname, bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    return fig


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
    uvfits = "/home/ilya/github/bam/data/0716+714/0716+714.u.2006_12_01.uvf"
    df = create_data_file(uvfits)
    fig = radplot(df, style="reim")
