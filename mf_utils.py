import os
import sys
import astropy.io.fits as pf
import matplotlib
import numpy as np
import pandas as pd
# matplotlib.use("TkAgg")
from data_utils import time_average, get_data_file_from_ehtim, add_noise, radplot, gaussian_circ_ft, mas_to_rad
sys.path.insert(0, 've/vlbi_errors')
from uv_data import UVData
from subprocess import Popen, PIPE
from scipy.optimize import newton
# matplotlib.use("TkAgg")


def load_from_fits(uvfits):
    hdus = pf.open(uvfits, mode="readonly")
    header = hdus[0].header
    freq = header["CRVAL{}".format([card for card in header.cards if "FREQ" in card][0][0][-1])]
    try:
        u = freq*hdus[0].data["UU"]
        v = freq*hdus[0].data["VV"]
    except KeyError:
        u = freq*hdus[0].data["UU--"]
        v = freq*hdus[0].data["VV--"]
    # STOKES = RR, COMPLEX = RE
    rr_re = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 0, 0], axis=1)
    rr_im = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 0, 1], axis=1)
    ll_re = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 1, 0], axis=1)
    ll_im = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 1, 1], axis=1)
    rr_we = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 0, 2], axis=1)
    ll_we = np.nanmean(hdus[0].data["DATA"][:, 0, 0, 0, :, 1, 2], axis=1)
    error_rr = np.nanmedian(1/np.sqrt(rr_we))
    error_ll = np.nanmedian(1/np.sqrt(ll_we))
    error = 0.5*(error_ll + error_rr)
    # error = 0.005
    print("Error = ", error)
    df = pd.DataFrame.from_dict({"u": u,
                                 "v": v,
                                 "vis_re": 0.5*(rr_re + ll_re),
                                 "vis_im": 0.5*(rr_im + ll_im),
                                 "error": error})
    return df


def save_to_uvfits(uvfits, re=None, im=None):
    hdus = pf.open(uvfits, mode="update")
    u = hdus[0].data["UU"]
    v = hdus[0].data["VV"]
    header = hdus[0].header
    freq = header["CRVAL{}".format([card for card in header.cards if "FREQ" in card][0][0][-1])]
    uv = np.vstack((u, v)).T
    if re is None or im is None:
        re, im = gaussian_circ_ft(flux=2.0, RA=2.0, DEC=5.0, bmaj=0.3, uv=uv*freq)
    # STOKES = RR, COMPLEX = RE
    hdus[0].data["DATA"][:, :, :, :, :, 0, 0] = np.array(re)[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    # STOKES = RR, COMPLEX = IM
    hdus[0].data["DATA"][:, :, :, :, :, 0, 1] = np.array(im)[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    # STOKES = LL, COMPLEX = RE
    hdus[0].data["DATA"][:, :, :, :, :, 1, 0] = np.array(re)[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    # STOKES = LL, COMPLEX = IM
    hdus[0].data["DATA"][:, :, :, :, :, 1, 1] = np.array(im)[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis]
    hdus.flush(verbose=True)


def create_data_files_from_real_uvfits(band_uvfits_files_dict, save_dir, time_average_sec_dict):
    import matplotlib.pyplot as plt
    matplotlib.use("TkAgg")
    for band, uvfits in band_uvfits_files_dict.items():
        uvdata = UVData(uvfits)
        freq_ghz = uvdata.frequency/1E+09
        out_fname = os.path.join(save_dir, "{}_{:.2f}".format(band, freq_ghz))
        df = get_data_file_from_ehtim(uvfits, out_fname, avg_time_sec=time_average_sec_dict[band], average_using="difmap", working_dir=None)
        fig = radplot(df, label="Data")
        plt.show()


if __name__ == "__main__":

    # Generation of real data set
    # # band_uvfits_files_dict = {"c1": "/home/ilya/Downloads/MF/0851+202.c1.2009_02_02.uvf",
    # #                           "c2": "/home/ilya/Downloads/MF/0851+202.c2.2009_02_02.uvf",
    # #                           "k1": "/home/ilya/Downloads/MF/0851+202.k1.2009_02_02.uvf",
    # #                           "q1": "/home/ilya/Downloads/MF/0851+202.q1.2009_02_02.uvf",
    # #                           "u1": "/home/ilya/Downloads/MF/0851+202.u1.2009_02_02.uvf",
    # #                           "x1": "/home/ilya/Downloads/MF/0851+202.x1.2009_02_02.uvf",
    # #                           "x2": "/home/ilya/Downloads/MF/0851+202.x2.2009_02_02.uvf"}
    # band_uvfits_files_dict = {"c1": "/home/ilya/Downloads/BK150_mdf/0827+243.c1.2009_04_09.puvf",
    #                           "c2": "/home/ilya/Downloads/BK150_mdf/0827+243.c2.2009_04_09.uvf_difmap",
    #                           "k1": "/home/ilya/Downloads/BK150_mdf/0827+243.k1.2009_04_09.puvf0",
    #                           "q1": "/home/ilya/Downloads/BK150_mdf/0827+243.q1.2009_04_09.uvf_difmap",
    #                           "u1": "/home/ilya/Downloads/BK150_mdf/0827+243.u1.2009_04_09.puvf0",
    #                           "x1": "/home/ilya/Downloads/BK150_mdf/0827+243.x1.2009_04_09.puvf",
    #                           "x2": "/home/ilya/Downloads/BK150_mdf/0827+243.x2.2009_04_09.puvf0"}
    # time_average_sec_dict = {"c1": 120, "c2": 120, "x1": 120, "x2": 120, "u1": 120, "k1": 120, "q1": 120}
    # save_dir = "/home/ilya/github/bam/mf/0827"
    # if not os.path.exists(save_dir):
    #     os.mkdir(save_dir)
    # create_data_files_from_real_uvfits(band_uvfits_files_dict, save_dir, time_average_sec_dict=time_average_sec_dict)
    #
    # sys.exit(0)


    simulate = True
    time_average_sec = 120
    save_dir = "/home/sonya/bam/Release"
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    # Coordinates relative to the true jet origin
    # a, p, c, PA, size_1GHz, k_r, S_1GHz, alpha
    core_component = (3.0, 1.5, 1.0, np.pi/6, 1.0, 1.0, 2.0, 0.0)
    # RA, DEC, Size, nu_max, S_nu_max, alpha_thick, alpha_thin
    jet_components = [(3.0, 4.5, 0.5, 10.0, 1.0, 1.5, -0.5),
                      (9.0, 10.0, 1.5, 5.0, 0.5, 2.0, -0.5)]
    band_uvfits_files_dict = {"c1": "data/0851+202.c1.2009_02_02.uvf",
                              # "c2": "/home/ilya/Downloads/MF/0851+202.c2.2009_02_02.uvf",
                              "k1": "data/0851+202.k1.2009_02_02.uvf",
                              "q1": "data/0851+202.q1.2009_02_02.uvf",
                              "u1": "data/0851+202.u1.2009_02_02.uvf",
                              "x1": "data/0851+202.x1.2009_02_02.uvf"}
                              # "x2": "/home/ilya/Downloads/MF/0851+202.x2.2009_02_02.uvf"}

    for band, uvfits in band_uvfits_files_dict.items():
        total_flux = 0.0
        uvdata = UVData(uvfits)
        freq_ghz = uvdata.frequency/1E+09
        out_fname = os.path.join(save_dir, "{}_{:.2f}".format(band, freq_ghz))
        # df = create_data_file_v2(uvfits, out_fname, time_average_sec=120)

        path, fname = os.path.split(uvfits)
        to_read_uvfits = os.path.join(path, "ta{}sec_{}".format(time_average_sec, fname))
        # save_uvfits = to_read_uvfits
        time_average(uvfits, to_read_uvfits, time_sec=time_average_sec)
        df = load_from_fits(to_read_uvfits)


        if simulate:
            u = df["u"].values
            v = df["v"].values
            uv = np.vstack((u, v)).T

            df["vis_re"] = 0
            df["vis_im"] = 0
            for gc in jet_components:
                print(gc)
                RA, DEC, Size, nu_max, S_nu_max, alpha_thick, alpha_thin = gc
                flux = S_nu_max*(freq_ghz/nu_max)**alpha_thick/(1 - np.exp(-1)) * (1 - np.exp(-(freq_ghz/nu_max)**(alpha_thin-alpha_thick)))
                print("Jet Component Flux at nu = {:.2f} is S = {:.4f}".format(freq_ghz, flux))
                # Add model
                re, im = gaussian_circ_ft(flux=flux, RA=RA, DEC=DEC, bmaj=Size, uv=uv)
                df["vis_re"] += re
                df["vis_im"] += im
                total_flux += flux

            # Core
            a, p, c, PA, size_1GHz, k_r, S_1GHz, alpha = core_component
            flux = S_1GHz*freq_ghz**alpha
            print("Core Flux at nu = {:.2f} is S = {:.4f}".format(freq_ghz, flux))

            def func(t, c, a, freq_ghz, k_r):
                return c * (np.arcsinh(t) + t * np.sqrt(np.power(t, 2) + 1)) - a * np.power(freq_ghz, -1/k_r)
            def func_prime(t, c):
                return 2 * c * np.sqrt(pow(t, 2) + 1)
            root = newton(lambda t: func(t, c, a, freq_ghz, k_r), x0=0, fprime=lambda t: func_prime(t, c), maxiter=30)
            Ra_before_rotation = p * np.power(root, 2)
            Dec_before_rotation = 2 * p * root
            RA = Ra_before_rotation * np.cos(PA) - Dec_before_rotation * np.sin(PA)
            DEC = Ra_before_rotation * np.sin(PA) + Dec_before_rotation * np.cos(PA)


            total_flux += flux
            Size = size_1GHz*freq_ghz**(-1/k_r)
            re, im = gaussian_circ_ft(flux=flux, RA=RA, DEC=DEC, bmaj=Size, uv=uv)
            df["vis_re"] += re
            df["vis_im"] += im

            print("Total flux at frequency {:.2f} is S = {:.2f}".format(freq_ghz, total_flux))

            # center_mass /= total_flux
            # print("Center mas (RA, DEC) [mas] = ", center_mass)
            # Shift to bring center mass to the phase center
            re = df["vis_re"]
            im = df["vis_im"]
            vis = re + 1j*im
            # "-" becaue we want to move in the opposite direction
            shift = [-RA*mas_to_rad, -DEC*mas_to_rad]
            result = np.exp(2.0*np.pi*1j*(uv @ shift))
            # print("In shifting center mass uv = ", uv)
            vis *= result
            df["vis_re"] = np.real(vis)
            df["vis_im"] = np.imag(vis)

        try:
            fig = radplot(df, label="Data")
        except AttributeError:
            fig = radplot(df, label="Data", style="reim")

        # Add noise and plot
        if simulate:
            df_updated = add_noise(df, use_global_median_noise=True, global_noise_scale=1)
            fig = radplot(df_updated, color="#ff7f0e", fig=fig, label="With noise")
        else:
            df_updated = df
        df_updated.to_csv(out_fname, sep=" ", index=False, header=False)
        # save_uvfits = os.path.join(save_dir, "0851+202.{}.2009_02_02.uvf".format(band))
        # save_to_uvfits(save_uvfits, re=df_updated["vis_re"], im=df_updated["vis_im"])
        os.unlink(to_read_uvfits)
