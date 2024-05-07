import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pickle
from postprocessing_utils import (plot_per_antenna_jitters_and_offsets, rj_plot_ncomponents_distribution,
                                  clusterization, plot_model_predictions, get_samples_for_each_n, sort_samples_by_r,
                                  plot_position_posterior, plot_position_posterior_clustered)
import sys
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from from_fits import create_clean_image_from_fits_file
from image import plot as iplot
from spydiff import find_bbox, find_image_std, import_difmap_model


def get_n_antennas_from_data_file(data_file):
    df = pd.read_csv(data_file)
    return len(set(list(df.t1.values) + list(df.t2.values)))


component_type = "eg"
n_max = 20
n_max_samples_to_plot = 1000
# Plot all samples - for easy handling component cluster membership
jitter_first = False
skip_hp = True
if component_type == "cg":
    comp_length = 4
elif component_type == "eg":
    comp_length = 6
else:
    raise Exception
pixsize_mas = 0.1
freq_ghz = 15.4

source = "0136+176"
base_save_dir = "/home/ilya/data/rjbam"
uvfits_dir = "/home/ilya/Downloads/mojave/0136+176/uvfits_dir"
ccfits_dir = "/home/ilya/Downloads/mojave/0136+176/ccfits_dir"
uvfits_files = None
ccfits_files = None
data_files = None

if uvfits_files is None:
    uvfits_files = sorted(glob.glob(os.path.join(uvfits_dir, f"ta60_{source}*.uvf")))

epochs = list()
for uvfits_file in uvfits_files:
    _, fn = os.path.split(uvfits_file)
    epoch = fn.split(".")[2]
    epochs.append(epoch)

print(epochs)

if ccfits_files is None:
    ccfits_files = list()
    for epoch in epochs:
        ccfits_file = os.path.join(ccfits_dir, f"{source}.u.{epoch}.icn.fits")
        ccfits_files.append(ccfits_file)

if data_files is None:
    data_files = list()
    for epoch in epochs:
        data_file = os.path.join(base_save_dir, source, epoch, f"{source}.u.{epoch}.csv")
        data_files.append(data_file)

results_dirs = list()
for epoch in epochs:
    results_dir = os.path.join(base_save_dir, source, epoch)
    results_dirs.append(results_dir)

print(uvfits_files)
print("===============")
print(ccfits_files)
print("===============")
print(data_files)
print("===============")
print(results_dirs)
print("===============")

# sys.exit(0)

for epoch, uvfits, original_ccfits, data_file, results_dir in zip(epochs, uvfits_files, ccfits_files, data_files, results_dirs):
    df = pd.read_csv(data_file)
    posterior_file = os.path.join(results_dir, f"posterior_sample_{source}.u.{epoch}.txt")
    save_dir = results_dir
    save_rj_ncomp_distribution_file = os.path.join(save_dir, "ncomponents_distribution.png")
    # FIXME: Find from data file

    n_antennas = get_n_antennas_from_data_file(data_file)
    if jitter_first:
        n_jitters = n_antennas
    else:
        n_jitters = 0

    posterior_samples = np.loadtxt(posterior_file)
    print("Read {:d} posterior samples".format(len(posterior_samples)))

    save_basename = f"{source}.u.{epoch}"


    if jitter_first:
        plot_per_antenna_jitters_and_offsets(posterior_samples, uvfits=uvfits, save_dir=save_dir,
                                             save_basename=save_basename, plot_offsets=False,
                                             n_antennas=n_antennas)

    fig = rj_plot_ncomponents_distribution(posterior_file, picture_fn=save_rj_ncomp_distribution_file,
                                           jitter_first=jitter_first, n_jitters=n_jitters, type=component_type,
                                           normed=True, show=False, skip_hyperparameters=skip_hp)
    plt.close()
    # samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
    #                                             n_jitters=n_jitters, n_max=n_max,
    #                                             skip_hyperparameters=skip_hp,
    #                                             type=component_type)

    # sys.exit(0)

    # labels_dict = clusterization(posterior_file, jitter_first=jitter_first, n_jitters=n_jitters, n_max=n_max,
    #                              skip_hyperparameters=skip_hp, component_type=component_type,
    #                              algorithm="hdbscan",
    #                              dbscan_min_core_samples_frac_of_posterior_size=0.3,
    #                              dbscan_eps=0.2,
    #                              hdbscan_min_cluster_frac_of_posterior_size=0.5)

    # sys.exit(0)

    plot_model_predictions(posterior_file, data_file, rj=True, n_jitters=n_jitters, component_type=component_type,
                           style="ap", n_samples_to_plot=1000, alpha_model=0.01, skip_hyperparameters=skip_hp,
                           savefname=os.path.join(save_dir, "radplot.png"), show=False)
    plt.close()



    # sys.exit(0)

    samples_for_each_n = get_samples_for_each_n(posterior_samples, jitter_first,
                                                n_jitters=n_jitters, n_max=n_max,
                                                skip_hyperparameters=skip_hp,
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

        # fig_out = plot_position_posterior_clustered(samples_to_plot,
        #                                             savefn=None, ra_lim=None, dec_lim=None,
        #                                             difmap_model_fn=None, type=component_type, s=1.0, figsize=None,
        #                                             sorted_componets=False, fig=fig_p2,
        #                                             inverse_xaxis=False,
        #                                             alpha_opacity=0.3,
        #                                             cluster_membership=labels_dict[n_component])
        # fig_out.savefig(os.path.join(save_dir, f"CLEAN_image_ncomp_{n_component}_clusters.png"), dpi=600)
        # plt.close(fig_out)

        # f = plot_size_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :],
        #                                  savefn=os.path.join(save_dir, f"r_R_ncomp_{n_component}.png"),
        #                                  type=component_type)
        # plt.close(f)
        #
        # f = plot_tb_distance_posterior(samples_to_plot[:n_max_samples_to_plot, :], freq_ghz, type=component_type,
        #                                savefn=os.path.join(save_dir, f"r_Tb_ncomp_{n_component}.png"))
        # plt.close(f)
