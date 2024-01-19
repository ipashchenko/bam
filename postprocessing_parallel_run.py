import os
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


component_type = "cg"
n_max = 20
n_max_samples_to_plot = 1000
# Plot all samples - for easy handling component cluster membership
jitter_first = True
skip_hp = True
if component_type == "cg":
    comp_length = 4
elif component_type == "eg":
    comp_length = 6
else:
    raise Exception
pixsize_mas = 0.1
freq_ghz = 15.4


for uvfits, original_ccfits, data_file, results_dir in zip():
    df = pd.read_csv(data_file)
    posterior_file = os.path.join(results_dir, "posterior_sample.txt")
    save_dir = results_dir
    save_rj_ncomp_distribution_file = os.path.join(save_dir, "ncomponents_distribution.png")
    # FIXME: Find from data file
    n_antennas = 10
    n_jitters = n_antennas

    posterior_samples = np.loadtxt(posterior_file)

    save_basename = os.path.split(uvfits)[-1].split(".uvf")[0]


    plot_per_antenna_jitters_and_offsets(posterior_samples, uvfits=uvfits, save_dir=save_dir,
                                         save_basename=save_basename, plot_offsets=False,
                                         n_antennas=n_antennas)

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
