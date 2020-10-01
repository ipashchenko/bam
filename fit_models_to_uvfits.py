import os
import sys
from string import Template
import numpy as np
from convert_uvfits_to_datafile import create_data_file, radplot, plot_model_predictions
sys.path.insert(0, '/home/ilya/github/dnest4post')
from postprocessing_utils import process_norj_samples
from postprocess import postprocess
import matplotlib.pyplot as plt


# uvfits = "/home/ilya/github/bam/data/0716+714/0716+714.u.2006_12_01.uvf"
uvfits = "/home/ilya/data/zhenya/1652+398.u1.2009_05_14.uvf_difmap"
time_average_sec = 120
n_components_min = 6
n_components_max = 6
# results_dir = "/home/ilya/github/bam/data/0716+714"
results_dir = "/home/ilya/data/zhenya/u1"

# For large models could be increased
maxnsaves = 4000


################# No need to change this #################
# How sort the posterior samples of components? By distance from phase center ("r"), by flux ("flux"), by RA ("ra") or
# by DEC ("dec")? Probably "r" will work. If not - try "flux".
sort_by = "r"
executables_dir = "/home/ilya/github/bam/Release"
excecutables = {i: "./bam{}".format(i) for i in range(1, 11)}
uvfits_fname = os.path.split(uvfits)[-1]
uvfits_basename = ".".join(uvfits_fname.split(".")[:-1])
data_file = os.path.join(results_dir, "{}.txt".format(uvfits_basename))
df = create_data_file(uvfits, outfile=data_file, time_average_sec=time_average_sec)
fig = radplot(df, style="reim", savefname=os.path.join(results_dir, "radplot_{}.png".format(uvfits_basename)), show=False)

n_components = range(n_components_min, n_components_max+1)

for ncomp in n_components:
    run_name = "{}_ncomp{}".format(uvfits_basename, ncomp)
    template_options = "/home/ilya/github/bam/OPTIONS_gen"
    result_options = os.path.join(results_dir, "OPTIONS_{}".format(run_name))
    options = {'nparticles': 1,
               'newlevel': 10000,
               'saveinterval': 10000,
               'threadsteps': 100,
               'maxnlevels': 0,
               'lambda': 30,
               'beta': 100,
               'maxnsaves': maxnsaves,
               'samplefile': 'sample_{}.txt'.format(run_name),
               'sampleinfo': 'sample_info_{}.txt'.format(run_name),
               'levelsfile': 'levels_{}.txt'.format(run_name)}

    filein = open(template_options)
    src = Template(filein.read())
    result = src.substitute(options)
    with open(result_options, "w") as fo:
        print(result, file=fo)

    os.chdir(executables_dir)
    os.system("{} -t 4 -o {} -d {}".format(excecutables[ncomp],
                                           result_options,
                                           data_file))
    print("POSTPROCESSING DNEST OUTPUT...")
    post_samples_file = os.path.join(results_dir, "posterior_sample_{}.txt".format(run_name))
    logZ, _, _ = postprocess(plot=False,
                             sample_file='sample_{}.txt'.format(run_name),
                             level_file='levels_{}.txt'.format(run_name),
                             sample_info_file='sample_info_{}.txt'.format(run_name),
                             post_sample_file=post_samples_file)
    os.system("touch {}/{}_logZ_{:.3f}.info".format(results_dir, run_name, logZ))
    os.system("rm sample_{}.txt".format(run_name))
    os.system("rm levels_{}.txt".format(run_name))
    os.system("rm sample_info_{}.txt".format(run_name))

    post_samples = np.loadtxt(post_samples_file)
    n_post_samples = len(post_samples)
    if n_post_samples < 100:
        os.system("touch {}/WARNING_SMALL_SAMPLE_SIZE_FOR_NCOMP{}_{}".format(results_dir, ncomp, run_name))

    fig = plot_model_predictions(post_samples, df, savefname=os.path.join(results_dir, "radplot_model_{}.png".format(run_name)),
                                 show=False)
    plt.close(fig)
    figs = process_norj_samples(post_samples_file, jitter_first=True, ra_lim=(-2, 2), dec_lim=(-2, 2),
                                freq_ghz=15.4, z=0.0, sort_by=sort_by,
                                savefn_position_post=os.path.join(results_dir, "position_posterior_{}.png".format(run_name)),
                                savefn_fluxsize_isot_post=os.path.join(results_dir, "flux_size_isoTb_{}.png".format(run_name)),
                                savefn_fluxsize_post=os.path.join(results_dir, "flux_size_{}.png".format(run_name)),
                                savefn_rtb_post=os.path.join(results_dir, "r_tb_{}.png".format(run_name)),
                                savefn_sizer_post=os.path.join(results_dir, "r_size_{}.png".format(run_name)),
                                savefn_sizetb_post=os.path.join(results_dir, "size_tb_{}.png".format(run_name)),
                                savefn_corner=os.path.join(results_dir, "corner_{}.png".format(run_name)),
                                dfm_outname=os.path.join(results_dir, "median_model_{}.dfm".format(run_name)))
    for fig in figs:
        plt.close(fig)
