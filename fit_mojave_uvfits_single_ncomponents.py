import os
import sys
import pathlib
from string import Template
import numpy as np
import argparse
sys.path.insert(0, '/home/ilya/github/dnest4post')
from postprocess import postprocess


def run(basename, maxnsaves, data_file, results_dir, template_options, excecutable):

    run_name = basename
    result_options = os.path.join(results_dir, "OPTIONS_{}".format(run_name))
    options = {'nparticles': 1,
               'newlevel': 30000,
               'saveinterval': 30000,
               'threadsteps': 100,
               'maxnlevels': 0,
               'lambda': 30,
               'beta': 100,
               'maxnsaves': maxnsaves,
               'samplefile': '{}/sample_{}.txt'.format(results_dir, run_name),
               'sampleinfo': '{}/sample_info_{}.txt'.format(results_dir, run_name),
               'levelsfile': '{}/levels_{}.txt'.format(results_dir, run_name)}

    filein = open(template_options)
    src = Template(filein.read())
    result = src.substitute(options)
    with open(result_options, "w") as fo:
        print(result, file=fo)

    executable_dir = os.path.split(executable)[0]
    os.chdir(executable_dir)
    os.system("{} -t 4 -o {} -d {}".format(excecutable,
                                           result_options,
                                           data_file))
    print("POSTPROCESSING DNEST OUTPUT...")
    post_samples_file = os.path.join(results_dir, "posterior_sample_{}.txt".format(run_name))

    logZ, _, _ = postprocess(plot=True, pics_prefix=f"{results_dir}/{basename}", show_pics=False,
                             sample_file='{}/sample_{}.txt'.format(results_dir, run_name),
                             level_file='{}/levels_{}.txt'.format(results_dir, run_name),
                             sample_info_file='{}/sample_info_{}.txt'.format(results_dir, run_name),
                             post_sample_file=post_samples_file)
    os.system("touch {}/{}_logZ_{:.3f}.info".format(results_dir, run_name, logZ))
    os.system("rm {}/sample_{}.txt".format(results_dir, run_name))
    os.system("rm {}/levels_{}.txt".format(results_dir, run_name))
    os.system("rm {}/sample_info_{}.txt".format(results_dir, run_name))

    post_samples = np.loadtxt(post_samples_file)
    n_post_samples = len(post_samples)
    if n_post_samples < 100:
        os.system("touch {}/WARNING_SMALL_SAMPLE_SIZE_{}".format(results_dir, run_name))


if __name__ == "__main__":

    CLI = argparse.ArgumentParser()
    CLI.add_argument("--basename",
                     type=str)
    CLI.add_argument("--data_file",
                     type=str)
    CLI.add_argument("--maxnsaves",
                     type=int)
    CLI.add_argument("--results_dir",
                     type=str)
    CLI.add_argument("--executable",
                     type=str)
    CLI.add_argument("--template_options",
                     type=str)
    args = CLI.parse_args()

    basename = args.basename
    data_file = args.data_file
    maxnsaves = args.maxnsaves
    results_dir = args.results_dir
    executable = args.executable
    template_options = args.template_options

    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)

    print("Basename = ", basename)
    print("Data file = ", data_file)
    print("Results dir = ", results_dir)
    run(basename, maxnsaves, data_file, results_dir, template_options, executable)
