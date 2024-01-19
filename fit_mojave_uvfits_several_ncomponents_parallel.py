import sys
import os
import glob
from data_utils import get_data_file_from_ehtim


dry_run = False
source = "0136+176"
# cg, eg, sp. Not used now.
component_type = "cg"
# Directory with UVFITS files for this source in MOJAVE format
uvfits_dir = "/home/ilya/Downloads/mojave/0136+176/uvfits_dir"
# Or iterable of UVFITS files
uvfits_files = None
# Directory to save the results
base_save_dir = "/home/ilya/data/rjbam"
n_jobs = 5

# FIXME: Choose executable depending on ``component_type``
executable = "/home/ilya/github/bam/Release/bam"

maxnsaves = 70000
# Average time in difmap. It uses weighted vector averaging and also re-calculate weightsю This is helpful when they are
# not reliable.
difmap_avg_time_sec = 60
template_options = "/home/ilya/github/bam/OPTIONS_gen"

if uvfits_files is None:
    uvfits_files = glob.glob(os.path.join(uvfits_dir, f"{source}*.uvf"))

epochs = list()
for uvfits_file in uvfits_files:
    _, fn = os.path.split(uvfits_file)
    epoch = fn.split(".")[2]
    epochs.append(epoch)

results_dir = os.path.join(base_save_dir, source)
if not os.path.exists(results_dir):
    os.mkdir(results_dir)
results_dirs = list()
for epoch in epochs:
    results_dir_epoch = os.path.join(results_dir, epoch)
    if not os.path.exists(results_dir_epoch):
        os.mkdir(results_dir_epoch)
    results_dirs.append(results_dir_epoch)

out_fnames = list()
basenames = list()
for uvfits_file, epoch in zip(uvfits_files, epochs):
    basename = f"{source}.u.{epoch}"
    basenames.append(basename)
    out_fname = os.path.join(base_save_dir, source, epoch, f"{basename}.csv")
    out_fnames.append(out_fname)
    df = get_data_file_from_ehtim(uvfits_file, out_fname, avg_time_sec=difmap_avg_time_sec, average_using="difmap")

args = " ".join(["{:s} {:s} {:s}".format(i, j, k) for (i, j, k) in zip(basenames, out_fnames, results_dirs)])

basenames = " ".join(basenames)
out_fnames = " ".join(out_fnames)
results_dirs = " ".join(results_dirs)

print("basenames\n", basenames)
print("outfnames\n", out_fnames)
print("results_dirs\n", results_dirs)
print("args\n", args)

if dry_run:
    os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --dryrun --link '
              '"python fit_mojave_uvfits_single_ncomponents.py --executable {} --template_options {} --basename {{1}} '
              '--data_file {{2}} --results_dir {{3}}" ::: {} ::: {} ::: {}'.format(results_dir, results_dir, n_jobs,
                                                                                   executable, template_options,
                                                                                   basenames, out_fnames, results_dirs))
else:
    os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --link '
              '"python fit_mojave_uvfits_single_ncomponents.py --executable {} --template_options {} --basename {{1}} '
              '--data_file {{2}} --results_dir {{3}}" ::: {} ::: {} ::: {}'.format(results_dir, results_dir, n_jobs,
                                                                                   executable, template_options,
                                                                                   basenames, out_fnames, results_dirs))

# This works in bash
# $ parallel -k  --link --dryrun "python fit_mojave_uvfits_single_ncomponents.py --uvfits \"uvfits\" --results_dir \"results_dir\" --ncomps {1} --maxnsaves {2}" ::: 1 2 ::: 1000 2000
