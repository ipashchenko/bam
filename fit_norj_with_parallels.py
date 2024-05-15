import pandas as pd
import sys
import os
import glob
import pathlib
from data_utils import get_data_file_from_ehtim


dry_run = True
# table_file = "/home/ilya/data/VLBI_Gaia/obs_pairs_good_8.csv"
table_file = "/home/ilya/data/VLBI_Gaia/test.csv"
# Path to UVFITS, path to FITS, n_components
df = pd.read_csv(table_file)
base_UVFITS_dir = "/mnt/jet1/yyk/VLBI/RFC/images"
# Directory to save the results
base_save_dir = "/mnt/storage/ilya/VLBI_Gaia"
n_jobs = 5
maxnsaves = 1000

# FIXME: Choose executable depending on ``component_type``
executable_dict = {1: "/home/ilya/github/bam/Release/bam_1",
                   2: "/home/ilya/github/bam/Release/bam_2",
                   3: "/home/ilya/github/bam/Release/bam_3",
                   4: "/home/ilya/github/bam/Release/bam_4",
                   5: "/home/ilya/github/bam/Release/bam_5",
                   6: "/home/ilya/github/bam/Release/bam_6",
                   7: "/home/ilya/github/bam/Release/bam_7"}
executable = "/home/ilya/github/bam/Release/bam"

# Average time in difmap. It uses weighted vector averaging and
# also re-calculate weights. This is helpful when they are not reliable.
difmap_avg_time_sec = 60
template_options = "/home/ilya/github/bam/OPTIONS_gen"


uvfits_files = list()
ccfits_file = list()
n_components = list()
results_dirs = list()
basenames = list()
out_fnames = list()
for row in df.itertuples():
    uvfits = row.uvfits
    ccfits = row.ccfits
    ncomps = row.ncomps


    uvfits_fn = os.path.split(uvfits)[-1]
    source, band, year, month, day, author, product = uvfits_fn.split("_")
    product = product.split(".")[0]
    assert product == "vis"
    save_basename = f"{source}_{band}_{year}_{month}_{day}_{author}"
    basenames.append(save_basename)
    results_dir = os.path.join(base_save_dir, save_basename)
    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
    out_fname = os.path.join(results_dir, f"{save_basename}.csv")
    df = get_data_file_from_ehtim(uvfits, out_fname)
    out_fnames.append(out_fname)


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
              '"python fit_mojave_uvfits_single_ncomponents.py --executable {} --maxnsaves {} --template_options {} --basename {{1}} '
              '--data_file {{2}} --results_dir {{3}}" ::: {} ::: {} ::: {}'.format(results_dir, results_dir, n_jobs,
                                                                                   executable, maxnsaves, template_options,
                                                                                   basenames, out_fnames, results_dirs))
else:
    os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --link '
              '"python fit_mojave_uvfits_single_ncomponents.py --executable {} --maxnsaves {} --template_options {} --basename {{1}} '
              '--data_file {{2}} --results_dir {{3}}" ::: {} ::: {} ::: {}'.format(results_dir, results_dir, n_jobs,
                                                                                   executable, maxnsaves, template_options,
                                                                                   basenames, out_fnames, results_dirs))

# This works in bash
# $ parallel -k  --link --dryrun "python fit_mojave_uvfits_single_ncomponents.py --uvfits \"uvfits\" --results_dir \"results_dir\" --ncomps {1} --maxnsaves {2}" ::: 1 2 ::: 1000 2000
