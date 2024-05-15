import pandas as pd
import sys
import os
import glob
import pathlib
from data_utils import get_data_file_from_ehtim


def choose_maxnsaves(n_comps, n_vis):
    if n_vis < 100:
        maxnsaves = 3000
    elif 100 <= n_vis < 300:
        maxnsaves = 5000
    elif 300 <= n_vis < 600:
        maxnsaves = 7000
    else:
        maxnsaves = 10000
    return int(n_comps*maxnsaves/2/100)


dry_run = False
# table_file = "/home/ilya/data/VLBI_Gaia/obs_pairs_good_8.csv"
table_file = "/home/ilya/data/VLBI_Gaia/test.csv"
# Path to UVFITS, path to FITS, n_components
df = pd.read_csv(table_file)
base_UVFITS_dir = "/mnt/jet1/yyk/VLBI/RFC/images"
# Directory to save the results
base_save_dir = "/mnt/storage/ilya/VLBI_Gaia"
n_jobs = 2

# FIXME: Choose executable depending on ``component_type``
executable_dict = {2: "/home/ilya/github/bam/Release/bam_2",
                   3: "/home/ilya/github/bam/Release/bam_3"}
executable = "/home/ilya/github/bam/Release/bam"

# Average time in difmap. It uses weighted vector averaging and
# also re-calculate weights. This is helpful when they are not reliable.
# FIXME: Use function from create_data_from_ehtim
difmap_avg_time_sec = 60
template_options = "/home/ilya/github/bam/OPTIONS_gen"


uvfits_files = list()
ccfits_file = list()
n_components = list()
results_dirs = list()
basenames = list()
data_files = list()
maxnsaves_list = list()
executables = list()
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
    results_dirs.append(results_dir)
    pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
    out_fname = os.path.join(results_dir, f"{save_basename}.csv")
    df = get_data_file_from_ehtim(uvfits, out_fname)
    n_vis = len(df)
    print("Number of vis = ", n_vis)
    print("Number of comps = ", ncomps)
    maxnsaves = choose_maxnsaves(ncomps, n_vis)
    print("Chosen maxnsaves = ", maxnsaves)
    maxnsaves_list.append(str(maxnsaves))
    n_components.append(ncomps)
    try:
        executables.append(executable_dict[ncomps])
    except KeyError:
        raise Exception(f"No executable for #components = {ncomps}!")
    data_files.append(out_fname)


args = " ".join(["{:s} {:s} {:s}".format(i, j, k) for (i, j, k) in zip(basenames, data_files, results_dirs)])

maxnsaves = " ".join(maxnsaves_list)
executables = " ".join(executables)
basenames = " ".join(basenames)
data_files = " ".join(data_files)
results_dirs = " ".join(results_dirs)

print("basenames\n", basenames)
print("data_files\n", data_files)
print("results_dirs\n", results_dirs)
print("maxnsaves\n", maxnsaves)
print("executables\n", executables)
print("args\n", args)

if dry_run:
    os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --dryrun --link '
              '"python fit_mojave_uvfits_single_ncomponents.py --template_options {} --basename {{1}} '
        '--data_file {{2}} --results_dir {{3}} --maxnsaves {{4}} --executable {{5}}" ::: {} ::: {} ::: {} ::: {} ::: {}'.format(base_save_dir, base_save_dir, n_jobs,
                                                                                   template_options,
                                                                                   basenames, data_files, results_dirs, maxnsaves, executables))
else:
    os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --link '
              '"python fit_mojave_uvfits_single_ncomponents.py --template_options {} --basename {{1}} '
        '--data_file {{2}} --results_dir {{3}} --maxnsaves {{4}} --executable {{5}}" ::: {} ::: {} ::: {} ::: {} ::: {}'.format(base_save_dir, base_save_dir, n_jobs,
                                                                                   template_options,
                                                                                   basenames, data_files, results_dirs, maxnsaves, executables))

# This works in bash
# $ parallel -k  --link --dryrun "python fit_mojave_uvfits_single_ncomponents.py --uvfits \"uvfits\" --results_dir \"results_dir\" --ncomps {1} --maxnsaves {2}" ::: 1 2 ::: 1000 2000
