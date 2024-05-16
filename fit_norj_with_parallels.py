import pandas as pd
import sys
import os
import glob
import pathlib
from data_utils import get_data_file_from_ehtim
from tempfile import TemporaryDirectory


def reformat_table(original_csv):
    df = pd.read_csv(original_csv)
    csv_dir, _ = os.path.split(original_csv)

    uvfits_files = list()
    ccfits_files = list()
    ncomps = list()
    base_UVFITS_dir = "/mnt/jet1/yyk/VLBI/RFC/images"
    for row in df.itertuples():
        source = row.source
        obs1 = row.obs1
        obs2 = row.obs2
        uvfits = f"{base_UVFITS_dir}/{source}/{source}_{obs1}_vis.fits"
        if not os.path.exists(uvfits):
            raise Exception(f"No UVFITS {uvfits}")
        ccfits = f"{base_UVFITS_dir}/{source}/{source}_{obs1}_map.fits"
        if not os.path.exists(ccfits):
            raise Exception(f"No CCFITS {ccfits}")
        uvfits_files.append(uvfits)
        ccfits_files.append(ccfits)
        ncomps.append(2)

        uvfits = f"{base_UVFITS_dir}/{source}/{source}_{obs2}_vis.fits"
        if not os.path.exists(uvfits):
            raise Exception(f"No UVFITS {uvfits}")
        ccfits = f"{base_UVFITS_dir}/{source}/{source}_{obs2}_map.fits"
        if not os.path.exists(ccfits):
            raise Exception(f"No CCFITS {ccfits}")
        uvfits_files.append(uvfits)
        ccfits_files.append(ccfits)
        ncomps.append(2)

    dic = {"uvfits": uvfits_files, "ccfits": ccfits_files, "ncomps": ncomps}
    df = pd.DataFrame.from_dict(dic)
    df.to_csv(os.path.join(csv_dir, "full_table_ncomp_2.csv"), header=True, index=False)
    return df


def choose_maxnsaves(n_comps, n_vis):
    if n_vis < 100:
        maxnsaves = 3000
    elif 100 <= n_vis < 300:
        maxnsaves = 5000
    elif 300 <= n_vis < 600:
        maxnsaves = 7000
    else:
        maxnsaves = 10000
    return int(n_comps*maxnsaves/2)


def run_parallels_on_df(df, base_save_dir, executable_dict, template_options, n_jobs, difmap_avg_time_sec, dry_run=False):



    uvfits_files = list()
    ccfits_file = list()
    n_components = list()
    results_dirs = list()
    basenames = list()
    data_files = list()
    maxnsaves_list = list()
    executables = list()
    failed_uvfits = list()
    for row in df.itertuples():
        uvfits = row.uvfits
        ccfits = row.ccfits
        ncomps = row.ncomps


        uvfits_fn = os.path.split(uvfits)[-1]
        source, band, year, month, day, author, product = uvfits_fn.split("_")
        product = product.split(".")[0]
        assert product == "vis"
        save_basename = f"{source}_{band}_{year}_{month}_{day}_{author}_ncomp_{ncomps}"
        results_dir = os.path.join(base_save_dir, f"{source}_{band}_{year}_{month}_{day}_{author}")
        pathlib.Path(results_dir).mkdir(parents=True, exist_ok=True)
        out_fname = os.path.join(results_dir, f"{source}_{band}_{year}_{month}_{day}_{author}_avg_{difmap_avg_time_sec}sec.csv")
        with TemporaryDirectory() as working_dir:
            try:
                df = get_data_file_from_ehtim(uvfits, out_fname, avg_time_sec=difmap_avg_time_sec, working_dir=working_dir)
            except:
                signal_file = os.path.join(base_save_dir, f"FAILED_CREATING_DATA_FILE_{source}_{band}_{year}_{month}_{day}_{author}.info")
                os.system(f"touch {signal_file}")
                failed_uvfits.append(uvfits)
                continue

        n_vis = len(df)
        print("Number of vis = ", n_vis)
        print("Number of comps = ", ncomps)
        maxnsaves = choose_maxnsaves(ncomps, n_vis)
        print("Chosen maxnsaves = ", maxnsaves)
        # Append only if its OK with data file!
        basenames.append(save_basename)
        results_dirs.append(results_dir)
        maxnsaves_list.append(str(maxnsaves))
        n_components.append(ncomps)
        data_files.append(out_fname)

        try:
            executables.append(executable_dict[ncomps])
        except KeyError:
            raise Exception(f"No executable for #components = {ncomps}!")


    with open(os.path.join(base_save_dir, "failed_UVFITS.txt"), "w") as fo:
        for fn in failed_uvfits:
            fo.write(fn + "\n")

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

    # 536 is OK!
    # 609 failed!
    print(f"Length of the argument is {len(maxnsaves_list)}")

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
    df_log = pd.read_csv(f"{base_save_dir}/joblog", sep="\t")
    return df_log
# This works in bash
# $ parallel -k  --link --dryrun "python fit_mojave_uvfits_single_ncomponents.py --uvfits \"uvfits\" --results_dir \"results_dir\" --ncomps {1} --maxnsaves {2}" ::: 1 2 ::: 1000 2000


# Path to UVFITS, path to FITS, n_components
table_file = "/home/ilya/data/VLBI_Gaia/full_table_ncomp_2.csv"
n_jobs = 10
# Average time in difmap. It uses weighted vector averaging and
# also re-calculate weights. This is helpful when they are not reliable.
difmap_avg_time_sec = 60
# Directory to save the results
base_save_dir = "/mnt/storage/ilya/VLBI_Gaia"
# FIXME: Choose executable depending on ``component_type``
executable_dict = {2: "/home/ilya/github/bam/Release/bam_2"}
                   #3: "/home/ilya/github/bam/Release/bam_3"}
template_options = "/home/ilya/github/bam/OPTIONS_gen"

logs_list = list() 
df_full = pd.read_csv(table_file)
# df = df_full.head(600)
df = df_full[:100]
df_log = run_parallels_on_df(df, base_save_dir, executable_dict, template_options, n_jobs, difmap_avg_time_sec, dry_run=False)
logs_list.append(df_log)
