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
    ccfits_files = list()
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
                signal_file_home = os.path.join(results_dir, f"FAILED_CREATING_DATA_FILE_{source}_{band}_{year}_{month}_{day}_{author}.info")
                os.system(f"touch {signal_file}")
                os.system(f"touch {signal_file_home}")
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
        ccfits_files.append(ccfits)

        try:
            executables.append(executable_dict[ncomps])
        except KeyError:
            raise Exception(f"No executable for #components = {ncomps}!")


    with open(os.path.join(base_save_dir, "failed_UVFITS.txt"), "w") as fo:
        for fn in failed_uvfits:
            fo.write(fn + "\n")

    maxnsaves_file = "maxnsaves.data"
    with open(maxnsaves_file, "w") as fo:
        for item in maxnsaves_list:
            fo.write(item + "\n")
    executables_file = "executables.data"
    with open(executables_file, "w") as fo:
        for item in executables:
            fo.write(item + "\n")
    basenames_file = "basenames.data"
    with open(basenames_file, "w") as fo:
        for item in basenames:
            fo.write(item + "\n")
    data_files_file = "data_files.data"
    with open(data_files_file, "w") as fo:
        for item in data_files:
            fo.write(item + "\n")
    results_dirs_file = "results_dirs.data"
    with open(results_dirs_file, "w") as fo:
        for item in results_dirs:
            fo.write(item + "\n")
    ccfits_files_file = "cfits_files.data"
    with open(ccfits_files_file, "w") as fo:
        for item in ccfits_files:
            fo.write(item + "\n")

    if dry_run:
        os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --dryrun --link '
                  '--xapply "python fit_mojave_uvfits_single_ncomponents.py --template_options {} --basename {{1}} '
            '--data_file {{2}} --results_dir {{3}} --maxnsaves {{4}} --executable {{5}}" :::: {} :::: {} :::: {} :::: {} :::: {}'.format(base_save_dir, base_save_dir, n_jobs,
                                                                                       template_options,
                                                                                       basenames_file, data_files_file, results_dirs_file, maxnsaves_file, executables_file))
    else:
        os.system('parallel --files --results {}/res_{{1}} --joblog {}/joblog --jobs {} --link '
                  '--xapply "python fit_mojave_uvfits_single_ncomponents.py --template_options {} --basename {{1}} '
            '--data_file {{2}} --results_dir {{3}} --maxnsaves {{4}} --executable {{5}} --ccfits_file {{6}}" :::: {} :::: {} :::: {} :::: {} :::: {} :::: {}'.format(base_save_dir, base_save_dir, n_jobs,
                                                                                       template_options,
                                                                                       basenames_file, data_files_file, results_dirs_file, maxnsaves_file, executables_file, ccfits_files_file))
    df_log = pd.read_csv(f"{base_save_dir}/joblog", sep="\t")
    return df_log
# This works in bash
# $ parallel -k  --link --dryrun "python fit_mojave_uvfits_single_ncomponents.py --uvfits \"uvfits\" --results_dir \"results_dir\" --ncomps {1} --maxnsaves {2}" ::: 1 2 ::: 1000 2000


if __name__ == "__main__":
    # Path to UVFITS, path to FITS, n_components
    table_file = "/home/ilya/data/VLBI_Gaia/full_table_ncomp_2.csv"
    n_jobs = 15
    # Average time in difmap. It uses weighted vector averaging and
    # also re-calculate weights. This is helpful when they are not reliable.
    difmap_avg_time_sec = 60
    # Directory to save the results
    base_save_dir = "/mnt/storage/ilya/VLBI_Gaia"
    # FIXME: Choose executable depending on ``component_type``
    executable_dict = {2: "/home/ilya/github/bam/Release/bam_2"}
                       #3: "/home/ilya/github/bam/Release/bam_3"}
    template_options = "/home/ilya/github/bam/OPTIONS_gen"

    df_full = pd.read_csv(table_file)
    # df = df_full.head(600)
    df = df_full
    df_log = run_parallels_on_df(df, base_save_dir, executable_dict, template_options, n_jobs, difmap_avg_time_sec, dry_run=False)
