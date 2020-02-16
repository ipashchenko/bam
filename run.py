import glob
import os
from string import Template
import pandas as pd


def convert_mojave_epoch_inverse(epoch):
    year = epoch.split('_')[0]
    month = epoch.split('_')[1]
    day = epoch.split('_')[2]
    return "{}-{}-{}".format(year, month, day)


n_components = None
source = "0716+714"
ncomp_file = "/home/ilya/github/bam/data/0716+714_comp_number.txt"
df = pd.read_csv(ncomp_file, sep="\s+", names=["source", "epoch", "ncomp"])
df_source = df.query("source == @source")
ncomp_med = int(df_source.ncomp.median())
data_files = sorted(glob.glob(os.path.join("/home/ilya/github/bam/data/{}".format(source),
                                             "{}_*.txt".format(source))))
for data_file in data_files:
    fname = os.path.split(data_file)[-1]
    epoch = fname[9:19]
    epoch_inv = convert_mojave_epoch_inverse(epoch)
    df_epoch = df_source.query("epoch == @epoch_inv")
    if df_epoch.empty:
        ncomp = ncomp_med
    else:
        ncomp = df_epoch.ncomp.values[0]
    run_name = "{}_{}_ncomp{}".format(source, epoch, ncomp)
    template_options = "/home/ilya/github/bam/OPTIONS_gen"
    result_options = "/home/ilya/github/bam/Release/OPTIONS_{}".format(run_name)
    options = {'nparticles': 1,
               'newlevel': 30000,
               'saveinterval': 30000,
               'threadsteps': 100,
               'maxnlevels': 0,
               'lambda': 30,
               'beta': 100,
               'maxnsaves': 10000,
               'samplefile': 'sample_{}.txt'.format(run_name),
               'sampleinfo': 'sample_info_{}.txt'.format(run_name),
               'levelsfile': 'levels_{}.txt'.format(run_name)}

    excecutables = {i: "./bam{}".format(i) for i in range(1, 11)}

    filein = open(template_options)
    src = Template(filein.read())
    result = src.substitute(options)
    with open(result_options, "w") as fo:
        print(result, file=fo)

    import os
    os.chdir("/home/ilya/github/bam/Release")
    # os.system("{} -t 4 -o {} -d {} >/dev/null 2>&1 &".format(excecutables[n_components],
    #                                          result_options,
    #                                          data_file))
    os.system("{} -t 5 -o {} -d {}".format(excecutables[ncomp],
                                                        result_options,
                                                        data_file))
    # print("Running {} -t 4 -o {} -d {}".format(excecutables[ncomp],
    #                                                     result_options,
    #                                                     data_file))