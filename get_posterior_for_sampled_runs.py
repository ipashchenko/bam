import os
import sys
import glob
sys.path.insert(0, '/home/ilya/github/dnest4post')
from postprocess import postprocess


samples_dir = "/home/ilya/github/bam/Release"
source = "2200+420"

sample_files = glob.globS(os.path.join(samples_dir, "sample_{}_*.txt".format(source)))
epochs = list()
for sample_file in sample_files:
    fn = os.path.split(sample_file)[-1]
    epoch = fn[16:26]
    epochs.append(epoch)


for epoch in epochs:
    sample_file = glob.glob(os.path.join(samples_dir, "sample_{}_{}_*txt".format(source, epoch)))[0]
    level_file = glob.glob(os.path.join(samples_dir, "levels_{}_{}_*txt".format(source, epoch)))[0]
    sample_info_file = glob.glob(os.path.join(samples_dir, "sample_info_{}_{}_*txt".format(source, epoch)))[0]
    print(sample_file, sample_info_file, level_file)
    assert sample_file[-5:-4] == sample_info_file[-5:-4] == level_file[-5:-4]
    n_comp = sample_file[-5:-4]
    logZ, _, _ = postprocess(plot=False,
                             sample_file=sample_file,
                             level_file=level_file,
                             sample_info_file=sample_info_file,
                             post_sample_file=os.path.join(samples_dir, "posterior_sample_{}_{}_ncomp_{}.txt".format(source, epoch, n_comp)))