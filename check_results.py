import os
import sys
sys.path.insert(0, '/home/ilya/github/dnest4post')
from postprocess import postprocess


source = "0552+398"
source = sys.argv[2]
epoch = "1995_01_20"
epoch = sys.argv[3]
n_components = sys.argv[1]
run_name = "{}_{}_ncomp{}".format(source, epoch, n_components)

os.system("rsync -P --rsh=ssh ilya@193.232.12.59:/home/ilya/github/bam/Release/sample_info_{}.txt $PWD".format(run_name))
os.system("rsync -P --rsh=ssh ilya@193.232.12.59:/home/ilya/github/bam/Release/sample_{}.txt $PWD".format(run_name))
os.system("rsync -P --rsh=ssh ilya@193.232.12.59:/home/ilya/github/bam/Release/levels_{}.txt $PWD".format(run_name))
os.system("cp /home/ilya/github/bam/sample_info_{}.txt sample_info.txt".format(run_name))
os.system("cp /home/ilya/github/bam/sample_{}.txt sample.txt".format(run_name))
os.system("cp /home/ilya/github/bam/levels_{}.txt levels.txt".format(run_name))
postprocess()
os.system("rm sample_info.txt levels.txt log_prior_weights.txt weights.txt")