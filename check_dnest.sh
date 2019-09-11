#!/bin/bash
rsync -P --rsh=ssh ilya@193.232.12.59:/home/ilya/github/bam/Release/\{sample_info.txt,sample.txt,levels.txt\} $PWD
python -c "import os; cwd = os.getcwd(); os.chdir(\"/home/ilya/github/dnest4post\"); from postprocess import postprocess; os.chdir(cwd); postprocess()"
rm sample_info.txt levels.txt log_prior_weights.txt weights.txt