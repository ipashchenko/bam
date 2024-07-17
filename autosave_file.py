from matplotlib import pyplot as plt
from postprocess import postprocess
import os
import time
import shutil
import glob

os.chdir('Release')
files = glob.glob('/home/sonya/bam/saved_data/*')
for f in files:
    os.remove(f)
i = 0
while True:
    i += 1
    postprocess(pics_prefix=f'/home/sonya/bam/saved_data/{i}', show_pics=False)
    shutil.copy('/home/sonya/bam/Release/posterior_sample.txt', f'/home/sonya/bam/saved_data/{i}_posterior_sample.txt')
    shutil.copy('/home/sonya/bam/Release/sample.txt', f'/home/sonya/bam/saved_data/{i}_sample.txt')
    time.sleep(180)
    plt.close()