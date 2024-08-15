from matplotlib import pyplot as plt
from postprocess import postprocess
import os
import time
import shutil
import glob


path = input('Enter base path: ')
if path == 'home':
    path = '/home/sonya/bam/'
elif path == 'server':
    path = '/home/skomkova/'

os.chdir('Release')
files = glob.glob(path + 'saved_data/*')
for f in files:
    os.remove(f)
i = 0
while True:
    i += 1
    postprocess(pics_prefix=path + f'saved_data/{i}', show_pics=False)
    shutil.copy(path + 'Release/posterior_sample.txt', path + f'saved_data/{i}_posterior_sample.txt')
    shutil.copy(path + 'Release/sample.txt', path + f'saved_data/{i}_sample.txt')
    time.sleep(180)
    plt.close()