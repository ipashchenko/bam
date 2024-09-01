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

folder = input('Enter folder name where you want to save your files in: ')
if folder == 'saved_data_1':
    folder = 'saved_data_1/'
elif folder == 'saved_data_2':
    folder = 'saved_data_2/'
elif folder == 'saved_data_3':
    folder = 'saved_data_3/'
elif folder == 'saved_data_4':
    folder = 'saved_data_4/'
elif folder == 'saved_data_5':
    folder = 'saved_data_5/'

os.chdir('Release')
files = glob.glob(path + folder + '*')
for f in files:
    os.remove(f)
i = 0
while True:
    i += 1
    postprocess(pics_prefix=path + folder + f'{i}', show_pics=False)
    shutil.copy(path + 'Release/posterior_sample.txt', path + folder + f'{i}_posterior_sample.txt')
    shutil.copy(path + 'Release/sample.txt', path + folder + f'{i}_sample.txt')
    time.sleep(180)
    plt.close()