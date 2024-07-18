import matplotlib.pyplot as plt
from postprocessing_mf_utils import convert_posterior_file_to_pandas_df
import corner
import numpy as np

df = convert_posterior_file_to_pandas_df("/home/sonya/bam/saved_data/94_posterior_sample.txt")
# print(df['a'])
a = 3
p = 0.7
PA = 0
PA_1 = np.pi/6
c = 0.1
nu = 4.6
k_r = 1.5
# print(len(df['p']))
# array_p = np.full((len(df['p']), 1), 1.5)

# data = np.vstack([array_p, df["p"]]).T
figure = corner.corner(df[['a', 'k_r', 'p', 'c', 'PA']], labels=['a', 'k_r', 'p', 'c', 'PA'], show_titles=True)
# plt.hist(df.iloc[:, 22], bins=30)
# figure = corner.corner(data, bins=30, labels=['a', 'p'], show_titles=True)
plt.show()