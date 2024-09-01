import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from postprocessing_mf_utils import convert_posterior_file_to_pandas_df

def func(t, c, a, nu, k_r):
    return c * (np.arcsinh(t) + t * np.sqrt(np.power(t, 2) + 1)) - a * np.power(nu, -1/k_r)
def func_prime(t, c):
    return 2 * c * np.sqrt(pow(t, 2) + 1)

df = convert_posterior_file_to_pandas_df("/home/sonya/bam/saved_data/28_posterior_sample.txt")
a = 3
a_1 = 3.22
PA = np.pi/6
PA_1 = 0.54
c = 0.1
c_1 = 0.86
# nu = 4.6
k_r = 1.0
k_r_1 = 0.95

Ra = []; Dec = []
Ra1 = []; Dec1 = []

for nu in [4.6, 8.1, 15.4, 23.79]:
    root = newton(lambda t: func(t, c, a, nu, k_r), x0=0, fprime=lambda t: func_prime(t, c), maxiter=30)
    print(root)
    Ra_before_rotation = c * np.power(root, 2)
    Dec_before_rotation = 2 * c * root
    Ra.append(Ra_before_rotation * np.cos(PA) - Dec_before_rotation * np.sin(PA))
    Dec.append(Ra_before_rotation * np.sin(PA) + Dec_before_rotation * np.cos(PA))
    # print(f"freq: {nu}, root: {root}, Ra: {RA}, Dec: {DEC}")

for nu in [4.6, 8.1, 15.4, 23.79]:
    root = newton(lambda t: func(t, c_1, a_1, nu, k_r_1), x0=0, fprime=lambda t: func_prime(t, c_1), maxiter=30)
    print(root)
    Ra_before_rotation_1 = c_1 * np.power(root, 2)
    Dec_before_rotation_1 = 2 * c_1 * root
    Ra1.append(- Ra_before_rotation_1 * np.cos(PA_1) + Dec_before_rotation_1 * np.sin(PA_1))
    Dec1.append(- Ra_before_rotation_1 * np.sin(PA_1) - Dec_before_rotation_1 * np.cos(PA_1))

# print(Ra, Dec)
# print(Ra1, Dec1)



plt.xlabel('Ra')
plt.ylabel('Dec')
fig, axes = plt.subplots(1,1)
axes.set_aspect('equal')
plt.scatter(df.iloc[:, 4], df.iloc[:, 8], color='blue')
plt.scatter(df.iloc[:, 5], df.iloc[:, 9], color='red')
plt.scatter(df.iloc[:, 6], df.iloc[:, 10], color='black')
plt.scatter(df.iloc[:, 7], df.iloc[:, 11], color='orange')
axes.scatter(0, 0, color='black')
axes.scatter(Ra, Dec, color='blue')
axes.scatter(Ra1, Dec1, color='green')
plt.grid()
plt.show()