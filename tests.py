import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.optimize import curve_fit

def func(t, c, a, nu, k_r):
    return c * (np.arcsinh(t) + t * np.sqrt(np.power(t, 2) + 1)) - a * np.power(nu, -1/k_r)
def func_prime(t, c):
    return 2 * c * np.sqrt(pow(t, 2) + 1)

a = 10
PA = np.pi/6
c = 1.0
# nu = 4.6
k_r = 1.0

Ra = []
Dec = []

for nu in [4.6, 5.0, 8.1, 8.43, 15.4, 23.8, 43.2]:
    root = newton(lambda t: func(t, c, a, nu, k_r), x0=0, fprime=lambda t: func_prime(t, c), maxiter=30)
    Ra_before_rotation = a * np.power(root, 2)
    Dec_before_rotation = 2 * a * root
    Ra.append(Ra_before_rotation * np.cos(PA) - Dec_before_rotation * np.sin(PA))
    Dec.append(Ra_before_rotation * np.sin(PA) + Dec_before_rotation * np.cos(PA))
    # print(f"freq: {nu}, root: {root}, Ra: {RA}, Dec: {DEC}")

def parabola(x, a, b, c):
    return a * x**2 + b * x + c

popt, pcov = curve_fit(parabola, Ra, Dec)
print(popt)

# plt.plot(Ra, parabola(np.array(Ra), *popt))
# plt.grid()
# plt.show()