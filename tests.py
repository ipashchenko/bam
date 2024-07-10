# import matplotlib.pyplot as plt
# import numpy as np


# def func_1(t, nu, k_r):
#     return (np.arcsinh(t) + t*np.sqrt(t**2+1)) - np.power(nu, -1/k_r)

# t = np.linspace(-10, 10, 1000)
# plt.plot(t, func_1(t, 4.6*10**9, 5))
# plt.show()

from sympy.solvers import solve, sqrt, arcsinh
from sympy import Symbol
import numpy as np

t = Symbol('t')
solve(arcsinh(t) + t*sqrt(t**2 + 1))