import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.optimize import curve_fit
from matplotlib.animation import FuncAnimation

def func(t, c, a, nu, k_r):
    return c * (np.arcsinh(t) + t * np.sqrt(np.power(t, 2) + 1)) - a * np.power(nu, -1/k_r)
def func_prime(t, c):
    return 2 * c * np.sqrt(pow(t, 2) + 1)
def parabola(x, a, b, c):
    return a * x**2 + b * x + c

# a = 1.0
PA = np.pi/6
c = 1.0
# nu = 4.6
k_r = 1.0

results = []

for a in np.linspace(-10, 10, 10):
    Ra = []
    Dec = []
    for nu in [4.6, 5.0, 8.1, 8.43, 15.4, 23.8, 43.2]:
        root = newton(lambda t: func(t, c, a, nu, k_r), x0=0, fprime=lambda t: func_prime(t, c), maxiter=30)
        Ra_before_rotation = a * np.power(root, 2)
        Dec_before_rotation = 2 * a * root
        Ra.append(Ra_before_rotation * np.cos(PA) - Dec_before_rotation * np.sin(PA))
        Dec.append(Ra_before_rotation * np.sin(PA) + Dec_before_rotation * np.cos(PA))
    # print(f"freq: {nu}, root: {root}, Ra: {RA}, Dec: {DEC}")
    popt, _ = curve_fit(parabola, Ra, Dec)
    results.append((Ra, Dec, popt))

fig, ax = plt.subplots()
colors = plt.cm.viridis(np.linspace(0, 1, len(results)))
sc = ax.scatter([], [], color='blue')
line, = ax.plot([], [], lw=2)
ax.grid()

def init():
    sc.set_offsets([])
    line.set_data([], [])
    return sc, line

def update(frame):
    Ra, Dec, popt = results[frame]
    sc.set_offsets(np.c_[Ra, Dec])
    x = np.linspace(min(Ra), max(Ra), 100)
    y = parabola(x, *popt)
    line.set_data(x, y)
    line.set_color(colors[frame])
    sc.set_color(colors[frame])
    return sc, line

ani = FuncAnimation(fig, update, frames=len(results), init_func=init, blit=True, repeat=False)
plt.show()