import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.optimize import curve_fit
from matplotlib.widgets import Button, Slider

def func(t, c, a, nu, k_r):
    return c * (np.arcsinh(t) + t * np.sqrt(np.power(t, 2) + 1)) - a * np.power(nu, -1/k_r)
def func_prime(t, c):
    return 2 * c * np.sqrt(pow(t, 2) + 1)

t = np.linspace(-10, 10, 1000)
a = 10
PA = np.pi/6
c = 1.0
nu = 4.6
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

# Create the figure and the line that we will manipulate
fig, ax = plt.subplots()
line, = ax.plot(t, parabola(np.array(Ra), *popt), lw=2)
ax.set_xlabel('Ra')

# adjust the main plot to make room for the sliders
fig.subplots_adjust(left=0.25, bottom=0.25)

# Make a horizontal slider to control parameter a.
axa = fig.add_axes([0.25, 0.1, 0.65, 0.03])
freq_slider = Slider(
    ax=axa,
    label='Parameter a',
    valmin=-10,
    valmax=10,
    valinit=a,
)

# Make a vertically oriented slider to control the amplitude
axPA = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
amp_slider = Slider(
    ax=axPA,
    label="PA",
    valmin=0,
    valmax=20,
    valinit=PA,
    orientation="vertical"
)

# The function to be called anytime a slider's value changes
def update(val):
    line.set_ydata(parabola(np.array(Ra), *popt))
    fig.canvas.draw_idle()

# plt.plot(Ra, parabola(np.array(Ra), *popt))
# plt.grid()
# plt.show()