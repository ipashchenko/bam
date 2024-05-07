import sys

import numpy as np
from scipy.linalg import eigh, cholesky
from scipy.stats import norm
from pylab import plot, show, axis, subplot, xlabel, ylabel, grid
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import corner
from tqdm import tqdm

num_samples = 10000


class Gaussian2D(object):
    def __init__(self, mu_1, mu_2, rho, sigma_1, sigma_2):
        self.mu_1 = mu_1
        self.mu_2 = mu_2
        self.rho = rho
        self.sigma_1 = sigma_1
        self.sigma_2 = sigma_2

    def cdf(self, x):
        norm_1 = norm(loc=self.mu_1, scale=self.sigma_1)
        u_1 = norm_1.cdf(x[0])
        norm_2 = norm(loc=self.mu_2 + self.rho*self.sigma_2/self.sigma_1*(x[0] - self.mu_1),
                      scale=np.sqrt(1 - self.rho**2)*self.sigma_2)
        u_2 = norm_2.cdf(x[1])
        return np.array([u_1, u_2])

    def cdf_inverse(self, u):
        norm_1 = norm(loc=self.mu_1, scale=self.sigma_1)
        x_1 = norm_1.ppf(u[0])
        norm_2 = norm(loc=self.mu_2 + self.rho*self.sigma_2/self.sigma_1*(x_1 - self.mu_1),
                      scale=np.sqrt(1 - self.rho**2)*self.sigma_2)
        x_2 = norm_2.ppf(u[1])
        return np.array([x_1, x_2])

    def generate(self):
        u = np.random.uniform(0, 1, size=2)
        return self.cdf_inverse(u)


# g2d = Gaussian2D(-1.74, -0.71, -0.8, 1.62, 1.03)
g2d = Gaussian2D(-1.74, -0.71, -0.8, 0.75*1.62, 0.75*1.03)
samples = list()
u_samples = list()
for i in tqdm(range(num_samples)):
    sample = g2d.generate()
    u_sample = g2d.cdf(sample)
    samples.append(sample)
    u_samples.append(u_sample)
samples = np.atleast_2d(samples)
u_samples = np.atleast_2d(u_samples)
#
#
# corner.corner(samples, labels=["logflux", "logsize"])
# plt.show()
#
# corner.corner(u_samples, labels=["u_logflux", "u_logsize"])
# plt.show()


# sys.exit(0)
method = 'cholesky'

# r = np.array([
#         [  1.00, -0.50],
#         [ -0.50,  1.00]
#         ])


logflux_orig = np.loadtxt("logflux.txt")
logsize_orig = np.loadtxt("logsize.txt")
fig, axes = plt.subplots(1,1)
axes.scatter(logflux_orig, logsize_orig, alpha=0.2, s=3)


# logflux_samples, logsize_samples = np.loadtxt("Release/samples.txt", unpack=True)
logflux_samples = samples[:, 0]
logsize_samples = samples[:, 1]
axes.scatter(logflux_samples, logsize_samples, alpha=0.2, s=3)
axes.set_xlabel("logflux")
axes.set_ylabel("logsize")
plt.show()

plt.hist(logflux_samples, bins=40)
plt.xlabel("flux")
plt.show()

plt.hist(logsize_samples, bins=40)
plt.xlabel("size")
plt.show()

# logflux_u_samples, logsize_u_samples = np.loadtxt("Release/u_samples.txt", unpack=True)
# x = np.vstack((logflux_u_samples, logsize_u_samples)).T
# corner.corner(x, labels=["logflux", "logsize"])
# plt.show()

sys.exit(0)

# import corner
# x = np.vstack((logflux_orig, logsize_orig)).T
# corner.corner(x, labels=["flux", "size"])
# plt.show()

gmm = GaussianMixture(1)
a = np.vstack((logflux_orig, logsize_orig)).T
gmm.fit(a)

# import sys
# sys.exit(0)

r_orig = np.array([[2.63036257, -0.69296081],
                   [-0.69296081,  1.06533235]])

mean_logflux_orig = -1.73741931
mean_logsize_orig = -0.71104578
sigma_logflux_orig = np.sqrt(2.63036257)
sigma_logsize_orig = np.sqrt(1.06533235)
ro_orig = -0.69296081
sigma_logsize_cond_orig = np.sqrt((1 - ro_orig ** 2)) * sigma_logsize_orig


mean_logflux = mean_logflux_orig
mean_logsize = mean_logsize_orig
sigma_logflux = sigma_logflux_orig
sigma_logsize = sigma_logsize_orig
ro = ro_orig
sigma_logsize_cond = 1.13*sigma_logsize_cond_orig

# N(-1.74, 1.62)
# logflux = np.random.normal(mean_logflux, sigma_logflux, 3000)
logflux = np.random.normal(-1.73, 1.12, num_samples)
# logflux = np.random.normal(-2.29, 1.62, num_samples)
# N(-0.71-0.44*logflux, 0.84)
# logsize = np.random.normal(mean_logsize + ro * sigma_logsize / sigma_logflux * logflux, sigma_logsize_cond, 3000)
# logsize = np.random.normal(-0.71 -0.44 * logflux, 0.84, 3000)
logsize = np.random.normal(-1.71 - 0.5 * logflux, 0.44, num_samples)

# axes.scatter(logflux_orig, logsize_orig, alpha=0.05)
axes.scatter(logflux, logsize, alpha=0.01, color="C2")
axes.set_xlabel("logflux")
axes.set_xlim([-7, 5])
axes.set_ylim([-7, 5])
axes.set_aspect('equal')
plt.show()

plt.hist(logflux, bins=40)
plt.xlabel("flux")
plt.show()

plt.hist(logsize, bins=40)
plt.xlabel("size")
plt.show()




# x = norm.rvs(mean_logflux_orig, 1, size=(2, num_samples))
# c = cholesky(r_orig, lower=True)
# y = np.dot(c, x)
#
# subplot(1,1,1)
# plot(y[0], y[1], 'b.')
# ylabel('y[1]')
# xlabel('y[0]')
# axis('equal')
# grid(True)
#
# show()