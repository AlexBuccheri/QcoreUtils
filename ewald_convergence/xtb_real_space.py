""" Look at the modified xTB Ewald potential """
import numpy as np
from scipy.special import erfc
import matplotlib.pyplot as plt

def xtb(r, eta):
    return np.sqrt(1 / (r * r + eta ** (-2)))

def delta_xtb(r, eta):
    return xtb(r, eta) - 1/np.abs(r)

def ewald(r, alpha):
    return erfc(np.sqrt(alpha) * np.abs(r)) / np.abs(r)

eta = 10
alpha = 1

r_values = np.linspace(0.001, 2.001, 500)
one_over_r = [1/r for r in r_values]
long_range = [1-ewald(r, alpha) for r in r_values]
ewald_real = [ewald(r, alpha) for r in r_values]
ewald_delta_xtb = [delta_xtb(r, eta) for r in r_values]
xtb_standard = [xtb(r, eta) for r in r_values]

xtb_real = np.asarray(ewald_real) + np.asarray(ewald_delta_xtb)

plt.plot(r_values, one_over_r, label="1/r")
plt.plot(r_values, ewald_real, label="Screened standard Ewald")
plt.plot(r_values, xtb_real, label="Xtb Correction + Screened Ewald")
plt.plot(r_values, long_range, label="Long-range part of Ewald")
plt.plot(r_values, ewald_delta_xtb, label="xTB Ewald Correction")
plt.plot(r_values, xtb_standard, label="Xtb sqrt(1/|r|^2+eta^-2)")


plt.legend()
plt.ylim(-10, 20)
plt.xlim(0, 1)
plt.show()

plt.plot(r_values, one_over_r, label="1/r")
plt.plot(r_values, ewald_real, label="Screened standard Ewald")
plt.plot(r_values, xtb_real, label="Xtb Correction + Screened Ewald")
plt.plot(r_values, xtb_standard, label="Xtb sqrt(1/|r|^2+eta^-2)")
plt.legend()
plt.ylim(-10, 20)
plt.xlim(0, 1)
plt.show()
plt.show()
""" Plot shows: If one cuts off the """