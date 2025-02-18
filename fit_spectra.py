import numpy as np
import matplotlib.pyplot as plt
from functions import wong
from functions import fit_c3
from functions import calc_nu
from functions import wong_eff

dat = np.loadtxt("1Kugel_25_corr_vers2.txt")
sig0 = 11.5e-3

w_exp = 2*np.pi*np.flip(dat[18:58,0])
sig_exp = (np.flip(dat[18:58,1])+1j*np.flip(dat[18:58,2]))

form = sig0*1e3/dat[61,1]/1.0092384355718929

F = 96485.309
mu = 5e-8
c0 = sig0/2/mu/F
a = 4.75e-3
nu = calc_nu(-21e-3, 26e-3, 0, 24e-3, a, 102e-3)
poro = 0.4015
T = 293
epsr = 80
alpha = 1e-10
beta = 1e-2*1000

w_mod = 2*np.pi*np.logspace(-2, 3, 500)

c3_m1 = fit_c3(c0, T, epsr, mu, a, alpha, beta, nu, w_exp, sig_exp*form*1e-3)
sig_m1, sig_m2 = wong_eff(w_mod, c0, c3_m1, T, epsr, mu, a, alpha, beta, nu,
                          dat[18:58,0]*2*np.pi, dat[18:58,1]+1j*dat[18:58,2],
                          mod="sig")/form

c3_m3 = fit_c3(c0*poro, T, epsr, mu/form/poro, a, alpha, beta, nu, w_exp,
               sig_exp*1e-3)
sig_m3 = wong(w_mod, c0*poro, c3_m3, T, epsr, mu/form/poro, a, alpha, beta, nu,
              mod="sig")

plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplots(2, 2, figsize=(12,9))
ax[0,0].semilogx(dat[:,0], dat[:,1], "o", markersize=5)
ax[1,0].semilogx(dat[:,0], dat[:,2], "o", markersize=5)
ax[0,1].semilogx(dat[:,0], dat[:,1], "o", markersize=5)
ax[1,1].semilogx(dat[:,0], dat[:,2], "o", markersize=5)

ax[0,0].semilogx(w_mod/2/np.pi, 1e3*sig_m1.real, "-.", linewidth=2.5,
                 color="black")
ax[1,0].semilogx(w_mod/2/np.pi, 1e3*sig_m1.imag, "-.", linewidth=2.5,
                 color="black")

ax[0,1].semilogx(w_mod/2/np.pi, 1e3*sig_m2.real, "-", linewidth=2.5,
                 color="tab:red")
ax[1,1].semilogx(w_mod/2/np.pi, 1e3*sig_m2.imag, "-", linewidth=2.5,
                 color="tab:red")
ax[0,1].semilogx(w_mod/2/np.pi, 1e3*sig_m3.real, "--", linewidth=2.5,
                 color="tab:orange")
ax[1,1].semilogx(w_mod/2/np.pi, 1e3*sig_m3.imag, "--", linewidth=2.5,
                 color="tab:orange")

ax[0,0].set_xlim(1e-2, 1e3)
ax[1,0].set_xlim(1e-2, 1e3)
ax[0,1].set_xlim(1e-2, 1e3)
ax[1,1].set_xlim(1e-2, 1e3)
ax[0,0].set_ylim(3.3, 3.5)
ax[1,0].set_ylim(0, 0.04)
ax[0,1].set_ylim(3.3, 3.5)
ax[1,1].set_ylim(0, 0.04)

ax[0,0].legend(["Measured data", "Original response, divided by $F$"], loc=9)
ax[0,1].legend(["Measured data",
                "Response with app. mobility, divided by $F$",
                "Response with turtuosity"], loc=9)

ax[0,0].set_xlabel("$\u03C9$ [rad/s]")
ax[0,1].set_xlabel("$\u03C9$ [rad/s]")
ax[1,0].set_xlabel("$\u03C9$ [rad/s]")
ax[1,1].set_xlabel("$\u03C9$ [rad/s]")
ax[0,0].set_ylabel("$\u03C3'$ [mS/m]")
ax[0,1].set_ylabel("$\u03C3'$ [mS/m]")
ax[1,0].set_ylabel("$\u03C3''$ [mS/m]")
ax[1,1].set_ylabel("$\u03C3''$ [mS/m]")

ax[0,0].set_title("(a)", loc="left")
ax[1,0].set_title("(b)", loc="left")
ax[0,1].set_title("(c)", loc="left")
ax[1,1].set_title("(d)", loc="left")

fig.tight_layout()
fig.savefig("figure_09", dpi=300)