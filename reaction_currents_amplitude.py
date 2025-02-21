import numpy as np
import matplotlib.pyplot as plt
from functions import wong
from functions import findmax
from functions import calc_nu
from functions import coleColeFit as ccf

w = np.logspace(-4, 4, 1600)

sig0 = 11.5e-3
F = 96485.309
mu = 5e-8
c0 = sig0/2/mu/F

T = 293
epsr = 80
a = 4.75e-3
alpha = 1e-10
beta = 1e-2
nu = calc_nu(-21e-3, 26e-3, 0, 24e-3, a, 102e-3)

sig = []

sigmax = []
tau = []
taucc = []
c = []


c3 = np.logspace(-4, -1, 100)
for c3_temp in c3:
    sig_temp = wong(w, c0, c0*c3_temp, T, epsr, mu, a, alpha, beta, nu)
    
    sigmax_temp, tau_temp = findmax(w, sig_temp.imag)
    sigmax.append(sigmax_temp)
    tau.append(tau_temp)
    
    sol = ccf(w, sig_temp)
    taucc.append(sol[2])
    c.append(sol[3])
    
    sig.append(sig_temp)

kB = 8.6173e-5
eps0 = 8.85e-12
D = mu*kB*T
lambdaD = 1/np.sqrt(2*c0*F/(eps0*epsr*kB*T))
tau_buecker_vd = a**2/(4*D*(1-2*(c0/c3))**2)
tau_buecker_dl = a*lambdaD/(2*D)*c3/c3

plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplot_mosaic([["specre", "specre", "amplitude"],
                                ["specre", "specre", "amplitude"],
                                ["specre", "specre", "tau"],
                                ["specim", "specim", "tau"],
                                ["specim", "specim", "c"],
                                ["specim", "specim", "c"]], figsize=(13,11))

ax["specre"].set_title("(a)", loc="left")
ax["specim"].set_title("(b)", loc="left")
ax["amplitude"].set_title("(c)", loc="left")
ax["tau"].set_title("(d)", loc="left")
ax["c"].set_title("(e)", loc="left")

ax["amplitude"].semilogx(c3, sigmax, "-", linewidth=2)
ax["amplitude"].set_xlabel("$c_3 / c_0$ [-]")
ax["amplitude"].set_ylabel("$\u03C3_{max}''/\u03C3_w$ [-]")

ax["tau"].loglog(c3, taucc, "-", linewidth=2)
ax["tau"].loglog(c3, tau_buecker_dl, "--", linewidth=2)
ax["tau"].loglog(c3, tau_buecker_vd, "-.", linewidth=2)
ax["tau"].set_xlabel("$c_3 / c_0$ [-]")
ax["tau"].set_ylabel("$\u03C4$ [s]")
ax["tau"].legend(["Num. model", "Eq. (25)", "Eq. (26)"])
ax["tau"].set_ylim([5e-3, 2e1])

ax["c"].semilogx(c3, c, "-", linewidth=2)
ax["c"].set_xlabel("$c_3 / c_0$ [-]")
ax["c"].set_ylabel("$c$ [-]")

col = plt.cm.viridis(np.linspace(0.1, 1, 10))

ii = 0
jj = 0

for sigma in sig:
    if ii % 10 == 0:
        ax["specre"].semilogx(w, sigma.real, linewidth=2, color=col[jj])
        ax["specim"].semilogx(w, sigma.imag, linewidth=2, color=col[jj])
        jj += 1
    ii += 1
ax["specim"].legend(["$c_3 = $0.1 mmol/m³", "$c_3 = $0.2 mmol/m³",
                      "$c_3 = $0.4 mmol/m³", "$c_3 = $0.8 mmol/m³",
                      "$c_3 = $1.6 mmol/m³", "$c_3 = $3.3 mmol/m³",
                      "$c_3 = $6.6 mmol/m³", "$c_3 = $13.2 mmol/m³",
                      "$c_3 = $26.6 mmol/m³", "$c_3 = $53.4 mmol/m³"],
                     loc="upper left")
ax["specre"].set_xlabel("$\u03C9$ [rad/s]")
ax["specre"].set_ylabel("$\u03C3'/\u03C3_w$ [-]")
ax["specim"].set_xlabel("$\u03C9$ [rad/s]")
ax["specim"].set_ylabel("$\u03C3''/\u03C3_w$ [-]")
fig.tight_layout()

fig.savefig("figure_07.pdf", dpi=300)