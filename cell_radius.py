import numpy as np
import matplotlib.pyplot as plt
from functions import load_data as ld
from functions import coleColeFit as ccf

def gurin(nu_list):
    m = []
    
    for nu in nu_list:
        m_temp = 1 - (2*(1-nu)**2) / ((2+nu) * (1 + 2*nu))
        m.append(m_temp)
        
    return m

def revil(nu_list):
    m = []
    
    for nu in nu_list:
        m_temp = 9/2 * nu
        m.append(m_temp)
        
    return m

plt.rcParams.update({'font.size': 13})

erg_nonreg = []
erg_reg = []

fname = ["8e-3", "10e-3", "12e-3", "14e-3", "16e-3", "18e-3", "20e-3", "22e-3",
         "24e-3", "26e-3", "28e-3", "30e-3", "32e-3"]

R = []

for name in fname:
    R_temp = float(name)
    R.append(R_temp)
    erg_nonreg.append(ld("cell_size", "R_"+name, -0.03, 0.03, 0, R_temp,
                         4.75e-3, 102e-3))
    erg_reg.append(ld("cell_size", "R_"+name, -0.03, 0.03, 0, R_temp, 4.75e-3,
                      102e-3, "reg"))
    
w_nonreg = []
sig_nonreg = []
m_nonreg = []
tau_nonreg = []
nu_nonreg = []

for obj in erg_nonreg:
    w_nonreg.append(obj[0])
    sig_nonreg.append(obj[1])
    nu_nonreg.append(obj[4])
    sol = ccf(obj[0], obj[1])
    m_nonreg.append(sol[1])
    tau_nonreg.append(sol[2])

w_reg = []
sig_reg = []
m_reg = []
tau_reg = []

for obj in erg_reg:
    w_reg.append(obj[0])
    sig_reg.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_reg.append(sol[1])
    tau_reg.append(sol[2])

col = plt.cm.viridis(np.linspace(0.1, 1, 5))

fig1, ax1 = plt.subplot_mosaic([["spec_re_nonreg", "spec_re_reg"],
                                ["spec_re_nonreg", "spec_re_reg"],
                                ["spec_im_nonreg", "spec_im_reg"],
                                ["spec_im_nonreg", "spec_im_reg"],
                                ["m", "tau"],
                                ["m", "tau"],
                                ["m", "tau"]])

fig1.set_figwidth(13)
fig1.set_figheight(12)

jj = 0

for ii in [0, 3, 6, 9, 12]:
    ax1["spec_re_nonreg"].semilogx(w_nonreg[ii].real, sig_nonreg[ii].real,
                                   color=col[jj], linewidth=2)
    ax1["spec_im_nonreg"].semilogx(w_nonreg[ii].real, sig_nonreg[ii].imag,
                                   color=col[jj], linewidth=2)
    ax1["spec_re_reg"].semilogx(w_reg[ii].real, sig_reg[ii].real,
                                color=col[jj], linewidth=2)
    ax1["spec_im_reg"].semilogx(w_reg[ii].real, sig_reg[ii].imag,
                                color=col[jj], linewidth=2)
    jj = jj+1
    
ax1["m"].plot(np.array(R), m_nonreg, "-s", fillstyle="none", color=col[1],
              linewidth=2, markersize=8)
ax1["m"].plot(np.array(R), m_reg, "-o", fillstyle="none", color=col[3],
              linewidth=2, markersize=8)

sax1m = ax1["m"].secondary_xaxis("top", functions=(lambda x: x/4.75e-3,
                                                   lambda x: x/4.75e-3))
sax1m.set_xlabel("$R/a$ [-]")

sax1t = ax1["tau"].secondary_xaxis("top", functions=(lambda x: x/4.75e-3,
                                                     lambda x: x/4.75e-3))
sax1t.set_xlabel("$R/a$ [-]")


ax1["m"].set_ylim([0, 0.15])

ax1["tau"].plot(np.array(R), tau_nonreg, "-s", fillstyle="none", color=col[1],
                linewidth=2, markersize=8)
ax1["tau"].plot(np.array(R), tau_reg, "-o", fillstyle="none", color=col[3],
                linewidth=2, markersize=8)

ax1["spec_re_nonreg"].set_title("(a)", loc="left")
ax1["spec_re_reg"].set_title("(c)", loc="left")
ax1["spec_im_nonreg"].set_title("(b)", loc="left")
ax1["spec_im_reg"].set_title("(d)", loc="left")
ax1["m"].set_title("(e)", loc="left")
ax1["tau"].set_title("(f)", loc="left")

ax1["spec_re_nonreg"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_nonreg"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["spec_re_reg"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_reg"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["m"].set_ylabel("$m$ [-]")
ax1["tau"].set_ylabel("$\u03C4$ [s]")

ax1["spec_re_nonreg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_nonreg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_re_reg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_reg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["m"].set_xlabel("$R$ [m]")
ax1["tau"].set_xlabel("$R$ [m]")

ax1["spec_im_nonreg"].legend(["$R = 0.8$ cm", "$R = 1.4$ cm", "$R = 2.0$ cm",
                              "$R = 2.6$ cm", "$R = 3.2$ cm"], loc=2,
                             fontsize=12)
ax1["m"].legend(["Original spectra", "Re-scaled spectra"])

fig1.tight_layout()

fig1.savefig("Figure_03.png",dpi=300,bbox_inches="tight")

nu = np.logspace(np.log10(nu_nonreg[0]), np.log10(nu_nonreg[-1]))

fig2, ax2 = plt.subplots(figsize=(7,5))
ax2.loglog(np.array(nu_nonreg), m_nonreg, "o", fillstyle="none",
             color="black", markersize=8, linewidth=2)
ax2.loglog(np.array(nu), gurin(nu), "-",color="tab:blue", linewidth=2)
ax2.loglog(np.array(nu), revil(nu), "--", color="tab:red", linewidth=2)
ax2.set_xticks([2e-3, 4e-3, 1e-2, 2e-2, 4e-2],
               ["0.2%", "0.4%", "1.0%", "2.0%", "4.0%"])
ax2.legend(["Numerical model", "Gurin et al. (2015)", "Revil et al. (2015)"])
ax2.set_xlabel("$\u03BD_{eff}$ [-]")
ax2.set_ylabel("$m [-]$")

fig2.tight_layout()

fig2.savefig("figure_04.pdf",dpi=300,bbox_inches="tight")