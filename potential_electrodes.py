import numpy as np
import matplotlib.pyplot as plt
from functions import load_data as ld
from functions import coleColeFit as ccf

plt.rcParams.update({'font.size': 13})

erg_nonreg = []
erg_reg = []

fname = ["025e-4", "050e-4", "075e-4", "100e-4", "125e-4", "150e-4", "175e-4",
         "200e-4", "225e-4", "250e-4", "275e-4", "300e-4", "325e-4", "350e-4",
         "375e-4", "400e-4", "425e-4", "450e-4", "475e-4", "500e-4"]

p0 = []

for name in fname:
    p0_temp = float(name)
    p0.append(p0_temp)
    erg_nonreg.append(ld("potential_electrodes", "p0_"+name, -p0_temp, p0_temp,
                         0, 24e-3, 4.75e-3, 102e-3))
    erg_reg.append(ld("potential_electrodes", "p0_"+name, -p0_temp, p0_temp,
                      0, 24e-3, 4.75e-3, 102e-3, "reg"))
    
w_nonreg = []
sig_nonreg = []
m_nonreg = []
tau_nonreg = []
c_nonreg = []

for obj in erg_nonreg:
    w_nonreg.append(obj[0])
    sig_nonreg.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_nonreg.append(sol[1])
    tau_nonreg.append(sol[2])
    c_nonreg.append(sol[3])

w_reg = []
sig_reg = []
m_reg = []
tau_reg = []
c_reg = []

for obj in erg_reg:
    w_reg.append(obj[0])
    sig_reg.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_reg.append(sol[1])
    tau_reg.append(sol[2])
    c_reg.append(sol[3])

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

for ii in [1, 5, 9, 13, 17]:
    ax1["spec_re_nonreg"].semilogx(w_nonreg[ii].real, sig_nonreg[ii].real,
                                   color=col[jj], linewidth=2)
    ax1["spec_im_nonreg"].semilogx(w_nonreg[ii].real, sig_nonreg[ii].imag,
                                   color=col[jj], linewidth=2)
    ax1["spec_re_reg"].semilogx(w_reg[ii].real, sig_reg[ii].real,
                                color=col[jj], linewidth=2)
    ax1["spec_im_reg"].semilogx(w_reg[ii].real, sig_reg[ii].imag,
                                color=col[jj], linewidth=2)
    jj = jj+1
    
ax1["m"].plot(2*np.array(p0), m_nonreg, "-s", fillstyle="none", color=col[1],
                     linewidth=2, markersize=8)
ax1["m"].plot(2*np.array(p0), m_reg, "-o", fillstyle="none", color=col[3],
                  linewidth=2, markersize=8)

ax1["m"].plot([4.75e-3*2, 4.75e-3*2],[0.005, 0.045], "--k", alpha=0.5)

ax1["tau"].plot(2*np.array(p0), tau_nonreg, "-s", fillstyle="none",
                color=col[1], linewidth=2, markersize=8)
ax1["tau"].plot(2*np.array(p0), tau_reg, "-o", fillstyle="none",
                color=col[3], linewidth=2, markersize=8)

ax1["spec_re_nonreg"].set_title("(a)", loc="left")
ax1["spec_im_nonreg"].set_title("(b)", loc="left")
ax1["m"].set_title("(e)", loc="left")
ax1["spec_re_reg"].set_title("(c)", loc="left")
ax1["spec_im_reg"].set_title("(d)", loc="left")
ax1["tau"].set_title("(f)", loc="left")

ax1["spec_re_nonreg"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_nonreg"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["spec_re_reg"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_reg"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["m"].set_ylabel("m [-]")
ax1["tau"].set_ylabel("$\u03C4$ [s]")

ax1["spec_re_nonreg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_nonreg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_re_reg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_reg"].set_xlabel("$\u03C9$ [rad/s]")
ax1["m"].set_xlabel("$\u0394 z_p$ [m]")
ax1["tau"].set_xlabel("$\u0394 z_p$ [m]")

ax1["spec_im_nonreg"].legend(["$\u0394 z_p = 1$ cm", "$\u0394 z_p = 3$ cm",
                              "$\u0394 z_p = 5$ cm", "$\u0394 z_p = 7$ cm",
                              "$\u0394 z_p = 9$ cm"], loc=2, fontsize=12)

ax1["m"].legend(["Original spectra", "Re-scaled spectra"])

fig1.tight_layout()

fig1.savefig("Figure_05.png",dpi=300,bbox_inches="tight")