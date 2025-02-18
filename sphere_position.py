import numpy as np
import matplotlib.pyplot as plt
from functions import load_data as ld
from functions import coleColeFit as ccf

plt.rcParams.update({'font.size': 13})

erg_1 = []
erg_3 = []
erg_5 = []

fname = ["00e-3", "05e-3", "10e-3", "15e-3", "20e-3", "25e-3", "30e-3",
         "35e-3", "40e-3", "45e-3"]

z = []

for name in fname:
    z.append(float(name))
    erg_1.append(ld("position_sphere/p0_1e-2", "z0_"+name, -1e-2, 1e-2, 0,
                    24e-3, 4.75e-3, 102e-3))
    erg_3.append(ld("position_sphere/p0_3e-2", "z0_"+name, -3e-2, 3e-2, 0,
                    24e-3, 4.75e-3, 102e-3))
    erg_5.append(ld("position_sphere/p0_5e-2", "z0_"+name, -5e-2, 5e-2, 0,
                    24e-3, 4.75e-3, 102e-3))
    
w_1 = []
sig_1 = []
m_1 = []
tau_1 = []
c_1 = []

for obj in erg_1:
    w_1.append(obj[0])
    sig_1.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_1.append(sol[1])
    tau_1.append(sol[2])
    c_1.append(sol[3])
    
w_3 = []
sig_3 = []
m_3 = []
tau_3 = []
c_3 = []

for obj in erg_3:
    w_3.append(obj[0])
    sig_3.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_3.append(sol[1])
    tau_3.append(sol[2])
    c_3.append(sol[3])
    
w_5 = []
sig_5 = []
m_5 = []
tau_5 = []
c_5 = []

for obj in erg_5:
    w_5.append(obj[0])
    sig_5.append(obj[1])
    sol = ccf(obj[0], obj[1])
    m_5.append(sol[1])
    tau_5.append(sol[2])
    c_5.append(sol[3])

col = plt.cm.viridis(np.linspace(0.1, 1, 5))

fig1, ax1 = plt.subplot_mosaic([["spec_re_1", "spec_re_3"],
                                ["spec_im_1", "spec_im_3"],
                                ["spec_re_5", "m"],
                                ["spec_im_5", "m"]])

fig1.set_figwidth(14)
fig1.set_figheight(12)

jj = 0

for ii in [0, 2, 4, 6, 8]:
    ax1["spec_re_1"].semilogx(w_1[ii].real, sig_1[ii].real, color=col[jj],
                              linewidth=2)
    ax1["spec_re_3"].semilogx(w_3[ii].real, sig_3[ii].real, color=col[jj],
                              linewidth=2)
    ax1["spec_re_5"].semilogx(w_5[ii].real, sig_5[ii].real, color=col[jj],
                              linewidth=2)
    ax1["spec_im_1"].semilogx(w_1[ii].real, sig_1[ii].imag, color=col[jj],
                              linewidth=2)
    ax1["spec_im_3"].semilogx(w_3[ii].real, sig_3[ii].imag, color=col[jj],
                              linewidth=2)
    ax1["spec_im_5"].semilogx(w_5[ii].real, sig_5[ii].imag, color=col[jj],
                              linewidth=2)
    jj = jj+1
    
ax1["m"].plot(np.array(z), m_1, "-o", fillstyle="none", color=col[0],
                linewidth=2, markersize=9)
ax1["m"].plot(np.array(z), m_3, "-s", fillstyle="none", color=col[2],
                linewidth=2, markersize=9)
ax1["m"].plot(np.array(z), m_5, "-d", fillstyle="none", color=col[4],
                linewidth=2, markersize=9)

ax1["spec_re_1"].set_title("(a)", loc="left")
ax1["spec_re_3"].set_title("(c)", loc="left")
ax1["spec_re_5"].set_title("(e)", loc="left")
ax1["spec_im_1"].set_title("(b)", loc="left")
ax1["spec_im_3"].set_title("(d)", loc="left")
ax1["spec_im_5"].set_title("(f)", loc="left")
ax1["m"].set_title("(g)", loc="left")

ax1["spec_re_1"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_1"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["spec_re_3"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_3"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["spec_re_5"].set_ylabel("Re($\u03C3/\u03C3_w$) [-]")
ax1["spec_im_5"].set_ylabel("Im($\u03C3/\u03C3_w$) [-]")
ax1["m"].set_ylabel("m [-]")

ax1["spec_re_1"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_1"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_re_3"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_3"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_re_5"].set_xlabel("$\u03C9$ [rad/s]")
ax1["spec_im_5"].set_xlabel("$\u03C9$ [rad/s]")
ax1["m"].set_xlabel("$\u0394 z_s$ [m]")

ax1["spec_im_1"].legend(["$\u0394 z_s = 0$ cm", "$\u0394z_s = 1$ cm",
                         "$\u0394 z_s = 2$ cm", "$\u0394 z_s = 3$ cm",
                         "$\u0394 z_s = 4$ cm"], loc=2, fontsize=12)
ax1["m"].legend(["$\u0394 z_p = 2$ cm", "$\u0394 z_p = 6$ cm",
                 "$\u0394 z_p = 10$ cm"])

fig1.tight_layout()

fig1.savefig("Figure_06.png",dpi=300,bbox_inches="tight")