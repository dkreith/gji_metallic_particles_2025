import numpy as np
import matplotlib.pyplot as plt
from functions import load_data as ld
from functions import calc_nu
from functions import wong as wg
from functions import findmax

w = np.logspace(-7, 5, 1200)
nu = calc_nu(-3e-2, 3e-2, 0, 24e-3, 4.75e-3, 102e-3)

erg_1 = []
erg_2 = []
erg_3 = []
erg_4 = []

fname = ["0", "100e-6", "178e-6", "316e-6", "562e-6",
         "100e-5", "178e-5", "316e-5", "562e-5",
         "100e-4", "178e-4", "316e-4", "562e-4",
         "100e-3", "178e-3", "316e-3"]

for name in fname:
    erg_1.append(ld("reaction_current_valence/z_1", "c3_"+name, -3e-2, 3e-2, 0,
                    24e-3, 4.75e-3, 102e-3, ind=2))
    erg_2.append(ld("reaction_current_valence/z_2", "c3_"+name, -3e-2, 3e-2, 0,
                    24e-3, 4.75e-3, 102e-3, ind=2))
    erg_3.append(ld("reaction_current_valence/z_3", "c3_"+name, -3e-2, 3e-2, 0,
                    24e-3, 4.75e-3, 102e-3, ind=2))
    erg_4.append(ld("reaction_current_valence/z_4", "c3_"+name, -3e-2, 3e-2, 0,
                    24e-3, 4.75e-3, 102e-3, ind=2))

idx_list = [0, 3, 6, 9, 12]
freq_list = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34,
             36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66,
             68, 70, 72, 74, 76, 78, 80, 82, 84]

col = plt.cm.viridis(np.linspace(0.1, 1, 5))

w_1 = []
sig_1 = []
sigmax_1 = []

for obj in erg_1:
    w_1.append(obj[0])
    sig_1.append(obj[1])
    sigmax_temp, tau_temp = findmax(obj[0], obj[1].imag)
    sigmax_1.append(sigmax_temp)

for ii in range(len(fname)):
    if float(fname[ii]) > 1:
        sigmax_1[ii] = np.nan
    
w_1_plot = []
sig_1_plot = []
sig_wg_1 = []

for ii in idx_list:
    w_temp = []
    sig_temp = []
    for jj in freq_list:
        w_temp.append(w_1[ii][jj+7])
        sig_temp.append(sig_1[ii][jj+7])
    w_1_plot.append(w_temp)
    sig_1_plot.append(sig_temp)
    c3 = float(fname[ii])
    if c3 < 1:
        sig_wg_1.append(wg(w, 1, c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2, nu))
    else:
        sig_wg_1.append(wg(w, 1, c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu)*np.nan)

w_2 = []
sig_2 = []
sigmax_2 = []

for obj in erg_2:
    w_2.append(obj[0])
    sig_2.append(obj[1])
    sigmax_temp, tau_temp = findmax(obj[0], obj[1].imag)
    sigmax_2.append(sigmax_temp)

for ii in range(len(fname)):
    if 4*float(fname[ii]) > 1:
        sigmax_2[ii] = np.nan
    
w_2_plot = []
sig_2_plot = []
sig_wg_2 = []

for ii in idx_list:
    w_temp = []
    sig_temp = []
    for jj in freq_list:
        w_temp.append(w_2[ii][jj+7])
        sig_temp.append(sig_2[ii][jj+7])
    w_2_plot.append(w_temp)
    sig_2_plot.append(sig_temp)
    c3 = float(fname[ii])
    if 4*c3 < 1:
        sig_wg_2.append(wg(w, 1, 4*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu))
    else:
        sig_wg_2.append(wg(w, 1, 4*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu)*np.nan)

w_3 = []
sig_3 = []
sigmax_3 = []

for obj in erg_3:
    w_3.append(obj[0])
    sig_3.append(obj[1])
    sigmax_temp, tau_temp = findmax(obj[0], obj[1].imag)
    sigmax_3.append(sigmax_temp)

for ii in range(len(fname)):
    if 9*float(fname[ii]) > 1:
        sigmax_3[ii] = np.nan
    
w_3_plot = []
sig_3_plot = []
sig_wg_3 = []

for ii in idx_list:
    w_temp = []
    sig_temp = []
    for jj in freq_list:
        w_temp.append(w_3[ii][jj+7])
        sig_temp.append(sig_3[ii][jj+7])
    w_3_plot.append(w_temp)
    sig_3_plot.append(sig_temp)
    c3 = float(fname[ii])
    if 9*c3 < 1:
        sig_wg_3.append(wg(w, 1, 9*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu))
    else:
        sig_wg_3.append(wg(w, 1, 9*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu)*np.nan)
        
w_4 = []
sig_4 = []
sigmax_4 = []

for obj in erg_4:
    w_4.append(obj[0])
    sig_4.append(obj[1])
    sigmax_temp, tau_temp = findmax(obj[0], obj[1].imag)
    sigmax_4.append(sigmax_temp)
    
for ii in range(len(fname)):
    if 16*float(fname[ii]) > 1:
        sigmax_4[ii] = np.nan
    
w_4_plot = []
sig_4_plot = []
sig_wg_4 = []

for ii in idx_list:
    w_temp = []
    sig_temp = []
    for jj in freq_list:
        w_temp.append(w_4[ii][jj+7])
        sig_temp.append(sig_4[ii][jj+7])
    w_4_plot.append(w_temp)
    sig_4_plot.append(sig_temp)
    c3 = float(fname[ii])
    if 16*c3 < 1:
        sig_wg_4.append(wg(w, 1, 16*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu))
    else:
        sig_wg_4.append(wg(w, 1, 16*c3, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2,
                           nu)*np.nan)

sigmaxwg_1 = []
sigmaxwg_2 = []
sigmaxwg_3 = []
sigmaxwg_4 = []

c3 = np.logspace(-4, np.log10(float(fname[-1])), 100)
for c3_temp in c3:
    
    if c3_temp < 1:
        sig_temp_1 = wg(w, 1, c3_temp, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2, nu)
        sigmax_temp_1, tau_temp = findmax(w, sig_temp_1.imag.squeeze())
        sigmaxwg_1.append(sigmax_temp_1)
    else:
        sigmaxwg_1.append(np.nan)
    
    if 4*c3_temp < 1:
        sig_temp_2 = wg(w, 1, 4*c3_temp, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2/2,
                        nu)
        sigmax_temp_2, tau_temp = findmax(w, sig_temp_2.imag.squeeze())
        sigmaxwg_2.append(sigmax_temp_2)
    else:
        sigmaxwg_2.append(np.nan)
    
    if 9*c3_temp < 1:
        sig_temp_3 = wg(w, 1, 9*c3_temp, 293, 80, 5e-8, 4.75e-3, 1e-10, 1e-2/3,
                        nu)
        sigmax_temp_3, tau_temp = findmax(w, sig_temp_3.imag.squeeze())
        sigmaxwg_3.append(sigmax_temp_3)
    else:
        sigmaxwg_3.append(np.nan)
    
    if 16*c3_temp < 1:
        sig_temp_4 = wg(w, 1, 16*c3_temp, 293, 80, 5e-8, 4.75e-3, 1e-10,
                        1e-2/4, nu)
        sigmax_temp_4, tau_temp = findmax(w, sig_temp_4.imag.squeeze())
        sigmaxwg_4.append(sigmax_temp_4)
    else:
        sigmaxwg_4.append(np.nan)

plt.rcParams.update({'font.size': 13})

fig, ax = plt.subplots(2,2,figsize=(12,9))

ax_scttr = []
ax_lines = []

col_idx = 0
for ii in range(len(idx_list)):
    
    scttr_temp, = ax[0,0].semilogx(np.real(w_2_plot[ii]),
                                   np.imag(sig_2_plot[ii]), "x",
                                   color=col[col_idx])
    ax_scttr.append(scttr_temp)
    lines_temp, = ax[0,0].semilogx(w, sig_wg_2[ii].imag, color=col[col_idx])
    ax_lines.append(lines_temp)
    
    ax[0,1].semilogx(np.real(w_3_plot[ii]), np.imag(sig_3_plot[ii]), "x",
                     color=col[col_idx])
    ax[0,1].semilogx(w, sig_wg_3[ii].imag, color=col[col_idx])
    
    ax[1,0].semilogx(np.real(w_4_plot[ii]), np.imag(sig_4_plot[ii]), "x",
                     color=col[col_idx])
    ax[1,0].semilogx(w, sig_wg_4[ii].imag, color=col[col_idx])
    
    col_idx = col_idx + 1

ax[0,0].legend([(ax_scttr[0], ax_lines[0]), (ax_scttr[1], ax_lines[1]),
                (ax_scttr[2], ax_lines[2]), (ax_scttr[3], ax_lines[3]),
                (ax_scttr[4], ax_lines[4])],
               ["$c_3$ = 0.00 mmol/m³", "$c_3$ = 0.32 mmol/m³",
                "$c_3$ = 1.78 mmol/m³", "$c_3$ = 10.00 mmol/m³",
                "$c_3$ = 56.20 mmol/m³"])

    
num1, = ax[1,1].semilogx([float(c3_temp) for c3_temp in fname], sigmax_1, "^", 
                         color="tab:blue")
num2, = ax[1,1].semilogx([float(c3_temp) for c3_temp in fname], sigmax_2, "o",
                         color="tab:orange")
num3, = ax[1,1].semilogx([float(c3_temp) for c3_temp in fname], sigmax_3, "s",
                         color="tab:green")
num4, = ax[1,1].semilogx([float(c3_temp) for c3_temp in fname], sigmax_4, "D",
                         color="tab:red")

ana1, = ax[1,1].semilogx(c3, sigmaxwg_1, "-", color="tab:blue", linewidth=2)
ana2, = ax[1,1].semilogx(c3, sigmaxwg_2, "--", color="tab:orange", linewidth=2)
ana3, = ax[1,1].semilogx(c3, sigmaxwg_3, ":", color="tab:green", linewidth=2)
ana4, = ax[1,1].semilogx(c3, sigmaxwg_4, "-.", color="tab:red", linewidth=2)

ax[1,1].legend([(num1, ana1), (num2, ana2), (num3, ana3), (num4, ana4)],
               ["$v_3$ = 1", "$v_3$ = 2", "$v_3$ = 3", "$v_3$ = 4"])

ax[0,0].set_title("(a)", loc="left")
ax[0,1].set_title("(b)", loc="left")
ax[1,0].set_title("(c)", loc="left")
ax[1,1].set_title("(d)", loc="left")

ax[0,0].set_xlabel("$\u03C9$ [rad/s]")
ax[0,0].set_ylabel("$\u03C3''/\u03C3_w$ [-]")
ax[0,1].set_xlabel("$\u03C9$ [rad/s]")
ax[0,1].set_ylabel("$\u03C3''/\u03C3_w$ [-]")
ax[1,0].set_xlabel("$\u03C9$ [rad/s]")
ax[1,0].set_ylabel("$\u03C3''/\u03C3_w$ [-]")
ax[1,1].set_xlabel("$c_3 / c_0$ [-]")
ax[1,1].set_ylabel("$\u03C3_{max}''/\u03C3_w$ [-]")
fig.tight_layout()

fig.savefig("figure_09.pdf", dpi=300)