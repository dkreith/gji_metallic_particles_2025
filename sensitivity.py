import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import uniform_filter

plt.rcParams.update({'font.size': 15})

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = "seismic"

dat = np.loadtxt("sensitivity.txt", skiprows=5)
x_dat = dat[:,0]
y_dat = dat[:,1]
r_dat = dat[:,2]

x_val = sorted(set(x_dat))
y_val = sorted(set(y_dat))

s_dat = ((r_dat)-(1))

s_mat = np.zeros((len(y_val),len(x_val)))

for ii in range(len(y_val)):
    for jj in range(len(x_val)):
        s_mat[ii,jj] = s_dat[ii+jj*len(y_val)]

maxi = max(s_dat)*1.1/1.1*0.2

sp = np.zeros(s_mat.shape)
sn = np.zeros(s_mat.shape)

for ii in range(s_mat.shape[0]):
    for jj in range(s_mat.shape[1]):
        if s_mat[ii,jj] >= 0:
            sp[ii,jj] = s_mat[ii,jj]
            sn[ii,jj] = np.nan
        else:
            sp[ii,jj] = np.nan
            sn[ii,jj] = s_mat[ii,jj]

minlog = -4.75
maxlog = -3.25

vmin = minlog
vmax = maxlog
dv = vmax-vmin
n_cbar = 4
ticks_old = np.linspace(vmin, vmax, n_cbar)
ticks = np.linspace(vmin-dv, vmax, n_cbar*2-1)

ticks_n = []
ticks_p = []
for tick in ticks_old:
    tick_str = str(tick)
    ticks_n.append("$-10^{"+tick_str+"}$")
    ticks_p.append("$10^{"+tick_str+"}$")

ticks_n.reverse()
ticks_n.pop()
ticks_p.pop(0)
ticks_new = ticks_n+["$0$"]+ticks_p

red = truncate_colormap(plt.get_cmap(cmap), 0.5, 1.0)
blue = truncate_colormap(plt.get_cmap(cmap+"_r"), 0.5, 1.0)

levels1p = [-4+1/3, -4+2/3, -3, -3+1/3, -3+2/3, -2, -2+1/3, -2+2/3, -1]
levels2p = [-6+2/3, -5, -5+1/3, -5+2/3, -4,]
levels1n = [-4, -4+1/3, -4+2/3, -3, -3+1/3, -3+2/3, -2, -2+1/3, -2+2/3, -1]
levels2n = [-6+2/3, -5, -5+1/3, -5+2/3]

s_mat_f = uniform_filter(s_mat, size=10)

sp_f = np.zeros(s_mat.shape)
sn_f = np.zeros(s_mat.shape)

for ii in range(s_mat_f.shape[0]):
    for jj in range(s_mat_f.shape[1]):
        if s_mat_f[ii,jj] >= 0:
            sp_f[ii,jj] = s_mat_f[ii,jj]
            sn_f[ii,jj] = np.nan
        else:
            sp_f[ii,jj] = np.nan
            sn_f[ii,jj] = s_mat_f[ii,jj]

fig, ax = plt.subplots(figsize=(13,5))
cpp = ax.pcolormesh(np.array(x_val)-51e-3, y_val, np.log10(sp), 
                    vmin=minlog, vmax=maxlog, cmap=red)
cpn = ax.pcolormesh(np.array(x_val)-51e-3, y_val, np.log10(np.abs(sn)),
                    vmin=minlog, vmax=maxlog, cmap=blue)
cp_dummy = ax.pcolormesh(np.array(x_val)-51e-3, y_val, s_mat*np.nan,
                         vmin=vmin-dv, vmax=vmax, cmap=cmap)

cb_dummy = fig.colorbar(cp_dummy, ticks=ticks)
cb_dummy.ax.set_yticklabels(ticks_new)

ax.contour(np.array(x_val)-51e-3, y_val, np.log10(sp),
           colors="black", linestyles="dashed", levels=levels1p)
ax.contour(np.array(x_val)-51e-3, y_val, np.log10(np.abs(sn)),
           colors="black", linestyles="dashed", levels=levels1n)
ax.contour(np.array(x_val)-51e-3, y_val, np.log10(sp_f),
           colors="black", linestyles="dashed", levels=levels2p)
ax.contour(np.array(x_val)-51e-3, y_val, np.log10(np.abs(sn_f)),
           colors="black", linestyles="dashed", levels=levels2n)

ax.set_xlim([-51e-3, 51e-3])
ax.set_ylim([-2.4e-2, 2.4e-2])
ax.set_xlabel("$z$ [m]")
ax.set_ylabel("$r$ [m]")
cb_dummy.set_label("$S$ [-]")
fig.tight_layout()
fig.savefig("figure_02.png",dpi=300,bbox_inches="tight")