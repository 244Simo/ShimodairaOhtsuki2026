### import library
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from matplotlib.colors import ListedColormap

### import data
data10 = np.loadtxt("data_q10.csv", delimiter=",")
data095 = np.loadtxt("data_q095.csv", delimiter=",")
data09 = np.loadtxt("data_q09.csv", delimiter=",")
data08 = np.loadtxt("data_q08.csv", delimiter=",")

### make colormap
cmap = ListedColormap([
    "white",  # 0
    "gray",  # 1
    "lightgreen",  # 2
    "paleturquoise",  # 3
    "darkcyan"   # 4
])

### prepare for ticks
beta_ticks_val = np.arange(0, 31, 5)   # 0,5,10,...30
gamma_ticks_val = np.arange(0, 31, 5)
beta_ticks_idx = beta_ticks_val / 0.1
gamma_ticks_idx = gamma_ticks_val / 0.1

### convert each data(q=1.0)
beta_idx_10 = data10[:,0].astype(int)
gamma_idx_10 = data10[:,1].astype(int)
values_10 = data10[:,8]
# obtain size
beta_max_10 = beta_idx_10.max() + 1
gamma_max_10 = gamma_idx_10.max() + 1
# convert to matrix
Z_10 = np.zeros((beta_max_10, gamma_max_10))
for b, g, v in zip(beta_idx_10, gamma_idx_10, values_10):
    Z_10[b, g] = v

### convert each data(q=0.95)
beta_idx_095 = data095[:,0].astype(int)
gamma_idx_095 = data095[:,1].astype(int)
values_095 = data095[:,8]
# obtain size
beta_max_095 = beta_idx_095.max() + 1
gamma_max_095 = gamma_idx_095.max() + 1
# convert to matrix
Z_095 = np.zeros((beta_max_095, gamma_max_095))
for b, g, v in zip(beta_idx_095, gamma_idx_095, values_095):
    Z_095[b, g] = v

### convert each data(q=0.9)
beta_idx_09 = data09[:,0].astype(int)
gamma_idx_09 = data09[:,1].astype(int)
values_09 = data09[:,8]
# obtain size
beta_max_09 = beta_idx_09.max() + 1
gamma_max_09 = gamma_idx_09.max() + 1
# convert to matrix
Z_09 = np.zeros((beta_max_09, gamma_max_09))
for b, g, v in zip(beta_idx_09, gamma_idx_09, values_09):
    Z_09[b, g] = v

### convert each data(q=0.8)
beta_idx_08 = data08[:,0].astype(int)
gamma_idx_08 = data08[:,1].astype(int)
values_08 = data08[:,8]
# obtain size
beta_max_08 = beta_idx_08.max() + 1
gamma_max_08 = gamma_idx_08.max() + 1
# convert to matrix
Z_08 = np.zeros((beta_max_08, gamma_max_08))
for b, g, v in zip(beta_idx_08, gamma_idx_08, values_08):
    Z_08[b, g] = v

### setting object for plot
fig = plt.figure(figsize=(8,8))
fig.subplots_adjust(hspace=0.4, wspace=0.4)

## q=1.0
ax1 = fig.add_subplot(2,2,1)
ax1.set_aspect('equal')
# caption
ax1.text(
    -0.3, 1.12, "(a)$\it{q}$=1",
    transform=ax1.transAxes,
    fontsize=14, va='top', ha='left'
)

ax1.imshow(Z_10, origin="lower", cmap=cmap)
ax1.set_yticks(beta_ticks_idx, beta_ticks_val)
ax1.set_xticks(gamma_ticks_idx, gamma_ticks_val)
ax1.set_xlabel('Effect of teaching, $γ$', font='arial')
ax1.set_ylabel('Efficiency of observation learning, $β$', font='arial')

## q=0.95
ax2 = fig.add_subplot(2,2,2)
ax2.set_aspect('equal')
# caption
ax2.text(
    -0.3, 1.12, "(b)$\it{q}$=0.95",
    transform=ax2.transAxes,
    fontsize=14, va='top', ha='left'
)
ax2.imshow(Z_095, origin="lower", cmap=cmap)
ax2.set_yticks(beta_ticks_idx, beta_ticks_val)
ax2.set_xticks(gamma_ticks_idx, gamma_ticks_val)
ax2.set_xlabel('Effect of teaching, $γ$', font='arial')
ax2.set_ylabel('Efficiency of observation learning, $β$', font='arial')

## q=0.9
ax3 = fig.add_subplot(2,2,3)
ax3.set_aspect('equal')
# caption
ax3.text(
    -0.3, 1.12, "(c)$\it{q}$=0.9",
    transform=ax3.transAxes,
    fontsize=14, va='top', ha='left'
)
ax3.imshow(Z_09, origin="lower", cmap=cmap)
ax3.set_yticks(beta_ticks_idx, beta_ticks_val)
ax3.set_xticks(gamma_ticks_idx, gamma_ticks_val)
ax3.set_xlabel('Effect of teaching, $γ$', font='arial')
ax3.set_ylabel('Efficiency of observation learning, $β$', font='arial')

## q=0.8
ax4 = fig.add_subplot(2,2,4)
ax4.set_aspect('equal')
# caption
ax4.text(
    -0.3, 1.12, "(d)$\it{q}$=0.8",
    transform=ax4.transAxes,
    fontsize=14, va='top', ha='left'
)
ax4.imshow(Z_08, origin="lower", cmap=cmap)
ax4.set_yticks(beta_ticks_idx, beta_ticks_val)
ax4.set_xticks(gamma_ticks_idx, gamma_ticks_val)
ax4.set_xlabel('Effect of teaching, $γ$', font='arial')
ax4.set_ylabel('Efficiency of observation learning, $β$', font='arial')

### legend(bottom of figures)
## handle for legend
colors = ["gray", "paleturquoise", "white", "darkcyan", "lightgreen"]
labels = ["Only asocial ESS", "Asocial and teaching ESS coexist", "Only non-teaching ESS", "Non-teaching and teaching ESS coexist", "Only teaching ESS"]

handles = [
    Patch(facecolor=c, edgecolor="black", label=l)
    for c, l in zip(colors, labels)
]

fig.legend(
    handles=handles,
    loc="lower center",
    ncol=3,
    frameon=False,
    fontsize=10,
    bbox_to_anchor=(0.50, 0)
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

## visualize
plt.show()

## save figure
fig.savefig('Fig6.eps', bbox_inches="tight")
plt.close(fig)