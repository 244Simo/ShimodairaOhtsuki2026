### import library
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### import data(csv file)
df = pd.read_csv(
    ## change the name of data file accordingly
    "timelapse_data_I1.csv",
    #"timelapse_data_I2.csv",
    #"timelapse_data_I3.csv",
    header=None,
    names=["step","SL","IL","EX","TE","Z","W"]
)
## pick up data which will be used in plotting(set step size for plotting)
df_plot = df[df["step"] % 600 == 0] #FigureI1
#df_plot = df[df["step"] % 50 == 0]  #FigureI2
#df_plot = df[df["step"] % 100 == 0] #FigureI3

## set parameters for plot
#Figure I1
x_max=60000
z_max=2.5
w_max=1.06

#Figure I2
#x_max=3000
#z_max=5
#w_max=4

#Figure I3
#x_max=10000
#z_max=5
#w_max=14

### plot
fig, ax = plt.subplots(3,2, figsize=(6.7,9.4))
# social learning
ax[0,0].text(
    -0.25, 1.12, "(a)$\\tau_{\\mathrm{SL}}$",
    transform=ax[0,0].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[0,0].scatter(df_plot.step, df_plot.SL, s=1, color="blue",marker=",", rasterized=False)
ax[0,0].set_xlabel("generations")
ax[0,0].set_ylabel("lifetime allocation")
ax[0,0].set_ylim(0,1)
ax[0,0].set_xlim(0,x_max)

# individual learning
ax[0,1].text(
    -0.25, 1.12, "(b)$\\tau_{\\mathrm{IL}}$",
    transform=ax[0,1].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[0,1].scatter(df_plot.step, df_plot.IL, s=0.64, color="blue",marker=",", rasterized=False)
ax[0,1].set_xlabel("generations")
ax[0,1].set_ylabel("lifetime allocation")
ax[0,1].set_ylim(0,1)
ax[0,1].set_xlim(0,x_max)

# exploitation
ax[1,0].text(
    -0.25, 1.12, "(c)$\\tau_{\\mathrm{EX}}$",
    transform=ax[1,0].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[1,0].scatter(df_plot.step, df_plot.EX, s=0.64, color="blue", marker=",",rasterized=False)
ax[1,0].set_xlabel("generations")
ax[1,0].set_ylabel("lifetime allocation")
ax[1,0].set_ylim(0,1)
ax[1,0].set_xlim(0,x_max)

# teaching
ax[1,1].text(
    -0.25, 1.12, "(d)$\\tau_{\\mathrm{TE}}$",
    transform=ax[1,1].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[1,1].scatter(df_plot.step, df_plot.TE, s=0.25, color="blue",marker=",", rasterized=False)
ax[1,1].set_xlabel("generations")
ax[1,1].set_ylabel("lifetime allocation")
ax[1,1].set_ylim(0,1)
ax[1,1].set_xlim(0,x_max)

# skill level
ax[2,0].text(
    -0.25, 1.12, "(e)$z$",
    transform=ax[2,0].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[2,0].scatter(df_plot.step, df_plot.Z, s=0.64, color="blue", marker=",",rasterized=False)
ax[2,0].set_xlabel("generations")
ax[2,0].set_ylabel("cultural level")
ax[2,0].set_ylim(0,z_max)
ax[2,0].set_xlim(0,x_max)

# fitness
ax[2,1].text(
    -0.25, 1.12, "(f)$w$",
    transform=ax[2,1].transAxes,
    fontsize=14, va='top', ha='left'
)
ax[2,1].scatter(df_plot.step, df_plot.W, s=0.64, color="blue", marker=",",rasterized=False)
ax[2,1].set_xlabel("generations")
ax[2,1].set_ylabel("fitness")
ax[2,1].set_ylim(1,w_max)
ax[2,1].set_xlim(0,x_max)

plt.subplots_adjust(
    left=0.1,
    right=0.95,
    bottom=0.08,
    top=0.95,
    wspace=0.3,
    hspace=0.35
)

# visualize
#plt.show()

# save figure（eps style）
# change file name accordingly
plt.savefig("SupFigI.eps")