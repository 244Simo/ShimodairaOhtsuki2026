### import library
import matplotlib.pyplot as plt
import numpy as np

### define fuction
## the top of q(t)
def q_top(t,m,q):
    return -(1+q)*(t**2)+np.sqrt(q)*m*t-1-np.exp(t**2)

## the bottom of q(t)
def q_bottom(t):
    return 1-np.exp(t**2)

## function of q(t)（top/bottom）
def f_q(t,m,q):
    return q_top(t,m,q)/q_bottom(t)

### paerameter setting for the plot
## horizontal axis
x_min=0.1
x_max=3
# array for the plot
x = np.arange(x_min, x_max, 0.01) 

## vertical axis
y_min=0
y_max=1.5

### set the object for plot
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(x_min,x_max)
ax1.set_ylim(y_min,y_max)
ax1.set_xlabel(r"$\xi$")
ax1.set_ylabel(r"$f_{q}(\xi)$")

## plot which depends on parameters
## q=1
ax1.plot(x,f_q(x,3.9,1),color='red',label='m=3.9,q=1') # m=3.9,q=1
ax1.plot(x,f_q(x,4,1),color='black',label='m=4,q=1') # m=4,q=1
ax1.plot(x,f_q(x,5,1),color='blue',label='m=5,q=1') # m=5,q=1
ax1.plot(x,f_q(x,6,1),color='green',label='m=6,q=1') # m=6,q=1
## q=0.9
ax1.plot(x,f_q(x,3.9,0.9),color='red',linestyle="dashdot",label='m=3.9,q=0.9') # m=3.9,q=0.9
ax1.plot(x,f_q(x,4,0.9),color='black',linestyle="dashdot",label='m=4,q=0.9') # m=4,q=0.9
ax1.plot(x,f_q(x,5,0.9),color='blue',linestyle="dashdot",label='m=5,q=0.9') # m=5,q=0.9
ax1.plot(x,f_q(x,6,0.9),color='green',linestyle="dashdot",label='m=6,q=0.9') # m=6,q=0.9

## plot the auxiliary line
ax1.axhline(1, linestyle='dotted', color ='gray') # y=1
ax1.axhline(0.8, linestyle='dotted', color ='gray') # y=0.8
ax1.axhline(0.6, linestyle='dotted', color ='gray') # y=0.6
ax1.axhline(0.4, linestyle='dotted', color ='gray') # y=0.4
ax1.axhline(0.2, linestyle='dotted', color ='gray') # y=0.2
ax1.axvline(1, linestyle='dotted', color ='lightgrey', alpha=0.5) # ξ=1

## plot the legend
ax1.legend()

## output the figure
plt.show()
#fig.savefig('SupFigC2.eps')
#plt.close(fig)