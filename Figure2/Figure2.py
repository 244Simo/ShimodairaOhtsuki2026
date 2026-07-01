### import library
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

### define function
## boundary(only teaching COS)
# upper
def f_up(gamma):
    return 0.5*(gamma+np.sqrt(gamma**2-8*gamma))
# lower(8<γ<9)
def f_sub(gamma):
    return 0.5*(gamma-np.sqrt(gamma**2-8*gamma))

# fitness on non-teaching COS（COS in Kobayashi et al 2015）
def f_w_noteaching(beta):
    return (1/(beta**2))*np.exp(beta-2)
# fitness on teaching COS
def f_w_teaching(beta,gamma):
    return (1/gamma)*(-1+((beta+gamma)**2/(4*gamma))*(0.5-np.sqrt(0.25-(4*gamma)/(beta+gamma)**2)))*np.exp(((beta+gamma)**2/(4*gamma))*(0.5-(4*gamma)/((beta+gamma)**2)+np.sqrt(0.25-(4*gamma)/(beta+gamma)**2)))

# to set the initial value of β(for COS comparison)
def f_beta_border(gamma):
    return 4*np.sqrt(gamma)-gamma

### settings for the plot
## horizontal axis(γ)
gamma_min=0
gamma_max=20
## vertical axis(β)
beta_min=0
beta_max=20
## stepsize for numerical calculation
delta=10**(-2)
## array for plot
gamma_array=[]
beta_upper_array=[]
beta_sub_array=[]

#### main program
## calculation of the boundary point between teaching COS and non-teaching COS on β=2 
## We can start the calculation from γ=12(γ=12 is not the boundary but it satisfies (β+γ)**2-16γ>0)
gamma_intersection_temp=12
while f_w_noteaching(2)>f_w_teaching(2,gamma_intersection_temp):
    gamma_intersection_temp += delta
gamma_intersection=gamma_intersection_temp

## calculate each boundary
for i in np.arange(8, gamma_max+delta/2, delta):
    gamma=np.round(i, 5)    # round the γ-value
    if gamma <= 9:
        gamma_array.append(gamma)
        beta_upper_array.append(f_up(gamma))
        beta_sub_array.append(f_sub(gamma))
    # if γ>9, compare teaching COS and non-teaching COS/asocial COS  
    else:
        beta_border_initial=f_beta_border(gamma) #set the initial value of β(obtained from numerical analysis)
        beta_border_temp=beta_border_initial
        # if β>2, compare teaching COS and non-teaching COS
        if gamma < gamma_intersection:
            if f_w_noteaching(beta_border_temp)>f_w_teaching(beta_border_temp,gamma):
                while f_w_noteaching(beta_border_temp)>f_w_teaching(beta_border_temp,gamma):
                    beta_border_temp += delta
                beta_border = beta_border_temp
        
                gamma_array.append(gamma)
                beta_upper_array.append(f_up(gamma))
                beta_sub_array.append(beta_border)
                beta_border=0
        # if β<2, compare teaching COS and asocial COS
        else:
            if 0.25>f_w_teaching(beta_border_temp,gamma):
                while 0.25>f_w_teaching(beta_border_temp,gamma):
                    beta_border_temp += delta
                beta_border = beta_border_temp
        
                gamma_array.append(gamma)
                beta_upper_array.append(f_up(gamma))
                beta_sub_array.append(beta_border)
                beta_border=0

### set the object for plot
fig = plt.figure()
ax = fig.add_subplot()

# fillin each region（asocial COS→teaching COS）
ax.axhspan(0, 2, color='gray') # asocial COS
ax.fill_between(gamma_array, beta_sub_array, beta_upper_array, color='lightgreen') # teaching COS
# set plot square
ax.set_aspect('equal')

## setting for axes
# range
ax.set_xlim(gamma_min, gamma_max)
ax.set_ylim(beta_min, beta_max)
# step for ticks
step=2 # value of step
ax.set_xticks(np.arange(gamma_min, gamma_max+1/step, step))
ax.set_yticks(np.arange(beta_min, beta_max+1/step, step))
# label of axes
ax.set_xlabel('Effect of teaching, $γ$', font='arial')
ax.set_ylabel('Efficiency of observation learning, $β$', font='arial')

# add text in figure(if plot within figure)
#ax.text(0.05,0.02,"Asocial COS",color="white",fontsize=16,transform=ax.transAxes)
#ax.text(0.05,0.9,"Non-teaching COS",color="black",fontsize=16,transform=ax.transAxes)
#ax.text(0.53,0.23,"Teaching COS",color="green",fontsize=16,transform=ax.transAxes)

# legend(outside figure)
handles = [
    Patch(
        facecolor='white',
        edgecolor='0.2',
        linewidth=1.5,
        label='Non-teaching COS'
    ),
    Patch(
        facecolor='lightgreen',
        edgecolor='none',
        label='Teaching COS'
    ),
    Patch(
        facecolor='gray',
        edgecolor='none',
        label='Asocial COS'
    ),
]
ax.legend(
    handles=handles,
    frameon=False,
    loc='upper left'
)

# visualize
plt.show()

## save figure（eps style）
fig.savefig('Fig2.eps')
plt.close(fig)




