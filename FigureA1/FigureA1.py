### import library
import matplotlib.pyplot as plt
import numpy as np

### function
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

### array for plot
## panel(a)
gamma_array_a=[]         #γ for teaching ESS
beta_upper_array_a=[]    #β for teaching ESS(upper)
beta_sub_array_a=[]      #β for teaching ESS(lower)
## panel(b)
gamma_array_b=[]         #γ for only teaching ESS
beta_upper_array_b=[]    #β for teaching ESS(upper)
beta_sub_array_b=[]      #β for teaching ESS(lower)
## panel(c)
gamma_array_c=[]         #γ for comparison region
beta_upper_array_c=[]    #β for comparison region(upper)
beta_sub_array_c=[]      #β for comparison region(lower)

### parameter setting
## horizontal axis(γ)
gamma_min=0
gamma_max=20
## vertical axis(β)
beta_min=0
beta_max=20
## stepsize for numerical calculation
delta=10**(-2)

#### main program
## calculation of the boundary point between teaching COS and non-teaching COS on β=2 
## We can start the calculation from γ=12(γ=12 is not the boundary but it satisfies (β+γ)**2-16γ>0)
gamma_intersection_temp=12
while f_w_noteaching(2)>f_w_teaching(2,gamma_intersection_temp):
    gamma_intersection_temp += delta
gamma_intersection=gamma_intersection_temp

## calculate each boundary
for i in np.arange(8, gamma_max+delta/2, delta):
    gamma=np.round(i, 5)    # ound the γ-value
    gamma_array_b.append(gamma)
    beta_upper_array_b.append(f_up(gamma))
    beta_sub_array_b.append(f_sub(gamma))
    if gamma <= 9:
        gamma_array_a.append(gamma)
        beta_upper_array_a.append(f_up(gamma))
        beta_sub_array_a.append(f_sub(gamma))
    else: 
        gamma_array_a.append(gamma)
        beta_upper_array_a.append(f_up(gamma))
        beta_sub_array_a.append(f_beta_border(gamma))

        gamma_array_c.append(gamma)
        beta_upper_array_c.append(f_sub(gamma))
        beta_sub_array_c.append(f_beta_border(gamma))

        
### set the object
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(left=0.15, right=0.95, bottom=0.1, top=0.95, hspace=0.4, wspace=0.4)
ax1 = fig.add_subplot(2,2,1)
# caption
ax1.text(
    -0.3, 1.12, "(a)",
    transform=ax1.transAxes,
    fontsize=18, va='top', ha='left'
)

# plot lines
ax1.plot(gamma_array_a, beta_upper_array_a, color ='red', linewidth=1)
ax1.plot(gamma_array_a, beta_sub_array_a, color ='red',linewidth=1)
# fillin region
ax1.fill_between(gamma_array_a, beta_sub_array_a, beta_upper_array_a, color='hotpink')
# set plot square
ax1.set_aspect('equal')

## setting for axes
# range
ax1.set_xlim(gamma_min, gamma_max)
ax1.set_ylim(beta_min, beta_max)
# step for ticks
step=2 # stepsize
ax1.set_xticks(np.arange(gamma_min, gamma_max+1/step, step))
ax1.set_yticks(np.arange(beta_min, beta_max+1/step, step))
# label
ax1.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=16)
ax1.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=16)

## set the object
ax2 = fig.add_subplot(2,2,2)
# caption
ax2.text(
    -0.3, 1.12, "(b)",
    transform=ax2.transAxes,
    fontsize=18, va='top', ha='left'
)
# plot lines
ax2.plot(gamma_array_b, beta_upper_array_b, color ='gold',linewidth=1)
ax2.plot(gamma_array_b, beta_sub_array_b, color ='gold',linewidth=1)
ax2.axhline(2, linestyle='dotted', color='gray')
# fillin region
ax2.fill_between(gamma_array_b, beta_sub_array_b, beta_upper_array_b, color='yellow')
# set plot square
ax2.set_aspect('equal')

## setting for axes
# range
ax2.set_xlim(gamma_min, gamma_max)
ax2.set_ylim(beta_min, beta_max)
# step for ticks
step=2 # stepsize
ax2.set_xticks(np.arange(gamma_min, gamma_max+1/step, step))
ax2.set_yticks(np.arange(beta_min, beta_max+1/step, step))
# label
ax2.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=16)
ax2.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=16)

## set the object
ax3 = fig.add_subplot(2,2,3)
# caption
ax3.text(
    -0.3, 1.12, "(c)",
    transform=ax3.transAxes,
    fontsize=18, va='top', ha='left'
)

# plot lines
ax3.plot(gamma_array_c, beta_upper_array_c, color ='orange',linewidth=1)
ax3.plot(gamma_array_c, beta_sub_array_c, color ='orange',linewidth=1)
ax3.axhline(2, linestyle='dotted', color='gray')
# fillin the region
ax3.fill_between(gamma_array_c, beta_sub_array_c, beta_upper_array_c, color='lightsalmon')
# set plot square
ax3.set_aspect('equal')

## setting for axes
# range
ax3.set_xlim(gamma_min, gamma_max)
ax3.set_ylim(beta_min, beta_max)
# step for ticks
step=2 # value of step
ax3.set_xticks(np.arange(gamma_min, gamma_max+1/step, step))
ax3.set_yticks(np.arange(beta_min, beta_max+1/step, step))
# label of axes
ax3.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=16)
ax3.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=16)

# visualize
plt.show()

## save figure(eps style)
fig.savefig('SuppFig1.eps')
plt.close(fig)




