### import library
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

### function
## calculate β base on Eq.13d given q and γ
def f1(m,x):
    return -1*x+m*np.sqrt(x)

## the boundary of only teaching ESS(Eq.C.24 in supplementary material)
def g1(beta,gamma):
    return beta-1-(beta**2)/(q*gamma)-q-(1-q)*np.exp((beta**2)/(q*gamma))

## calculate f(ξ) value(Eq.13a/Eq.C15a)
# top of f(ξ)
def q_top(t,m):
    return -(1+q)*(t**2)+np.sqrt(q)*m*t-1-np.exp(t**2)
# bottom of f(ξ)
def q_bottom(t):
    return 1-np.exp(t**2)
# f(ξ)（top/bottom）
def f_q(t,m):
    return q_top(t,m)/q_bottom(t)

### panel(a)q=1
## parameter setting
q=1
# range of γ（x-axis）
xmax=30
xmin=0
# range of β（y-axis）
ymax=30
ymin=0

## array for plot
# teaching-ESS(light green/dark cyan/light cyan)
beta_array_1=[]
gamma_array_1=[]
border_array_1=[]
# only teaching-ESS(light green)
beta_array_2=[]
gamma_array_2=[]

### obtain the border of teaching-ESS(Eq.13a/Eq.C15a)
## aim: calculate the value of m such that the line y=q is tangent to the graph of f(ξ) at its inflection point
## (q: given)
## t in this code is ξ in the manuscript
#-----
# 1.set starting point, m_temp & ξ_temp(=t)
# 2.calculate ξ-value of the inflection point
# 3.if f(ξ) at inflex point > q, update m_temp value
# 4:continue 2 and 3
# note: the setting of calculation uses the property of y= f(ξ)(see also Figure C.2 in supplementary material)
#-----

## parameter for calcutation
delta=10**(-8)
epsilon=10**(-6)
m_value_temp=4
t_inflex_temp=1
q_inflex_temp=f_q(t_inflex_temp,m_value_temp)

## calculate m-value
while q_inflex_temp>q:
    while q_inflex_temp>f_q(t_inflex_temp-delta,m_value_temp):
        t_inflex_temp-=delta
        q_inflex_temp=f_q(t_inflex_temp,m_value_temp)
    m_value_temp+=epsilon
# record the value
t_inflex=t_inflex_temp
m_value=m_value_temp
q_inflex=f_q(t_inflex,m_value)

## calculate the γ value of intersection of teaching-ESS region and β=2(need for plotting)
# γ can be caluculate by the definition of m(m=(β+γ)/sqrt(γ)), see Eq.13c
gamma_intersection=((m_value**2-4)+np.sqrt(m_value**4-8*m_value**2))/2

### obtain the border of only-teaching ESS(light green)
## solve the Eq.C.24(supplementary material)
# set initial value
beta=2
gamma=40
loopcount=0 # for judge the inflection point
delta_ess=10**(-3)

# calculate the border
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
    loopcount+=1  # increase loop count(If the condition is not met no matter how many times the loop is repeated, it is determined that the inflection point has been passed)
  # upper half
  if loopcount>30:
    gamma+=(loopcount-30)*delta_ess
    beta+=delta_ess
  # record the border value
  if border_value1>0:
    beta_array_2.append(beta)
    gamma_array_2.append(gamma)
    border_value2=g1(beta,gamma)
    loopcount=0   # reset the loop count
    gamma-=delta_ess
  # end the while loop
  if beta>40:
    beta =2
    gamma =40
    loopcount =0
    break

## obtain the boundary of only-teaching ESS region(lightgreen) and non-teaching ESS coexist region(dark cyan)
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
  if border_value1>0:
    beta_array_1.append(beta)
    gamma_array_1.append(gamma)
    border_array_1.append(f1(m_value,gamma))
    border_value2=g1(beta,gamma)
    gamma-=delta_ess
  # end the wile loop(β=3 is the possibole maximum)
  if beta>3:
    break

### plot
## legend will be set at the end of the program
# setting object
fig = plt.figure(figsize=(8,8))
fig.subplots_adjust(hspace=0.4, wspace=0.4)
ax1 = fig.add_subplot(2,2,1)
ax1.set_aspect('equal')
# caption
ax1.text(
    -0.3, 1.12, "(a)$\it{q}$=1",
    transform=ax1.transAxes,
    fontsize=14, va='top', ha='left'
)
# fill each region
ax1.axhspan(0, 2, color='gray') # only asocial ESS
ax1.fill_between(gamma_array_1,beta_array_1,border_array_1, fc="darkcyan")  # non-teaching ESS and teaching ESS coexist
ax1.fill_between(gamma_array_1,2,border_array_1,where=gamma_array_1>gamma_intersection,fc="paleturquoise")  # asocial ESS and teaching ESS coexist
ax1.fill_between(gamma_array_2,beta_array_2,0,fc="lightgreen", alpha=1) # only teaching ESS
# range of each axis
ax1.set_xlim(0, 30)
ax1.set_ylim(0, 30)
# label
ax1.set_xlabel('Effect of teaching, $γ$', font='arial')
ax1.set_ylabel('Efficiency of observation learning, $β$', font='arial')


### panel(b)q=0.95
## parameter setting
q=0.95
# range of γ（x-axis）
xmax=30
xmin=0
# range of β（y-axis）
ymax=30
ymin=0

## array for plot
# teaching-ESS(light green/dark cyan/light cyan)
beta_array_1=[]
gamma_array_1=[]
border_array_1=[]
# only teaching-ESS(light green)
beta_array_2=[]
gamma_array_2=[]

### obtain the border of teaching-ESS(Eq.13a/Eq.C15a)
## aim: calculate the value of m such that the line y=q is tangent to the graph of f(ξ) at its inflection point
## (q: given)
## t in this code is ξ in the manuscript
#-----
# 1.set starting point, m_temp & ξ_temp(=t)
# 2.calculate ξ-value of the inflection point
# 3.if f(ξ) at inflex point > q, update m_temp value
# 4:continue 2 and 3
# note: the setting of calculation uses the property of y= f(ξ)(see also Figure C.2 in supplementary material)
#-----

## parameter for calcutation
delta=10**(-8)
epsilon=10**(-6)
m_value_temp=4
t_inflex_temp=1
q_inflex_temp=f_q(t_inflex_temp,m_value_temp)

## calculate m-value
while q_inflex_temp>q:
    while q_inflex_temp>f_q(t_inflex_temp-delta,m_value_temp):
        t_inflex_temp-=delta
        q_inflex_temp=f_q(t_inflex_temp,m_value_temp)
    m_value_temp+=epsilon
# record the value
t_inflex=t_inflex_temp
m_value=m_value_temp
q_inflex=f_q(t_inflex,m_value)

## calculate the γ value of intersection of teaching-ESS region and β=2(need for plotting)
# γ can be caluculate by the definition of m(m=(β+γ)/sqrt(γ)), see Eq.13c
gamma_intersection=((m_value**2-4)+np.sqrt(m_value**4-8*m_value**2))/2

### obtain the border of only-teaching ESS(light green)
## solve the Eq.C.24(supplementary material)
# set initial value
beta=2
gamma=40
loopcount=0 # for judge the inflection point
delta_ess=10**(-3)

# calculate the border
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
    loopcount+=1  # increase loop count(If the condition is not met no matter how many times the loop is repeated, it is determined that the inflection point has been passed)
  # upper half
  if loopcount>30:
    gamma+=(loopcount-30)*delta_ess
    beta+=delta_ess
  # record the border value
  if border_value1>0:
    beta_array_2.append(beta)
    gamma_array_2.append(gamma)
    border_value2=g1(beta,gamma)
    loopcount=0   # reset the loop count
    gamma-=delta_ess
  # end the while loop
  if beta>40:
    beta =2
    gamma =40
    loopcount =0
    break

## obtain the boundary of only-teaching ESS region(lightgreen) and non-teaching ESS coexist region(dark cyan)
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
  if border_value1>0:
    beta_array_1.append(beta)
    gamma_array_1.append(gamma)
    border_array_1.append(f1(m_value,gamma))
    border_value2=g1(beta,gamma)
    gamma-=delta_ess
  # end the wile loop(β=3 is the possibole maximum)
  if beta>3:
    break

### plot
## legend will be set at the end of the program
# setting object
ax2 = fig.add_subplot(2,2,2)
ax2.set_aspect('equal')  # グラフの縦横比を同じにするための関数
# caption
ax2.text(
    -0.3, 1.12, "(b)$\it{q}$=0.95",
    transform=ax2.transAxes,
    fontsize=14, va='top', ha='left'
)
# fill each region
ax2.axhspan(0, 2, color='gray') # only asocial ESS
ax2.fill_between(gamma_array_1,beta_array_1,border_array_1, fc="darkcyan")  # non-teaching ESS and teaching ESS coexist
ax2.fill_between(gamma_array_1,2,border_array_1,where=gamma_array_1>gamma_intersection,fc="paleturquoise")  # asocial ESS and teaching ESS coexist
ax2.fill_between(gamma_array_2,beta_array_2,0,fc="lightgreen", alpha=1) # only teaching ESS
# range of each axis
ax2.set_xlim(0, 30)
ax2.set_ylim(0, 30)
# label
ax2.set_xlabel('Effect of teaching, $γ$', font='arial')
ax2.set_ylabel('Efficiency of observation learning, $β$', font='arial')


### panel(c)q=0.9
## parameter setting
q=0.9
# range of γ（x-axis）
xmax=30
xmin=0
# range of β（y-axis）
ymax=30
ymin=0

## array for plot
# teaching-ESS(light green/dark cyan/light cyan)
beta_array_1=[]
gamma_array_1=[]
border_array_1=[]
# only teaching-ESS(light green)
beta_array_2=[]
gamma_array_2=[]

### obtain the border of teaching-ESS(Eq.13a/Eq.C15a)
## aim: calculate the value of m such that the line y=q is tangent to the graph of f(ξ) at its inflection point
## (q: given)
## t in this code is ξ in the manuscript
#-----
# 1.set starting point, m_temp & ξ_temp(=t)
# 2.calculate ξ-value of the inflection point
# 3.if f(ξ) at inflex point > q, update m_temp value
# 4:continue 2 and 3
# note: the setting of calculation uses the property of y= f(ξ)(see also Figure C.2 in supplementary material)
#-----

## parameter for calcutation
delta=10**(-8)
epsilon=10**(-6)
m_value_temp=4
t_inflex_temp=1
q_inflex_temp=f_q(t_inflex_temp,m_value_temp)

## calculate m-value
while q_inflex_temp>q:
    while q_inflex_temp>f_q(t_inflex_temp-delta,m_value_temp):
        t_inflex_temp-=delta
        q_inflex_temp=f_q(t_inflex_temp,m_value_temp)
    m_value_temp+=epsilon
# record the value
t_inflex=t_inflex_temp
m_value=m_value_temp
q_inflex=f_q(t_inflex,m_value)

## calculate the γ value of intersection of teaching-ESS region and β=2(need for plotting)
# γ can be caluculate by the definition of m(m=(β+γ)/sqrt(γ)), see Eq.13c
gamma_intersection=((m_value**2-4)+np.sqrt(m_value**4-8*m_value**2))/2

### obtain the border of only-teaching ESS(light green)
## solve the Eq.C.24(supplementary material)
# set initial value
beta=2
gamma=40
loopcount=0 # for judge the inflection point
delta_ess=10**(-3)

# calculate the border
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
    loopcount+=1  # increase loop count(If the condition is not met no matter how many times the loop is repeated, it is determined that the inflection point has been passed)
  # upper half
  if loopcount>30:
    gamma+=(loopcount-30)*delta_ess
    beta+=delta_ess
  # record the border value
  if border_value1>0:
    beta_array_2.append(beta)
    gamma_array_2.append(gamma)
    border_value2=g1(beta,gamma)
    loopcount=0   # reset the loop count
    gamma-=delta_ess
  # end the while loop
  if beta>40:
    beta =2
    gamma =40
    loopcount =0
    break

## obtain the boundary of only-teaching ESS region(lightgreen) and non-teaching ESS coexist region(dark cyan)
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
  if border_value1>0:
    beta_array_1.append(beta)
    gamma_array_1.append(gamma)
    border_array_1.append(f1(m_value,gamma))
    border_value2=g1(beta,gamma)
    gamma-=delta_ess
  # end the wile loop(β=3 is the possibole maximum)
  if beta>3:
    break

### plot
## legend will be set at the end of the program
# setting object
ax3 = fig.add_subplot(2,2,3)
ax3.set_aspect('equal')  # グラフの縦横比を同じにするための関数
# caption
ax3.text(
    -0.3, 1.12, "(c)$\it{q}$=0.9",
    transform=ax3.transAxes,
    fontsize=14, va='top', ha='left'
)
# fill each region
ax3.axhspan(0, 2, color='gray') # only asocial ESS
ax3.fill_between(gamma_array_1,beta_array_1,border_array_1, fc="darkcyan")  # non-teaching ESS and teaching ESS coexist
ax3.fill_between(gamma_array_1,2,border_array_1,where=gamma_array_1>gamma_intersection,fc="paleturquoise")  # asocial ESS and teaching ESS coexist
ax3.fill_between(gamma_array_2,beta_array_2,0,fc="lightgreen", alpha=1) # only teaching ESS
# range of each axis
ax3.set_xlim(0, 30)
ax3.set_ylim(0, 30)
# label
ax3.set_xlabel('Effect of teaching, $γ$', font='arial')
ax3.set_ylabel('Efficiency of observation learning, $β$', font='arial')

### panel(d)q=0.8
## parameter setting
q=0.8
# range of γ（x-axis）
xmax=30
xmin=0
# range of β（y-axis）
ymax=30
ymin=0

## array for plot
# teaching-ESS(light green/dark cyan/light cyan)
beta_array_1=[]
gamma_array_1=[]
border_array_1=[]
# only teaching-ESS(light green)
beta_array_2=[]
gamma_array_2=[]

### obtain the border of teaching-ESS(Eq.13a/Eq.C15a)
## aim: calculate the value of m such that the line y=q is tangent to the graph of f(ξ) at its inflection point
## (q: given)
## t in this code is ξ in the manuscript
#-----
# 1.set starting point, m_temp & ξ_temp(=t)
# 2.calculate ξ-value of the inflection point
# 3.if f(ξ) at inflex point > q, update m_temp value
# 4:continue 2 and 3
# note: the setting of calculation uses the property of y= f(ξ)(see also Figure C.2 in supplementary material)
#-----

## parameter for calcutation
delta=10**(-8)
epsilon=10**(-6)
m_value_temp=4
t_inflex_temp=1
q_inflex_temp=f_q(t_inflex_temp,m_value_temp)

## calculate m-value
while q_inflex_temp>q:
    while q_inflex_temp>f_q(t_inflex_temp-delta,m_value_temp):
        t_inflex_temp-=delta
        q_inflex_temp=f_q(t_inflex_temp,m_value_temp)
    m_value_temp+=epsilon
# record the value
t_inflex=t_inflex_temp
m_value=m_value_temp
q_inflex=f_q(t_inflex,m_value)

## calculate the γ value of intersection of teaching-ESS region and β=2(need for plotting)
# γ can be caluculate by the definition of m(m=(β+γ)/sqrt(γ)), see Eq.13c
gamma_intersection=((m_value**2-4)+np.sqrt(m_value**4-8*m_value**2))/2

### obtain the border of only-teaching ESS(light green)
## solve the Eq.C.24(supplementary material)
# set initial value
beta=2
gamma=40
loopcount=0 # for judge the inflection point
delta_ess=10**(-3)

# calculate the border
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
    loopcount+=1  # increase loop count(If the condition is not met no matter how many times the loop is repeated, it is determined that the inflection point has been passed)
  # upper half
  if loopcount>30:
    gamma+=(loopcount-30)*delta_ess
    beta+=delta_ess
  # record the border value
  if border_value1>0:
    beta_array_2.append(beta)
    gamma_array_2.append(gamma)
    border_value2=g1(beta,gamma)
    loopcount=0   # reset the loop count
    gamma-=delta_ess
  # end the while loop
  if beta>40:
    beta =2
    gamma =40
    loopcount =0
    break

## obtain the boundary of only-teaching ESS region(lightgreen) and non-teaching ESS coexist region(dark cyan)
while beta>=2:
  border_value1=g1(beta,gamma)
  # lower half
  if border_value1<0:
    beta+=delta_ess
  if border_value1>0:
    beta_array_1.append(beta)
    gamma_array_1.append(gamma)
    border_array_1.append(f1(m_value,gamma))
    border_value2=g1(beta,gamma)
    gamma-=delta_ess
  # end the wile loop(β=3 is the possibole maximum)
  if beta>3:
    break

### plot
## legend will be set at the end of the program
# setting object
ax4 = fig.add_subplot(2,2,4)
ax4.set_aspect('equal')
# caption
ax4.text(
    -0.3, 1.12, "(d)$\it{q}$=0.8",
    transform=ax4.transAxes,
    fontsize=14, va='top', ha='left'
)
# fill each region
ax4.axhspan(0, 2, color='gray') # only asocial ESS
ax4.fill_between(gamma_array_1,beta_array_1,border_array_1, fc="darkcyan")  # non-teaching ESS and teaching ESS coexist
ax4.fill_between(gamma_array_1,2,border_array_1,where=gamma_array_1>gamma_intersection,fc="paleturquoise")  # asocial ESS and teaching ESS coexist
ax4.fill_between(gamma_array_2,beta_array_2,0,fc="lightgreen", alpha=1) # only teaching ESS
# range of each axis
ax4.set_xlim(0, 30)
ax4.set_ylim(0, 30)
# label
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
    bbox_to_anchor=(0.5, 0)
)

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

## plot the figure
plt.show()

## save figure(eps style)
fig.savefig('Fig3.eps', bbox_inches="tight")
plt.close(fig)




















