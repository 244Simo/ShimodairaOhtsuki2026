### import library
import matplotlib.pyplot as plt
import numpy as np

### function
## calculate m-value
def f_m(beta,gamma):
    return (beta+gamma)/np.sqrt(gamma)
## calculate f(ξ) value(Eq.13a/Eq.C15a)
# bottom of f(ξ)
def q_top(t,m,q):
    return -(1+q)*(t**2)+np.sqrt(q)*m*t-1-np.exp(t**2)
# top of f(ξ)
def q_bottom(t):
    return 1-np.exp(t**2)
# f(ξ)（top/bottom）
def f_q(t,m,q):
    return q_top(t,m,q)/q_bottom(t)

## calculate values for teaching ESS from ξ(t)
def f_teaching_ess_x(t):
    return 1/((t**2)+1)
def f_teaching_ess_v(x,q,gamma):
    return np.sqrt(1/(q*gamma*x*(1-x)))
def f_teaching_ess_y(x,v,beta,gamma):
    return (1/(gamma*(1-v)))*((1/(x*v))-beta)

## calculate values for non-teaching ESS
# v-value
def f_kwo_v(beta,q,v):
    return beta*(1-v)-q-(1-q)*np.exp(beta*v-1)
# x-value
def f_kwo_x(beta,v):
    return 1/(beta*v)

### caluculate z & w values(at the equilibrium)
def f_z(beta,gamma,x,y,v):
    return x*v*np.exp((beta+gamma*y*(1-v))*v*(1-x))
def f_w(beta,gamma,x,y,v, epsilon):
    return 1+epsilon*((1-y)*(1-v)*f_z(beta,gamma,x,y,v))

### calculate z value during social learning
def f_z_integral(beta,gamma,t,x_r,y_r,v_r):
    return x_r*v_r*np.exp((beta+gamma*y_r*(1-v_r))*v_r*(1-x_r))*(1-np.exp(-(beta+gamma*y_r*(1-v_r))*t))
def f_z_integral_nonteaching(beta,t,x_r,v_r):
    return x_r*v_r*np.exp(beta*v_r*(1-x_r))*(1-np.exp(-(beta*t)))

### calculate the solution of COS(=global ESS when q=1)
def f_cos_x(beta,gamma):
    return 0.5-np.sqrt(0.25-(4*gamma/((beta+gamma)**2)))
def f_cos_y(beta,gamma):
    return ((gamma+beta)/(gamma-beta))*(0.5+np.sqrt(0.25-(4*gamma/((beta+gamma)**2))))-2*beta/(gamma-beta)
def f_cos_v(beta,gamma):
    return 0.5*(1+beta/gamma)

### parameters
beta_const=4
gamma_const=20
q_var0=1
q_var1=0.9
q_var2=0.8
epsilon=0.1
t_inflex_temp=1
# stepsize for numerical calculation
delta=10**(-4)
delta_q=10**(-4)

### Summary
## 1:calculate each value for teaching ESS/no-teaching ESS
## 2:put values to the array for plot
## 3:visualize

### when q=1
## teaching ESS(COS is a global ESS when q=1)
x0=f_cos_x(beta_const,gamma_const)
y0=f_cos_y(beta_const,gamma_const)
v0=f_cos_v(beta_const,gamma_const)
sl0=(1-x0)*v0
il0=x0*v0
# put values for lifestage
x_axis_sl0=np.arange(0, sl0, 0.001)
x_axis_il0=np.arange(sl0, v0, 0.001)
x_axis_ma0=np.arange(v0, 1, 0.001)

## non-teaching ESS
x0_kwo=1/(beta_const-1)
v0_kwo=1-1/beta_const
sl0_kwo=(1-x0_kwo)*v0_kwo
il0_kwo=x0_kwo*v0_kwo
# put values for lifestage
x_axis_sl0_kwo=np.arange(0, sl0_kwo, 0.001)
x_axis_il0_kwo=np.arange(sl0_kwo, v0_kwo, 0.001)
x_axis_ma0_kwo=np.arange(v0_kwo, 1, 0.001)

### when q=0.9
## teaching ESS
# calculateb m-value
m_value=f_m(beta_const,gamma_const)
# inflex point of f(ξ)
while f_q(t_inflex_temp,m_value,q_var1)>f_q(t_inflex_temp-delta,m_value,q_var1):
        t_inflex_temp-=delta_q
t_inflex=t_inflex_temp
# initialize temporarily variant for next calculationt
t_inflex_temp=1

## calculate the solution by incrementally adjusting the ξ(=t)-value
t_temp=t_inflex
while abs(q_var1-f_q(t_temp,m_value,q_var1))>10**(-4):
    t_temp += delta
t_ess_at_const = t_temp
x1=f_teaching_ess_x(t_ess_at_const)
v1=f_teaching_ess_v(x1,q_var1,gamma_const)
y1=f_teaching_ess_y(x1, v1, beta_const, gamma_const)
sl1=(1-x1)*v1
il1=x1*v1
# put values for lifestage
x_axis_sl1=np.arange(0, sl1, 0.001)
x_axis_il1=np.arange(sl1, v1, 0.001)

## no-teaching ESS
v_temp=0
while f_kwo_v(beta_const,q_var1,v_temp)>10**(-6):
    v_temp += delta
v1_kwo=v_temp
x1_kwo=f_kwo_x(beta_const,v1_kwo)
sl1_kwo=(1-x1_kwo)*v1_kwo
il1_kwo=x1_kwo*v1_kwo
# put values for lifestage
x_axis_sl1_kwo=np.arange(0, sl1_kwo, 0.001)
x_axis_il1_kwo=np.arange(sl1_kwo, v1_kwo, 0.001)

### when q=0.8
## teaching ESS
# calculate m-value
m_value=f_m(beta_const,gamma_const)
# inflex point of f(ξ)
while f_q(t_inflex_temp,m_value,q_var2)>f_q(t_inflex_temp-delta,m_value,q_var2):
        t_inflex_temp-=delta_q
t_inflex=t_inflex_temp
# initialize temporarily variant for next calculationt
t_inflex_temp=1

## calculate the solution by incrementally adjusting the ξ(=t)-value
t_temp=t_inflex
while abs(q_var2-f_q(t_temp,m_value,q_var2))>10**(-4):
    t_temp += delta
t_ess_at_const = t_temp
x2=f_teaching_ess_x(t_ess_at_const)
v2=f_teaching_ess_v(x2,q_var2,gamma_const)
y2=f_teaching_ess_y(x2, v2, beta_const, gamma_const)
sl2=(1-x2)*v2
il2=x2*v2
# put values for lifestage
x_axis_sl2=np.arange(0, sl2, 0.001)
x_axis_il2=np.arange(sl2, v2, 0.001)

## no-teaching ESS
v_temp=0
while f_kwo_v(beta_const,q_var2,v_temp)>10**(-6):
    v_temp += delta
v2_kwo=v_temp
x2_kwo=f_kwo_x(beta_const,v2_kwo)
sl2_kwo=(1-x2_kwo)*v2_kwo
il2_kwo=x2_kwo*v2_kwo
# put values for lifestage
x_axis_sl2_kwo=np.arange(0, sl2_kwo, 0.001)
x_axis_il2_kwo=np.arange(sl2_kwo, v2_kwo, 0.001)


### plot
## setting object
fig = plt.figure(figsize=(8,10))

## q=1
ax0 = fig.add_subplot(3,1,1)
# caption
ax0.text(
    -0.1, 1.2, "(a)$\it{q}$=1",
    transform=ax0.transAxes,
    fontsize=14, va='top', ha='left'
)
# lines
ax0.plot(x_axis_sl0,f_z_integral(beta_const,gamma_const,x_axis_sl0,x0,y0,v0),color='blue')
ax0.plot(x_axis_il0,f_z_integral(beta_const,gamma_const,sl0,x0,y0,v0)+x_axis_il0-sl0,color='red')
ax0.axhline(f_z_integral(beta_const,gamma_const,sl0,x0,y0,v0)+il0, xmin=v0, xmax=1, color='black')
ax0.plot(x_axis_sl0_kwo,f_z_integral_nonteaching(beta_const,x_axis_sl0_kwo,x0_kwo,v0_kwo),color='blue', alpha=0.5, linestyle='dotted')
ax0.plot(x_axis_il0_kwo,f_z_integral_nonteaching(beta_const,sl0_kwo,x0_kwo,v0_kwo)+x_axis_il0_kwo-sl0_kwo,color='red', alpha=0.5, linestyle='dotted')
ax0.axhline(f_z_integral_nonteaching(beta_const,sl0_kwo,x0_kwo,v0_kwo)+il0_kwo, xmin=v0_kwo, xmax=1, color='black', alpha=0.5, linestyle='dotted')
# auxiliary lines
ax0.axvline(sl0, linestyle='dotted', color ='gray') # end of SL(teaching exsits)
ax0.axvline(v0, linestyle='dotted', color ='gray') # end of IL(teaching exists)
ax0.axvline(sl0_kwo, linestyle='dotted', color ='darkgray') # end of SL(teaching does not exist)
ax0.axvline(v0_kwo, linestyle='dotted', color ='darkgray') # end of IL(teaching does not exist)
# range
ax0.set_xlim(0, 1)
ax0.set_ylim(0, 16)
# label
ax0.set_xlabel('Lifetime, $t$', font='arial')
ax0.set_ylabel('Cultural level, $z$', font='arial')

## q=0.9
ax1 = fig.add_subplot(3,1,2)
# caption
ax1.text(
    -0.1, 1.2, "(b)$\it{q}$=0.9",
    transform=ax1.transAxes,
    fontsize=14, va='top', ha='left'
)
ax1.plot(x_axis_sl1,f_z_integral(beta_const,gamma_const,x_axis_sl1,x1,y1,v1),color='blue')
ax1.plot(x_axis_il1,f_z_integral(beta_const,gamma_const,sl1,x1,y1,v1)+x_axis_il1-sl1,color='red')
ax1.axhline(f_z_integral(beta_const,gamma_const,sl1,x1,y1,v1)+il1, xmin=v1, xmax=1, color='black')
ax1.plot(x_axis_sl1_kwo,f_z_integral_nonteaching(beta_const,x_axis_sl1_kwo,x1_kwo,v1_kwo),color='blue', alpha=0.5, linestyle='dotted')
ax1.plot(x_axis_il1_kwo,f_z_integral_nonteaching(beta_const,sl1_kwo,x1_kwo,v1_kwo)+x_axis_il1_kwo-sl1_kwo,color='red', alpha=0.5, linestyle='dotted')
ax1.axhline(f_z_integral_nonteaching(beta_const,sl1_kwo,x1_kwo,v1_kwo)+il1_kwo, xmin=v1_kwo, xmax=1, color='black', alpha=0.5, linestyle='dotted')
# auxiliary lines
ax1.axvline(sl1, linestyle='dotted', color ='gray') # end of SL(teaching exsits)
ax1.axvline(v1, linestyle='dotted', color ='gray') # end of IL(teaching exsits)
ax1.axvline(sl1_kwo, linestyle='dotted', color ='darkgray') # end of SL(teaching does not exist)
ax1.axvline(v1_kwo, linestyle='dotted', color ='darkgray') # # end of IL(teaching does not exist)
# range
ax1.set_xlim(0, 1)
ax1.set_ylim(0, 2.5)
# label
ax1.set_xlabel('Lifetime, $t$', font='arial')
ax1.set_ylabel('Cultural level, $z$', font='arial')

## q=0.8
ax2 = fig.add_subplot(3,1,3)
# caption
ax2.text(
    -0.1, 1.2, "(c)$\it{q}$=0.8",
    transform=ax2.transAxes,
    fontsize=14, va='top', ha='left'
)
ax2.plot(x_axis_sl2,f_z_integral(beta_const,gamma_const,x_axis_sl2,x2,y2,v2),color='blue')
ax2.plot(x_axis_il2,f_z_integral(beta_const,gamma_const,sl2,x2,y2,v2)+x_axis_il2-sl2,color='red')
ax2.axhline(f_z_integral(beta_const,gamma_const,sl2,x2,y2,v2)+il2, xmin=v2, xmax=1, color='black')
ax2.plot(x_axis_sl2_kwo,f_z_integral_nonteaching(beta_const,x_axis_sl2_kwo,x2_kwo,v2_kwo),color='blue', alpha=0.5, linestyle='dotted')
ax2.plot(x_axis_il2_kwo,f_z_integral_nonteaching(beta_const,sl2_kwo,x2_kwo,v2_kwo)+x_axis_il2_kwo-sl2_kwo,color='red', alpha=0.5, linestyle='dotted')
ax2.axhline(f_z_integral_nonteaching(beta_const,sl2_kwo,x2_kwo,v2_kwo)+il2_kwo, xmin=v2_kwo, xmax=1, color='black', alpha=0.5, linestyle='dotted')
# auxiliary lines
ax2.axvline(sl2, linestyle='dotted', color ='gray') # end of SL(teaching exsits)
ax2.axvline(v2, linestyle='dotted', color ='gray') # end of IL(teaching exsits)
ax2.axvline(sl2_kwo, linestyle='dotted', color ='darkgray') # end of SL(teaching does not exist)
ax2.axvline(v2_kwo, linestyle='dotted', color ='darkgray') # end of IL(teaching does not exist)
# range
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1.5)
# label
ax2.set_xlabel('Lifetime, $t$', font='arial')
ax2.set_ylabel('Cultural level, $z$', font='arial')

# arange layout
plt.subplots_adjust(hspace=0.6)

## visualize
plt.show()
# save figure(eps style)
fig.savefig('figure5.eps')
plt.close(fig)


