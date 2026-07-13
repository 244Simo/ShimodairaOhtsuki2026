### import library
import matplotlib.pyplot as plt
import numpy as np

### function
## calculate m-value
def f_m(beta,gamma):
    return (beta+gamma)/np.sqrt(gamma)

## calculate f(ξ) value(Eq.13a/Eq.C15a)
# top of f(ξ)
def q_top(t,m,q):
    return -(1+q)*(t**2)+np.sqrt(q)*m*t-1-np.exp(t**2)
# bottom of f(ξ)
def q_bottom(t):
    return 1-np.exp(t**2)
# f(ξ)（top/bottom）
def f_q(t,m,q):
    return q_top(t,m,q)/q_bottom(t)

## calculate (x,y,v)-values using ξ(=t in program)
# calculate x
def f_teaching_ess_x(t):
    return 1/((t**2)+1)
# calculate v(using x)
def f_teaching_ess_v(x,q,gamma):
    return np.sqrt(1/(q*gamma*x*(1-x)))
# calculate y(using x and v)
def f_teaching_ess_y(x,v,beta,gamma):
    return (1/(gamma*(1-v)))*((1/(x*v))-beta)

## calculate non-teaching ESS(ESS in Kobayashi et al 2015)
# calculate v-value on non-teaching ESS
def f_kwo_v(beta,q,v):
    return beta*(1-v)-q-(1-q)*np.exp(beta*v-1)
# calculate x-value using v-value
def f_kwo_x(beta,v):
    return 1/(beta*v)

## calculate z-value(cultural level)
def f_z(beta,gamma,x,y,v):
    return x*v*np.exp((beta+gamma*y*(1-v))*v*(1-x))
## calculate w-value(fitness)
def f_w(beta,gamma,x,y,v, epsilon):
    return 1+epsilon*((1-y)*(1-v)*f_z(beta,gamma,x,y,v))

### parameters
beta_const1=2.5
beta_const2=4
beta_const3=6
gamma_const=20
q_const=0.9
epsilon=0.1
# initial value of t(=ξ in the paper) when adjusting it to verify the existence of a solution
t_inflex_temp=1
# difference in numerical calculation
delta=10**(-4)
delta_q=10**(-4)

#### β=2.5
## parameters related to q(0<q<1)
q_min=0
q_max=1
delta_plot=0.001

### array
## q-value
q_array=[]
### lifetime allocation(life schedule)
## teaching ESS
# x,y,v
x_array_at_var_q=[]
y_array_at_var_q=[]
v_array_at_var_q=[]
# Social Learning(=(1-x)v)
sl_array_at_var_q=[]
# Individual learning(=xv)
il_array_at_var_q=[]
# SocialLearning+IndividulaLearning(v)
learning_array_at_var_q=[]
# SocialLearning+IndividulaLearning+exploitation(=v+(1-v)(1-y))
learning_and_exploit_array_at_var_q=[]
# SocialLearning+IndividualLearning+Teaching(=v+(1-v)y)
learning_and_teaching_array_at_var_q=[]
# z-value(cultural level)
z_array_at_var_q=[]
# w-value(fitness)
w_array_at_var_q=[]

## non-teaching ESS(Kobayashi et al 2015; hereafter, call it "KWO")
# x,y,v(non-teaching ESS)
x_kwo_array_at_var_q=[]
y_kwo_array_at_var_q=[]
v_kwo_array_at_var_q=[]
# Social Learning(=(1-x)v)(non-teaching ESS)
sl_kwo_array_at_var_q=[]
# Individual Learning(=xv)(non-teaching ESS)
il_kwo_array_at_var_q=[]
# Social Learning+Individual Learning(=v)(non-teaching ESS)
learning_kwo_array_at_var_q=[]
# z-value(cultural level)(non-teaching ESS)
z_kwo_array_at_var_q=[]
# w-value(fitnesss)(non-teaching ESS)
w_kwo_array_at_var_q=[]

### calculate each value based on q-value
## q-value is changed within the for loop(0<q<1)
for i in np.arange(q_min, q_max+delta_plot/2, delta_plot):
    # round q-value
    q=np.round(i, 6)
    # calculate m-value
    m_value=f_m(beta_const1,gamma_const)
    # inflex point of f(ξ)
    while f_q(t_inflex_temp,m_value,q)>f_q(t_inflex_temp-delta,m_value,q):
        t_inflex_temp -= delta_q
    t_inflex=t_inflex_temp
    # initialize temporarily variant for next loop
    t_inflex_temp=1

    ## when teaching ESS does not exist
    if f_q(t_inflex,m_value,q)>q:
    # when β<2, asocial ESS(no social learning and no teaching)
        if beta_const1 <= 2:
            x_ess=1
            y_ess=0
            v_ess=0.5
    # when β>2, non-teaching ESS(social learning exists)
        elif beta_const1 > 2:
            v_temp=0
            while f_kwo_v(beta_const1,q,v_temp)>10**(-8):
                v_temp += (delta * 0.01)
                v_ess=v_temp
                x_ess=f_kwo_x(beta_const1,v_ess)
                y_ess=0
    ## when teaching ESS exists
    else:
        # calculate the solution by incrementally adjusting the ξ(=t)-value
        t_temp=t_inflex
        while (q-f_q(t_temp,m_value,q))>10**(-8):
            t_temp += delta
        t_ess_temp= t_temp
        # calculate the candidates for teaching ESS using ξ(=t)-value
        x_ess_temp=f_teaching_ess_x(t_ess_temp)
        v_ess_temp=f_teaching_ess_v(x_ess_temp, q, gamma_const)
        y_ess_temp=f_teaching_ess_y(x_ess_temp, v_ess_temp, beta_const1, gamma_const)

        ## check domain
        # If there is a solution within the domain, then it is an ESS solution
        if x_ess_temp > 0 and x_ess_temp < 1 and y_ess_temp > 0 and y_ess_temp < 1 and v_ess_temp > 0 and v_ess_temp < 1:
            x_ess = x_ess_temp
            y_ess = y_ess_temp
            v_ess = v_ess_temp
        # If not, non-teaching ESS is the ESS
        else:
            v_temp=0
            while f_kwo_v(beta_const1,q,v_temp)>10**(-8):
                v_temp += delta
            v_ess=v_temp
            x_ess=f_kwo_x(beta_const1,v_ess)
            y_ess=0

    z_ess=f_z(beta_const1, gamma_const, x_ess, y_ess, v_ess)
    w_ess=f_w(beta_const1, gamma_const, x_ess, y_ess, v_ess,epsilon)
    
    ## put values into each array(teaching-ESS)
    ## q-value
    q_array.append(q)
    ## x,y,v
    x_array_at_var_q.append(x_ess)
    y_array_at_var_q.append(y_ess)
    v_array_at_var_q.append(v_ess)
    ## lifetime allocation
    # Social Learning(=(1-x)v)
    sl_array_at_var_q.append((1-x_ess)*v_ess)
    # Individual Learning(=xv)
    il_array_at_var_q.append(x_ess*v_ess)
    # Social Learning+Individual Learning(=v)
    learning_array_at_var_q.append(v_ess)
    # Social Learning+Individual Learning+exploitation(=v+(1-v)(1-y))
    learning_and_exploit_array_at_var_q.append(v_ess+(1-v_ess)*(1-y_ess))
    # Social Learning+Individual Learning+Teaching(=v+(1-v)y)
    learning_and_teaching_array_at_var_q.append(v_ess+(1-v_ess)*y_ess)
    ## z-value(cultural level)
    z_array_at_var_q.append(z_ess)
    ## w-value(fitness)
    w_array_at_var_q.append(w_ess)

    ## calculate non-teaching ESS
    if beta_const1 <= 2:
        x_kwo_ess=1
        y_kwo_ess=0
        v_kwo_ess=0.5
    elif beta_const1 >2:
        v_temp=0
        while f_kwo_v(beta_const1,q,v_temp)>10**(-7):
            v_temp += delta * 0.01
            v_kwo_ess=v_temp
            x_kwo_ess=f_kwo_x(beta_const1,v_kwo_ess)
            y_kwo_ess=0
    z_kwo_ess=f_z(beta_const1, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess)
    w_kwo_ess=f_w(beta_const1, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess,epsilon)

    ## put values into each array(non-teaching ESS)
    ## x,y,v values
    x_kwo_array_at_var_q.append(x_kwo_ess)
    y_kwo_array_at_var_q.append(y_kwo_ess)
    v_kwo_array_at_var_q.append(v_kwo_ess)
    ## lifetime allocation
    # Social Learning(=(1-x)v)
    sl_kwo_array_at_var_q.append((1-x_kwo_ess)*v_kwo_ess)
    # Individual Learning(=xv)
    il_kwo_array_at_var_q.append(x_kwo_ess*v_kwo_ess)
    # Social Learning+Individual Learning(=v)
    learning_kwo_array_at_var_q.append(v_kwo_ess)
    ## z-value(cultural level)
    z_kwo_array_at_var_q.append(z_kwo_ess)
    ## w-value(fitness)
    w_kwo_array_at_var_q.append(w_kwo_ess)

### visualization
## set object for plot
fig=plt.figure(figsize=(12,12))
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, hspace=0.4, wspace=0.4)

## lifetime allocation
ax11 = fig.add_subplot(3,3,1)
# caption
ax11.text(
    -0.3, 1.12, "(a) Life schedule ($β$=2.5)",
    transform=ax11.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax11.plot(q_array, sl_array_at_var_q, color='black')
ax11.plot(q_array, learning_array_at_var_q, color='black')
ax11.plot(q_array, learning_and_teaching_array_at_var_q, color='black')
# paint in different colors
ax11.fill_between(q_array, 0, sl_array_at_var_q, color='white')
ax11.fill_between(q_array, sl_array_at_var_q, learning_array_at_var_q,color='lightgray')
ax11.fill_between(q_array, learning_array_at_var_q, learning_and_teaching_array_at_var_q, color='lightgreen')
ax11.fill_between(q_array, learning_and_teaching_array_at_var_q,1,color='lightyellow')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax11.plot(q_array, sl_kwo_array_at_var_q, linestyle='dotted',color='gray')
ax11.plot(q_array, learning_kwo_array_at_var_q, linestyle='dotted', color='gray')
# range
ax11.set_xlim(q_min,q_max)
ax11.set_ylim(0,1)
# label
ax11.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax11.set_ylabel('Each investment, $τ$', font='arial',fontsize=14)
ax11.tick_params(axis='both', labelsize=10)
# legend(within panel)
ax11.text(0.05,0.02,"SL",color="black",fontsize=16,transform=ax11.transAxes)
ax11.text(0.05,0.27,"IL",color="white",fontsize=16,transform=ax11.transAxes)
ax11.text(0.05,0.85,"EX",color="gold",fontsize=16,transform=ax11.transAxes)
ax11.text(0.88,0.62,"TE",color="green",fontsize=16,transform=ax11.transAxes)

## z(cultural level)
ax12 = fig.add_subplot(3,3,2)
# caption
ax12.text(
    -0.3, 1.12, "(b) Cultural level, $\it{z}$ ($β$=2.5)",
    transform=ax12.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax12.plot(q_array, z_array_at_var_q,color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax12.plot(q_array, z_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax12.set_xlim(q_min,q_max)
ax12.set_ylim(0,6)
# label
ax12.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax12.set_ylabel('Cultural level, $\it{z}$', font='arial', fontsize=14)
ax12.tick_params(axis='both', labelsize=10)

## w(fitness)
ax13 = fig.add_subplot(3,3,3)
# caption
ax13.text(
    -0.3, 1.12, "(c) Fitness, $\it{w}$ ($β$=2.5)",
    transform=ax13.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax13.plot(q_array,w_array_at_var_q, color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax13.plot(q_array,w_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax13.set_xlim(q_min,q_max)
ax13.set_ylim(1,1.2)
# label
ax13.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax13.set_ylabel('Fitness, $\it{w}$', font='arial', fontsize=14)
ax13.tick_params(axis='both', labelsize=10)


#### β=4
## parameters related to q(0<q<1)
q_min=0
q_max=1
delta_plot=0.001

#### array
## q-value
q_array=[]

### lifetime allocation(life schedule)
## teaching ESS
# x,y,v
x_array_at_var_q=[]
y_array_at_var_q=[]
v_array_at_var_q=[]
# Social Learning(=(1-x)v)
sl_array_at_var_q=[]
# Individual learning(=xv)
il_array_at_var_q=[]
# SocialLearning+IndividulaLearning(v)
learning_array_at_var_q=[]
# SocialLearning+IndividulaLearning+exploitation(=v+(1-v)(1-y))
learning_and_exploit_array_at_var_q=[]
# SocialLearning+IndividualLearning+Teaching(=v+(1-v)y)
learning_and_teaching_array_at_var_q=[]
# z-value(cultural level)
z_array_at_var_q=[]
# w-value(fitness)
w_array_at_var_q=[]

## non-teaching ESS(Kobayashi et al 2015; hereafter, call it "KWO")
# x,y,v(non-teaching ESS)
x_kwo_array_at_var_q=[]
y_kwo_array_at_var_q=[]
v_kwo_array_at_var_q=[]
# Social Learning(=(1-x)v)(non-teaching ESS)
sl_kwo_array_at_var_q=[]
# Individual Learning(=xv)(non-teaching ESS)
il_kwo_array_at_var_q=[]
# Social Learning+Individual Learning(=v)(non-teaching ESS)
learning_kwo_array_at_var_q=[]
# z-value(cultural level)(non-teaching ESS)
z_kwo_array_at_var_q=[]
# w-value(fitnesss)(non-teaching ESS)
w_kwo_array_at_var_q=[]

### calculate each value based on q-value
## q-value is changed within the for loop(0<q<1)
for i in np.arange(q_min, q_max+delta_plot/2, delta_plot):
    # round q-value
    q=np.round(i, 5)
    # calculate m-value
    m_value=f_m(beta_const2,gamma_const)
    # inflex point of f(ξ)
    while f_q(t_inflex_temp,m_value,q)>f_q(t_inflex_temp-delta,m_value,q):
        t_inflex_temp-=delta_q
    t_inflex=t_inflex_temp
    # initialize temporarily variant for next loop
    t_inflex_temp=1

    ## when teaching ESS does not exist
    if f_q(t_inflex,m_value,q)>q:
    # when β<2, asocial ESS(no social learning and no teaching)
        if beta_const2 <= 2:
            x_ess=1
            y_ess=0
            v_ess=0.5
    # when β>2, non-teaching ESS(social learning exists)
        elif beta_const2 > 2:
            v_temp=0
            while f_kwo_v(beta_const2,q,v_temp)>10**(-6):
                v_temp += delta
                v_ess=v_temp
                x_ess=f_kwo_x(beta_const2,v_ess)
                y_ess=0
    ## when teaching ESS exists
    else:
        # calculate the solution by incrementally adjusting the ξ(=t)-value
        t_temp=t_inflex
        while (q-f_q(t_temp,m_value,q))>10**(-4):
            t_temp += delta
        t_ess_temp= t_temp
        # calculate the candidates for teaching ESS using ξ(=t)-value
        x_ess_temp=f_teaching_ess_x(t_ess_temp)
        v_ess_temp=f_teaching_ess_v(x_ess_temp, q, gamma_const)
        y_ess_temp=f_teaching_ess_y(x_ess_temp, v_ess_temp, beta_const2, gamma_const)

        ## check domain
        # If there is a solution within the domain, then it is an ESS solution
        if x_ess_temp > 0 and x_ess_temp < 1 and y_ess_temp > 0 and y_ess_temp < 1 and v_ess_temp > 0 and v_ess_temp < 1:
            x_ess = x_ess_temp
            y_ess = y_ess_temp
            v_ess = v_ess_temp
        # If not, non-teaching ESS is the ESS
        else:
            v_temp=0
            while f_kwo_v(beta_const2,q,v_temp)>10**(-6):
                v_temp += delta
            v_ess=v_temp
            x_ess=f_kwo_x(beta_const2,v_ess)
            y_ess=0

    z_ess=f_z(beta_const2, gamma_const, x_ess, y_ess, v_ess)
    w_ess=f_w(beta_const2, gamma_const, x_ess, y_ess, v_ess,epsilon)
    
    ## put values into each array(teaching ESS)
    ## q-value
    q_array.append(q)
    ## x,y,v
    x_array_at_var_q.append(x_ess)
    y_array_at_var_q.append(y_ess)
    v_array_at_var_q.append(v_ess)
    ## lifetime allocation
    # Social Learning(=(1-x)v)
    sl_array_at_var_q.append((1-x_ess)*v_ess)
    # Individual Learning(=xv)
    il_array_at_var_q.append(x_ess*v_ess)
    # Social Learning+Individual Learning(=v)
    learning_array_at_var_q.append(v_ess)
    # Social Learning+Individual Learning+exploitation(=v+(1-v)(1-y))
    learning_and_exploit_array_at_var_q.append(v_ess+(1-v_ess)*(1-y_ess))
    # Social Learning+Individual Learning+Teaching(=v+(1-v)y)
    learning_and_teaching_array_at_var_q.append(v_ess+(1-v_ess)*y_ess)
    # z-value(cultural level)
    z_array_at_var_q.append(z_ess)
    # w-value(fitness)
    w_array_at_var_q.append(w_ess)

    ## calculate non-teaching ESS
    if beta_const2 <= 2:
        x_kwo_ess=1
        y_kwo_ess=0
        v_kwo_ess=0.5
    elif beta_const2 >2:
        v_temp=0
        while f_kwo_v(beta_const2,q,v_temp)>10**(-4):
            v_temp += delta
            v_kwo_ess=v_temp
            x_kwo_ess=f_kwo_x(beta_const2,v_kwo_ess)
            y_kwo_ess=0
    z_kwo_ess=f_z(beta_const2, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess)
    w_kwo_ess=f_w(beta_const2, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess,epsilon)

    ## put values into each array(non-teaching ESS)
    ## x,y,v values
    x_kwo_array_at_var_q.append(x_kwo_ess)
    y_kwo_array_at_var_q.append(y_kwo_ess)
    v_kwo_array_at_var_q.append(v_kwo_ess)
    ## lifetime allocation
    # Social Learning(=(1-x)v)
    sl_kwo_array_at_var_q.append((1-x_kwo_ess)*v_kwo_ess)
    # Individual Learning(=xv)
    il_kwo_array_at_var_q.append(x_kwo_ess*v_kwo_ess)
    # Social Learning+Individual Learning(=v)
    learning_kwo_array_at_var_q.append(v_kwo_ess)
    ## z-value(cultural level)
    z_kwo_array_at_var_q.append(z_kwo_ess)
    ## w-value(fitness)
    w_kwo_array_at_var_q.append(w_kwo_ess)

### visualization
## lifetime allocation
ax21 = fig.add_subplot(3,3,4)
# caption
ax21.text(
    -0.3, 1.12, "(d) Life schedule ($β$=4)",
    transform=ax21.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax21.plot(q_array, sl_array_at_var_q, color='black')
ax21.plot(q_array, learning_array_at_var_q, color='black')
ax21.plot(q_array, learning_and_teaching_array_at_var_q, color='black')
# paint in different colors
ax21.fill_between(q_array, 0, sl_array_at_var_q, color='white')
ax21.fill_between(q_array, sl_array_at_var_q, learning_array_at_var_q,color='lightgray')
ax21.fill_between(q_array, learning_array_at_var_q, learning_and_teaching_array_at_var_q, color='lightgreen')
ax21.fill_between(q_array, learning_and_teaching_array_at_var_q,1,color='lightyellow')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax21.plot(q_array, sl_kwo_array_at_var_q, linestyle='dotted',color='gray')
ax21.plot(q_array, learning_kwo_array_at_var_q, linestyle='dotted', color='gray')
# range
ax21.set_xlim(q_min,q_max)
ax21.set_ylim(0,1)
# label
ax21.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial',fontsize=14)
ax21.set_ylabel('Each investment, $τ$', font='arial',fontsize=14)
ax21.tick_params(axis='both', labelsize=10)

## z(cultural level)
ax22 = fig.add_subplot(3,3,5)
# caption
ax22.text(
    -0.3, 1.12, "(e) Cultural level, $\it{z}$ ($β$=4)",
    transform=ax22.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax22.plot(q_array, z_array_at_var_q,color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax22.plot(q_array, z_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax22.set_xlim(q_min,q_max)
ax22.set_ylim(0,6)
# label
ax22.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax22.set_ylabel('Cultural level, $\it{z}$', font='arial', fontsize=14)
ax22.tick_params(axis='both', labelsize=10)

## w(fitness)
ax23 = fig.add_subplot(3,3,6)
# caption
ax23.text(
    -0.3, 1.12, "(f) Fitness, $\it{w}$ ($β$=4)",
    transform=ax23.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax23.plot(q_array,w_array_at_var_q, color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax23.plot(q_array,w_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax23.set_xlim(q_min,q_max)
ax23.set_ylim(1,1.2)
# label
ax23.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax23.set_ylabel('Fitness, $\it{w}$', font='arial', fontsize=14)
ax23.tick_params(axis='both', labelsize=10)


### β=6
## parameters related to q(0<q<1)
q_min=0
q_max=1
delta_plot=0.001

#### array
## q-value
q_array=[]

### lifetime allocation(life schedule)
## teaching ESS
# x,y,v
x_array_at_var_q=[]
y_array_at_var_q=[]
v_array_at_var_q=[]
# Social Learning(=(1-x)v)
sl_array_at_var_q=[]
# Individual learning(=xv)
il_array_at_var_q=[]
# SocialLearning+IndividulaLearning(v)
learning_array_at_var_q=[]
# SocialLearning+IndividulaLearning+exploitation(=v+(1-v)(1-y))
learning_and_exploit_array_at_var_q=[]
# SocialLearning+IndividualLearning+Teaching(=v+(1-v)y)
learning_and_teaching_array_at_var_q=[]
# z-value(cultural level)
z_array_at_var_q=[]
# w-value(fitness)
w_array_at_var_q=[]

## non-teaching ESS(Kobayashi et al 2015; hereafter, call it "KWO")
# x,y,v(non-teaching ESS)
x_kwo_array_at_var_q=[]
y_kwo_array_at_var_q=[]
v_kwo_array_at_var_q=[]
# Social Learning(=(1-x)v)(non-teaching ESS)
sl_kwo_array_at_var_q=[]
# Individual Learning(=xv)(non-teaching ESS)
il_kwo_array_at_var_q=[]
# Social Learning+Individual Learning(=v)(non-teaching ESS)
learning_kwo_array_at_var_q=[]
# z-value(cultural level)(non-teaching ESS)
z_kwo_array_at_var_q=[]
# w-value(fitnesss)(non-teaching ESS)
w_kwo_array_at_var_q=[]

### calculate each value based on q-value
## q-value is changed within the for loop(0<q<1)
for i in np.arange(q_min, q_max+delta_plot/2, delta_plot):
    # round q-value
    q=np.round(i, 5)
    # calculate m-value
    m_value=f_m(beta_const3,gamma_const)
    # inflex point of f(ξ)
    while f_q(t_inflex_temp,m_value,q)>f_q(t_inflex_temp-delta,m_value,q):
        t_inflex_temp-=delta_q
    t_inflex=t_inflex_temp
    # initialize temporarily variant for next loop
    t_inflex_temp=1

    ## when teaching ESS does not exist
    if f_q(t_inflex,m_value,q)>q:
    # when β<2, asocial ESS(no social learning and no teaching)
        if beta_const3 <= 2:
            x_ess=1
            y_ess=0
            v_ess=0.5
    # when β>2, non-teaching ESS(social learning exists)
        elif beta_const3 > 2:
            v_temp=0
            while f_kwo_v(beta_const3,q,v_temp)>10**(-6):
                v_temp += delta
                v_ess=v_temp
                x_ess=f_kwo_x(beta_const3,v_ess)
                y_ess=0
    ## when teaching ESS exists
    else:
        # calculate the solution by incrementally adjusting the ξ(=t)-value
        t_temp=t_inflex
        while (q-f_q(t_temp,m_value,q))>10**(-4):
            t_temp += delta
        t_ess_temp= t_temp
        # calculate the candidates for teaching ESS using ξ(=t)-value
        x_ess_temp=f_teaching_ess_x(t_ess_temp)
        v_ess_temp=f_teaching_ess_v(x_ess_temp, q, gamma_const)
        y_ess_temp=f_teaching_ess_y(x_ess_temp, v_ess_temp, beta_const3, gamma_const)

        ## check domain
        # If there is a solution within the domain, then it is an ESS solution
        if x_ess_temp > 0 and x_ess_temp < 1 and y_ess_temp > 0 and y_ess_temp < 1 and v_ess_temp > 0 and v_ess_temp < 1:
            x_ess = x_ess_temp
            y_ess = y_ess_temp
            v_ess = v_ess_temp
        # If not, non-teaching ESS is the ESS
        else:
            v_temp=0
            while f_kwo_v(beta_const3,q,v_temp)>10**(-6):
                v_temp += delta
            v_ess=v_temp
            x_ess=f_kwo_x(beta_const3,v_ess)
            y_ess=0

    z_ess=f_z(beta_const3, gamma_const, x_ess, y_ess, v_ess)
    w_ess=f_w(beta_const3, gamma_const, x_ess, y_ess, v_ess,epsilon)
    
    ## put values into each array(teaching ESS)
    ## q-value
    q_array.append(q)
    ## x,y,v
    x_array_at_var_q.append(x_ess)
    y_array_at_var_q.append(y_ess)
    v_array_at_var_q.append(v_ess)
    ## lifetime allocation
    # Social Learning(=(1-x)v)
    sl_array_at_var_q.append((1-x_ess)*v_ess)
    # Individual=Learning(=xv)
    il_array_at_var_q.append(x_ess*v_ess)
    # Social Learning+Individual Learning(=v)
    learning_array_at_var_q.append(v_ess)
    # Social Learning+Individual Learning+exploitation(=v+(1-v)(1-y))
    learning_and_exploit_array_at_var_q.append(v_ess+(1-v_ess)*(1-y_ess))
    # Social Learning+Individual Learning+Teaching(=v+(1-v)y)
    learning_and_teaching_array_at_var_q.append(v_ess+(1-v_ess)*y_ess)
    ## z-value(cultural level)
    z_array_at_var_q.append(z_ess)
    ## w-value(fitness)
    w_array_at_var_q.append(w_ess)

    ## calculate non-teaching ESS
    if beta_const3 <= 2:
        x_kwo_ess=1
        y_kwo_ess=0
        v_kwo_ess=0.5
    elif beta_const3 >2:
        v_temp=0
        while f_kwo_v(beta_const3,q,v_temp)>10**(-4):
            v_temp += delta
            v_kwo_ess=v_temp
            x_kwo_ess=f_kwo_x(beta_const3,v_kwo_ess)
            y_kwo_ess=0
    z_kwo_ess=f_z(beta_const3, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess)
    w_kwo_ess=f_w(beta_const3, gamma_const, x_kwo_ess, y_kwo_ess, v_kwo_ess,epsilon)

    ## put values into each array(non-teaching ESS)
    ## x,y,v values
    x_kwo_array_at_var_q.append(x_kwo_ess)
    y_kwo_array_at_var_q.append(y_kwo_ess)
    v_kwo_array_at_var_q.append(v_kwo_ess)
    ## lifetime allcation
    # Social Learning(=(1-x)v)
    sl_kwo_array_at_var_q.append((1-x_kwo_ess)*v_kwo_ess)
    # Individual Learning(=xv)
    il_kwo_array_at_var_q.append(x_kwo_ess*v_kwo_ess)
    # Social Learning+Individual Learning(=v)
    learning_kwo_array_at_var_q.append(v_kwo_ess)
    ## z-value(cultural level)
    z_kwo_array_at_var_q.append(z_kwo_ess)
    ## w-value(fitness)
    w_kwo_array_at_var_q.append(w_kwo_ess)

### visualization
## lifetime allocation
ax31 = fig.add_subplot(3,3,7)
# caption
ax31.text(
    -0.3, 1.12, "(g) Life schedule ($β$=6)",
    transform=ax31.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax31.plot(q_array, sl_array_at_var_q, color='black')
ax31.plot(q_array, learning_array_at_var_q, color='black')
ax31.plot(q_array, learning_and_teaching_array_at_var_q, color='black')
# paint in different colors
ax31.fill_between(q_array, 0, sl_array_at_var_q, color='white')
ax31.fill_between(q_array, sl_array_at_var_q, learning_array_at_var_q,color='lightgray')
ax31.fill_between(q_array, learning_array_at_var_q, learning_and_teaching_array_at_var_q, color='lightgreen')
ax31.fill_between(q_array, learning_and_teaching_array_at_var_q,1,color='lightyellow')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax31.plot(q_array, sl_kwo_array_at_var_q, linestyle='dotted',color='gray')
ax31.plot(q_array, learning_kwo_array_at_var_q, linestyle='dotted', color='gray')
# range
ax31.set_xlim(q_min,q_max)
ax31.set_ylim(0,1)
# label
ax31.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax31.set_ylabel('Each investment, $τ$', font='arial', fontsize=14)
ax31.tick_params(axis='both', labelsize=10)

## z(cultural level)
ax32 = fig.add_subplot(3,3,8)
# caption
ax32.text(
    -0.3, 1.12, "(h) Cultural level, $\it{z}$ ($β$=6)",
    transform=ax32.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax32.plot(q_array, z_array_at_var_q,color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax32.plot(q_array, z_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax32.set_xlim(q_min,q_max)
ax32.set_ylim(0,6)
# label
ax32.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax32.set_ylabel('Cultural level, $\it{z}$', font='arial', fontsize=14)
ax32.tick_params(axis='both', labelsize=10)

## w(fitness)
ax33 = fig.add_subplot(3,3,9)
# caption
ax33.text(
    -0.3, 1.12, "(i) Fitness, $\it{w}$ ($β$=6)",
    transform=ax33.transAxes,
    fontsize=16, va='top', ha='left'
)
# teaching ESS
ax33.plot(q_array,w_array_at_var_q, color='black')
# non-teaching ESS(ESS in Kobayashi et al 2015)
ax33.plot(q_array,w_kwo_array_at_var_q, linestyle='dotted',color='black')
# range
ax33.set_xlim(q_min,q_max)
ax33.set_ylim(1,1.2)
# label
ax33.set_xlabel('Prob. of vertical transmission, $\it{q}$', font='arial', fontsize=14)
ax33.set_ylabel('Fitness, $\it{w}$', font='arial', fontsize=14)
ax33.tick_params(axis='both', labelsize=10)

# visualize
plt.show()

# save figure(eps style)
fig.savefig('Fig4.eps')
plt.close(fig)