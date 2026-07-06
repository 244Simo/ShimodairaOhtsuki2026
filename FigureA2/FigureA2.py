### import library
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

## parameters
beta_max=20
gamma_max=20
delta_d=0.2 #deviation of calculation (β and γ)
# array for plotting
label_range=np.arange(0,beta_max+1,4)
label_range_calc=np.arange(0,(beta_max+1)/delta_d,4/delta_d)

# variant
w_calc=0  #calculation for fitness
w_max=0   #max fitness
z_opti=0  #z value(cultural level) when fitness is max
x_opti=0  #x value when fitness is max
y_opti=0  #y value when fitness is max
v_opti=0  #v value when fitness is max

## matrix for each (β,γ)
# beta is row, gamma is column
# 2D matrix
matrix_bg=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]  #動作確認用
matrix_x=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_y=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_v=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_z=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_w=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_sl=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_il=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_ex=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]
matrix_te=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]

### caluculation for each β & γ（i,j）
for j in range(int(gamma_max/delta_d)+1):
 for i in range(int(beta_max/delta_d)+1):
  matrix_bg[j][i]=0.5*(i+1)
  beta=delta_d*i
  gamma=delta_d*j
  ## calculate fitness value
  ## at the parameters where two candidates for COS exist, compare fitness value of plausible pairs of (x,y,v)[stepsize:0.05]
  for x in range(0,21,1):
    for y in range(0,21,1):
      for v in range(0,21,1):
        # variants for calculation（0~1,stepsize:0.05）
        x_calc=x/20
        y_calc=y/20
        v_calc=v/20

        #3 calculate fitness value
        A=(beta+gamma)**2-16*gamma #boundary the condition for where teaching COS exists
        B=beta**2-beta*gamma+2*gamma #boundary where only teaching COS exist
        ## When teaching COS does not exsit
        if A<0:
          if beta <= 2:
            x_opti=1
            y_opti=0
            v_opti=0.5
            z_opti=0.5
            w_opti=0.25
          elif beta>2:
            x_opti=1/(beta-1)
            y_opti=0
            v_opti=1-1/beta
            z_opti=(1/(beta))*np.exp(beta-2)
            w_opti=z_opti*(1-v_opti)
        ## When teaching COS exsits
        elif A>=0:
          if beta>=3:
            ## only teaching COS exist
            if B<=0:
              x_opti=0.5-np.sqrt(0.25-(4*gamma)/((beta+gamma)**2))
              y_opti=((gamma+beta)/(gamma-beta))*(0.5+np.sqrt(0.25-(4*gamma)/((beta+gamma)**2)))-((2*beta)/(gamma-beta))
              v_opti=0.5*(1+beta/gamma)
              z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
              w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
            ## teaching COS does not exist
            elif B>0:
              x_opti=1/(beta-1)
              y_opti=0
              v_opti=1-1/beta
              z_opti=(1/(beta))*np.exp(beta-2)
              w_opti=z_opti*(1-v_opti)
          elif beta<3:
            ## teaching COS does not exist
            if gamma<9:
              if beta>2:
                x_opti=1/(beta-1)
                y_opti=0
                v_opti=1-1/beta
                z_opti=(1/(beta))*np.exp(beta-2)
                w_opti=z_opti*(1-v_opti)
              elif beta<=2:
                x_opti=1
                y_opti=0
                v_opti=0.5
                z_opti=0.5
                w_opti=0.25
            elif gamma>=9:
              ## boundary
              if A ==0:
                x_opti=1
                y_opti=0
                v_opti=0.5
                z_opti=0.5
                w_opti=0.25
              ## teaching COS exist
              else:
                x_opti=0.5-np.sqrt(0.25-(4*gamma)/((beta+gamma)**2))
                y_opti=((gamma+beta)/(gamma-beta))*(0.5+np.sqrt(0.25-(4*gamma)/((beta+gamma)**2)))-((2*beta)/(gamma-beta))
                v_opti=0.5*(1+beta/gamma)
                z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
                w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))

            ## calculate non-teaching COS and compare teaching COS
            w_calc=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
            if beta>2:
              w_kwo=(1/(beta**2))*np.exp(beta-2)
              if w_kwo>w_calc:
                x_opti=1/(beta-1)
                y_opti=0
                v_opti=1-1/beta
                z_opti=(1/(beta))*np.exp(beta-2)
                w_opti=z_opti*(1-v_opti)
              else:
                z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
                w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
            elif beta<=2:
              w_kwo=0.25
              if w_kwo>w_calc:
                x_opti=1
                y_opti=0
                v_opti=0.5
                z_opti=0.5
                w_opti=0.25
              else:
                z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
                w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))

        # put results into arrays
        matrix_x[j][i]=x_opti
        matrix_y[j][i]=y_opti
        matrix_v[j][i]=v_opti
        matrix_sl[j][i]=(1-x_opti)*v_opti
        matrix_il[j][i]=x_opti*v_opti
        matrix_ex[j][i]=(1-y_opti)*(1-v_opti)
        matrix_te[j][i]=y_opti*(1-v_opti)
        matrix_z[j][i]=z_opti
        matrix_w[j][i]=w_opti

  # reset varinants before moving on to next loop
  w_max=0
  w_calc=0
  z_opti=0
  x_opti=0
  y_opti=0
  v_opti=0

## prepare for visualization
# transposed matrix(for plot)
matrix_x_t=np.array(matrix_x).T
matrix_y_t=np.array(matrix_y).T
matrix_v_t=np.array(matrix_v).T
matrix_z_t=np.array(matrix_z).T
matrix_w_t=np.array(matrix_w).T
matrix_sl_t=np.array(matrix_sl).T
matrix_il_t=np.array(matrix_il).T
matrix_ex_t=np.array(matrix_ex).T
matrix_te_t=np.array(matrix_te).T

### plot
## setting object
figS = plt.figure(figsize=(9,12))
figS.subplots_adjust(wspace=0.2,hspace=0.4)

# plot for SL
ax1 = figS.add_subplot(3,2,1)
ax1.text(
    -0.25, 1.15, "(a) $τ_{\mathrm{SL}}^{\circ}$",
    transform=ax1.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_sl2 = ax1.imshow(matrix_sl_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax1.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax1.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_sl2)

# plot for IL
ax2 = figS.add_subplot(3,2,2)
ax2.text(
    -0.25, 1.15, "(b) $τ_{\mathrm{IL}}^{\circ}$",
    transform=ax2.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_il2 = ax2.imshow(matrix_il_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax2.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax2.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_il2)

# plot for EX
ax3 = figS.add_subplot(3,2,3)
ax3.text(
    -0.25, 1.15, "(c) $τ_{\mathrm{EX}}^{\circ}$",
    transform=ax3.transAxes,
    font='DejaVu Sans',  fontsize=16, va='top', ha='left'
)
im_ex2 = ax3.imshow(matrix_ex_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax3.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax3.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_ex2)

# plot for TE
ax4 = figS.add_subplot(3,2,4)
ax4.text(
    -0.25, 1.15, "(d) $τ_{\mathrm{TE}}^{\circ}$",
    transform=ax4.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_te2 = ax4.imshow(matrix_te_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax4.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax4.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_te2)

# plot for z
ax5 = figS.add_subplot(3,2,5)
ax5.text(
    -0.25, 1.15, "(e) ${\mathit{z}}^{\circ}$",
    transform=ax5.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_z2 = ax5.imshow(matrix_z_t, origin='lower',cmap='plasma',norm=LogNorm(vmin=0.1, vmax=1000000))
ax5.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax5.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_z2)

# plot for w
ax6 = figS.add_subplot(3,2,6)
ax6.text(
    -0.25, 1.15, "(f) $w^{\circ}$",
    transform=ax6.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_w2 = ax6.imshow(matrix_w_t, origin='lower',cmap='plasma',norm=LogNorm(vmin=0.1, vmax=100000))
ax6.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax6.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_w2)

## visualization
plt.show()

# save figure(eps style)
figS.savefig('figureA2.eps')
plt.close(figS)
