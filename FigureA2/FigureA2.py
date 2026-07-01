### import library
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv

# 変数の準備
beta_max=20
gamma_max=20
delta_d=0.2 #β、γを変化させる際の差分

label_range=np.arange(0,beta_max+1,4)
label_range_calc=np.arange(0,(beta_max+1)/delta_d,4/delta_d)

# (β,γ)の組に対応したマトリクスをそれぞれ作成
# betaが列、gammaが行

# 変数の準備
w_calc=0  #適応度の計算
w_max=0 #最大の適応度の値を格納する変数
z_opti=0  #最大適応度をとるときのzの値を格納する変数
x_opti=0  #最大適応度をとるときのxの値を格納する変数
y_opti=0  #最大適応度をとるときのyの値を格納する変数
v_opti=0  #最大適応度をとるときのvの値を格納する変数

#二次元配列の準備
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
matrix_te_act=[[0 for j in range(int(gamma_max/delta_d)+1)] for i in range(int(beta_max/delta_d)+1)]

#β、γごとに総当たりのループ（i,j）
for j in range(int(beta_max/delta_d)+1):
 for i in range(int(beta_max/delta_d)+1):
  matrix_bg[j][i]=0.5*(i+1) #動作確認
  beta=delta_d*i
  gamma=delta_d*j
  # (x,y,v)を0~1の間で総当たり的に代入して適応度を計算する
  for x in range(0,21,1):
    for y in range(0,21,1):
      for v in range(0,21,1):
        # 計算用の変数を用意（0~1,0.05刻み）
        x_calc=x/20
        y_calc=y/20
        v_calc=v/20

        # 適応度の計算(解析解の値に応じて解を与える)
        A=(beta+gamma)**2-16*gamma #ルート内部の非負条件判定
        B=beta**2-beta*gamma+2*gamma #上部の境界条件(B≦0の時に解を持つ)
        if A<0:#ルート内部の非負条件を満たさない場合はKWO,2015解
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
        elif A>=0:
          if beta>=3:
            if B<=0:
              x_opti=0.5-np.sqrt(0.25-(4*gamma)/((beta+gamma)**2))
              y_opti=((gamma+beta)/(gamma-beta))*(0.5+np.sqrt(0.25-(4*gamma)/((beta+gamma)**2)))-((2*beta)/(gamma-beta))
              v_opti=0.5*(1+beta/gamma)
              z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
              w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
            elif B>0:
              x_opti=1/(beta-1)
              y_opti=0
              v_opti=1-1/beta
              z_opti=(1/(beta))*np.exp(beta-2)
              w_opti=z_opti*(1-v_opti)
          elif beta<3:
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
              if A ==0:
                x_opti=1
                y_opti=0
                v_opti=0.5
                z_opti=0.5
                w_opti=0.25
              else:
                x_opti=0.5-np.sqrt(0.25-(4*gamma)/((beta+gamma)**2))
                y_opti=((gamma+beta)/(gamma-beta))*(0.5+np.sqrt(0.25-(4*gamma)/((beta+gamma)**2)))-((2*beta)/(gamma-beta))
                v_opti=0.5*(1+beta/gamma)
                z_opti=v_opti*x_opti*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
                w_opti=v_opti*(1-v_opti)*x_opti*(1-y_opti)*np.exp((beta*v_opti*(1-x_opti)+gamma*y_opti*(1-v_opti)*v_opti*(1-x_opti)))
            #端点解の判定(指数関数の数値誤差の影響が大きいので、非負条件の境界のみ)

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

        # 求めた値をマトリクスに格納していく
        matrix_x[j][i]=x_opti
        matrix_y[j][i]=y_opti
        matrix_v[j][i]=v_opti
        matrix_sl[j][i]=(1-x_opti)*v_opti
        matrix_il[j][i]=x_opti*v_opti
        matrix_ex[j][i]=(1-y_opti)*(1-v_opti)
        matrix_te[j][i]=y_opti*(1-v_opti)
        matrix_te_act[j][i]=gamma*y_opti*(1-v_opti)
        matrix_z[j][i]=z_opti
        matrix_w[j][i]=w_opti

  # 格納変数のリセット（これをしないと次のループに前のbeta,gammaの時の値が持ち越されてしまう）
  w_max=0
  w_calc=0
  z_opti=0
  x_opti=0
  y_opti=0
  v_opti=0

## Visualization
# 配列の転置
matrix_x_t=np.array(matrix_x).T
matrix_y_t=np.array(matrix_y).T
matrix_v_t=np.array(matrix_v).T
matrix_z_t=np.array(matrix_z).T
matrix_w_t=np.array(matrix_w).T
matrix_sl_t=np.array(matrix_sl).T
matrix_il_t=np.array(matrix_il).T
matrix_ex_t=np.array(matrix_ex).T
matrix_te_t=np.array(matrix_te).T
matrix_te_act_t=np.array(matrix_te_act).T

# Supplementary Figure3(SL,IL,EX, TE,z,w)
figS3 = plt.figure(figsize=(9,12))
figS3.subplots_adjust(wspace=0.2,hspace=0.4)
# SLのプロット（転置）
ax1 = figS3.add_subplot(3,2,1)
ax1.text(
    -0.25, 1.15, "(a) $τ_{\mathrm{SL}}^{\circ}$",
    transform=ax1.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_sl2 = ax1.imshow(matrix_sl_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax1.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax1.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
# axS1.set_title("SL")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_sl2)
# ILのプロット（転置）
ax2 = figS3.add_subplot(3,2,2)
ax2.text(
    -0.25, 1.15, "(b) $τ_{\mathrm{IL}}^{\circ}$",
    transform=ax2.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_il2 = ax2.imshow(matrix_il_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax2.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax2.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
#axS2.set_title("IL")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_il2)
# EXのプロット（転置）
ax3 = figS3.add_subplot(3,2,3)
ax3.text(
    -0.25, 1.15, "(c) $τ_{\mathrm{EX}}^{\circ}$",
    transform=ax3.transAxes,
    font='DejaVu Sans',  fontsize=16, va='top', ha='left'
)
im_ex2 = ax3.imshow(matrix_ex_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax3.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax3.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
#axS3.set_title("EX")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_ex2)
# TEのプロット（転置）
ax4 = figS3.add_subplot(3,2,4)
ax4.text(
    -0.25, 1.15, "(d) $τ_{\mathrm{TE}}^{\circ}$",
    transform=ax4.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_te2 = ax4.imshow(matrix_te_t, origin='lower',cmap='viridis',vmin=0, vmax=1)
ax4.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax4.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
#ax4.set_title("TE")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_te2)
# zのプロット（転置・対数）
ax5 = figS3.add_subplot(3,2,5)
ax5.text(
    -0.25, 1.15, "(e) ${\mathit{z}}^{\circ}$",
    transform=ax5.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_z2 = ax5.imshow(matrix_z_t, origin='lower',cmap='plasma',norm=LogNorm(vmin=0.1, vmax=1000000))
ax5.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax5.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
#ax.set_title("z at COS")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_z2)
# wのプロット（転置・対数）
ax6 = figS3.add_subplot(3,2,6)
ax6.text(
    -0.25, 1.15, "(f) $w^{\circ}$",
    transform=ax6.transAxes,
    font='DejaVu Sans', fontsize=16, va='top', ha='left'
)
im_w2 = ax6.imshow(matrix_w_t, origin='lower',cmap='plasma',norm=LogNorm(vmin=0.1, vmax=100000))
ax6.set_xlabel('Effect of teaching, $γ$', font='arial', fontsize=10)
ax6.set_ylabel('Efficiency of observation learning, $β$', font='arial', fontsize=10)
#ax.set_title("w at COS")
plt.xticks(label_range_calc,label_range)
plt.yticks(label_range_calc,label_range)
plt.colorbar(im_w2)



#スペースの調整

figS3.savefig('figureS2_COS_value2.eps')


#plt.show()

#plt.close(figS2)
plt.close(figS3)
