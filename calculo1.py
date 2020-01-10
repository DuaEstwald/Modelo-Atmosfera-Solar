# Estoy hasta el conio y este es el nosecuantos intentos ya de que esto salga bien

import numpy as np
from poblaciones import *

def data(txt,ldata,lheader):
    d = open(txt,'r')
    dat = np.loadtxt(d.read().splitlines()[ldata:], unpack = True)
    return dat

header = 'k','lgTauR','lgTau5','T','Pe','Pg','Prad','Pturb'
t5000 = data('t5000.dat',7,6)
t8000 = data('t8000.dat',7,6)
nH = 12.00
nHe = 10.93


Pe = t5000[4], t8000[4] #[erg/cm]
T = t5000[3], t8000[3] #[K]
logtau = t5000[1], t8000[1]



gHmin = 1.
XHmin = 0.0

gHI = 2., 8., 18., 32.
XHI = 0.0, 10.1988357, 12.0875051,12.7485392

gHII = 1.
XHII = 0.0

gHeI = 1., 3., 1., 9., 3.
XHeI = 0.0, 19.81961468, 20.61577496, 20.96409651, 21.21802284

gHeII = 2., 8.
XHeII = 0.0, 40.8137567


gHeIII = 1.
XHeIII = 0.0 



X = 0.754, 13.598, 0.0, 24.587, 54.416, 0.0 # fundamentales H-, HI, HII, HeI, HeII, HeIII en eV
g = gHmin, gHI, gHII, gHeI, gHeII, gHeIII
XX = XHmin, XHI, XHII, XHeI, XHeII, XHeIII


# =====================================================================================================
# ============================ TABLAS DE POBLACIONES ==================================================
# =====================================================================================================


U = []
n_e = []
NN_rel = []
NN = [] # H-, HI, HII, HeI, HeII, HeIII, H, He en cm-3

# Para los niveles excitados
position = 1, 3, 4 #Position in NN for HI, HeI, HeIIi
n_ex = []



for ii in 0,1:
    uu = np.empty([len(Pe[ii]),len(g)])
    nee = np.empty(len(Pe[ii]))
    nnrel = np.empty([len(Pe[ii]),4])
    nn = np.empty([len(Pe[ii]),8])
    n0 = np.empty([len(T[ii]),len(position),5])

    p = Pe[ii]
    t = T[ii]
    for tau in range(len(logtau[ii])):
        uu[tau] = [U_jk(g[jk],XX[jk],t[tau]) for jk in range(len(g))]
        nee[tau] = ne(p[tau],t[tau],True)
        ll = 0
        for j in 0,1,3,4:
            nnrel[tau][ll] = Saha(nee[tau],t[tau],uu[tau][j],uu[tau][j+1],X[j])
            ll +=1
        nn[tau] = system(nH,nHe,nee[tau],nnrel[tau][0],nnrel[tau][1],nnrel[tau][2],nnrel[tau][3])

        for g1 in range(len(position)):
            for g2 in range(len(XX[position[g1]])):
                n0[tau][g1][g2] = Boltzmann(nn[tau][position[g1]],g[position[g1]][g2],uu[tau][position[g1]],XX[position[g1]][g2],t[tau])


    n_e.append(nee)
    NN_rel.append(nnrel)
    NN.append(nn)
    U.append(uu)
    n_ex.append(n0)

# ===============================================================================================
# ===============================================================================================
# ===============================================================================================

# Para calcular ahora las opacidades

from opacities import *

def lamb(Z,E,n,H):
    hc = 4.135667731e-15*3e8 #eV/s x m/s
    Rh = 1.0967758341e7 #1/m
    if H == True:
        return (n**2/(Rh*(Z**2)))*1e10 #[A]
    else:
        return (hc/E)*1e10 #[A]


# Vamos ahora a calcular por separado, empezando por H-, el cual es un atomo no hidrogenoide, por lo cual:

Z_H = 1.
Z_He = 2.

wl_Hion = lamb(Z_H,np.abs(X[0]),1.,False)

wl_HI1 = lamb(Z_H, X[1],1.,True) # [2]
wl_HI2 = lamb(Z_H, XHI[1],2.,True) # [3]
wl_HI3 = lamb(Z_H, XHI[2],3.,True)
wl_HI4 = lamb(Z_H, XHI[3],4.,True)

wl_HeI1 = lamb(Z_He, X[3], 1., False) # [1] Estas son las que piden la tabla
wl_HeI2 = lamb(Z_He, XHeI[1], 2.,False)
wl_HeI3 = lamb(Z_He, XHeI[2], 3.,False)
wl_HeI4 = lamb(Z_He, XHeI[3], 4.,False)
wl_HeI5 = lamb(Z_He, XHeI[4], 5.,False)

wl_HeII1 = lamb(Z_He, X[4], 1., True)
wl_HeII2 = lamb(Z_He, XHeII[1], 2., True)

deltawl = 1e-4

# Longitudes de onda que piden en la tabla

wl1 = [wl_HeI1-deltawl*wl_HeI1, wl_HeI1+deltawl*wl_HeI1]
wl2 = [wl_HI1-deltawl*wl_HI1, wl_HI1+deltawl*wl_HI1]
wl3 = [wl_HI2-deltawl*wl_HI2, wl_HI2+deltawl*wl_HI2]

WL = [wl1, wl2, wl3]

# =======================================================================================================
# =============================== TABLAS DE OPACIDAD ====================================================
# =======================================================================================================

kfHion = []
kfHI = []
kfHeI = []
kfHeII = []
kbfHion = []
kbfHIn1 = []
kbfHIn2 =[]
kbfHIn3 = []
kbfHIn4 = []
kbfHeIn1 = []
kbfHeIn2 = []
kbfHeIn3 = []
kbfHeIn4 = []
kbfHeIn5 = []
kbfHeIIn1 = []
kbfHeIIn2 = []

op_e = []
ktotal = []

for ii in 0,1:
    oHion = np.array([])
    oHI = np.array([])
    oHeI = np.array([])
    oHeII = np.array([])
    bHion = np.array([])
    bHI1 = np.array([])
    bHI2 = np.array([])
    bHI3 = np.array([])
    bHI4 = np.array([])
    bHeI1 = np.array([])
    bHeI2 = np.array([])
    bHeI3 = np.array([])
    bHeI4 = np.array([])
    bHeI5 = np.array([])
    bHeII1 = np.array([])
    bHeII2 = np.array([])
    
    oe = np.array([])
    for ww in WL:
        for kk in 0,1:
            oHion = np.append(oHion,Hion_op_ff(Pe[ii][40],ww[kk],T[ii][40],NN[ii][40][1]))
            oHI = np.append(oHI,HI_op_ff(ww[kk],NN[ii][40][2],T[ii][40],n_e[ii][40]))
            oHeI = np.append(oHeI,HeI_op_ff(ww[kk],NN[ii][40][4],T[ii][40],n_e[ii][40]))
            oHeII = np.append(oHeII,HeII_op_ff(ww[kk],NN[ii][40][5],T[ii][40],n_e[ii][40]))
            bHion = np.append(bHion,Hion_op_bf(Pe[ii][40],ww[kk],wl_Hion, T[ii][40],NN[ii][40][0]))
            bHI1 = np.append(bHI1,HI_op_bf(ww[kk],wl_HI1,1.,n_ex[ii][40][0][0],T[ii][40]))
            bHI2 = np.append(bHI2,HI_op_bf(ww[kk],wl_HI2,2.,n_ex[ii][40][0][1],T[ii][40]))
            bHI3 = np.append(bHI3,HI_op_bf(ww[kk],wl_HI3,3.,n_ex[ii][40][0][2],T[ii][40]))
            bHI4 = np.append(bHI4,HI_op_bf(ww[kk],wl_HI4,4.,n_ex[ii][40][0][3],T[ii][40]))
            bHeI1 = np.append(bHeI1,HeI_op_bf(ww[kk],wl_HeI1,1.,n_ex[ii][40][1][0],T[ii][40]))
            bHeI2 = np.append(bHeI2,HeI_op_bf(ww[kk],wl_HeI2,2.,n_ex[ii][40][1][1],T[ii][40]))
            bHeI3 = np.append(bHeI3,HeI_op_bf(ww[kk],wl_HeI3,3.,n_ex[ii][40][1][2],T[ii][40]))
            bHeI4 = np.append(bHeI4,HeI_op_bf(ww[kk],wl_HeI4,4.,n_ex[ii][40][1][3],T[ii][40]))
            bHeI5 = np.append(bHeI5,HeI_op_bf(ww[kk],wl_HeI5,5.,n_ex[ii][40][1][4],T[ii][40]))
            bHeII1 = np.append(bHeII1,HeII_op_bf(ww[kk],wl_HeII1,1.,n_ex[ii][40][2][0],T[ii][40]))
            bHeII2 = np.append(bHeII2,HeII_op_bf(ww[kk],wl_HeII2,2.,n_ex[ii][40][2][1],T[ii][40]))
            oe = np.append(oe,e_op(n_e[ii][40]))
    kfHion.append(oHion)
    kfHI.append(oHI)
    kfHeI.append(oHeI)
    kfHeII.append(oHeII)
    kbfHion.append(bHion)
    kbfHIn1.append(bHI1)
    kbfHIn2.append(bHI2)
    kbfHIn3.append(bHI3)
    kbfHIn4.append(bHI4)
    kbfHeIn1.append(bHeI1)
    kbfHeIn2.append(bHeI2)
    kbfHeIn3.append(bHeI3)
    kbfHeIn4.append(bHeI4)
    kbfHeIn5.append(bHeI5)
    kbfHeIIn1.append(bHeII1)
    kbfHeIIn2.append(bHeII2)
    op_e.append(oe)

    ktotal.append(oHion+oHI+oHeI+oHeII+bHion+bHI1+bHI2+bHI3+bHI4+bHeI1+bHeI2+bHeI3+bHeI4+bHeI5+bHeII1+bHeII2+oe)
# =======================================================================================================
# =======================================================================================================
# =======================================================================================================



# Realizamos ahora los plot para el array de la longitud de onda

wl = np.linspace(500,20000,10000) #Longitud de onda en A

Hionff = []
HIff = []
HeIff = []
HeIIff = []

ff_total = []

Hionbf = []
HIbf1 = []
HIbf2 =[]
HIbf3 = []
HIbf4 = []
HeIbf1 = []
HeIbf2 = []
HeIbf3 = []
HeIbf4 = []
HeIbf5 = []
HeIIbf1 = []
HeIIbf2 = []

bf_total = []

elec_op = []

k_total = []

for ii in 0,1:
    Hionff.append(Hion_op_ff(Pe[ii][40],wl,T[ii][40],NN[ii][40][1]))
    HIff.append(HI_op_ff(wl,NN[ii][40][2],T[ii][40],n_e[ii][40]))
    HeIff.append(HeI_op_ff(wl,NN[ii][40][4],T[ii][40],n_e[ii][40]))
    HeIIff.append(HeII_op_ff(wl,NN[ii][40][5],T[ii][40],n_e[ii][40]))

    ff_total.append(Hionff[ii]+HIff[ii]+HeIff[ii]+HeIIff[ii])
    bHion = np.array([])
    bHI1 = np.array([])
    bHI2 = np.array([])
    bHI3 = np.array([])
    bHI4 = np.array([])
    bHeI1 = np.array([])
    bHeI2 = np.array([])
    bHeI3 = np.array([])
    bHeI4 = np.array([])
    bHeI5 = np.array([])
    bHeII1 = np.array([])
    bHeII2 = np.array([])
    oe = np.array([])

    for ll in wl:
        bHion = np.append(bHion,Hion_op_bf(Pe[ii][40],ll,wl_Hion, T[ii][40],NN[ii][40][0]))
        bHI1 = np.append(bHI1,HI_op_bf(ll,wl_HI1,1.,n_ex[ii][40][0][0],T[ii][40]))
        bHI2 = np.append(bHI2,HI_op_bf(ll,wl_HI2,2.,n_ex[ii][40][0][1],T[ii][40]))
        bHI3 = np.append(bHI3,HI_op_bf(ll,wl_HI3,3.,n_ex[ii][40][0][2],T[ii][40]))
        bHI4 = np.append(bHI4,HI_op_bf(ll,wl_HI4,4.,n_ex[ii][40][0][3],T[ii][40]))
        bHeI1 = np.append(bHeI1,HeI_op_bf(ll,wl_HeI1,1.,n_ex[ii][40][1][0],T[ii][40]))
        bHeI2 = np.append(bHeI2,HeI_op_bf(ll,wl_HeI2,2.,n_ex[ii][40][1][1],T[ii][40]))
        bHeI3 = np.append(bHeI3,HeI_op_bf(ll,wl_HeI3,3.,n_ex[ii][40][1][2],T[ii][40]))
        bHeI4 = np.append(bHeI4,HeI_op_bf(ll,wl_HeI4,4.,n_ex[ii][40][1][3],T[ii][40]))
        bHeI5 = np.append(bHeI5,HeI_op_bf(ll,wl_HeI5,5.,n_ex[ii][40][1][4],T[ii][40]))
        bHeII1 = np.append(bHeII1,HeII_op_bf(ll,wl_HeII1,1.,n_ex[ii][40][2][0],T[ii][40]))
        bHeII2 = np.append(bHeII2,HeII_op_bf(ll,wl_HeII2,2.,n_ex[ii][40][2][1],T[ii][40]))
        oe = np.append(oe,e_op(n_e[ii][40]))

    Hionbf.append(bHion)
    HIbf1.append(bHI1)
    HIbf2.append(bHI2)
    HIbf3.append(bHI3)
    HIbf4.append(bHI4)
    HeIbf1.append(bHeI1)
    HeIbf2.append(bHeI2)
    HeIbf3.append(bHeI3)
    HeIbf4.append(bHeI4)
    HeIbf5.append(bHeI5)
    HeIIbf1.append(bHeII1)
    HeIIbf2.append(bHeII2)

    bf_total.append(Hionbf[ii]+HIbf1[ii]+HIbf2[ii]+HIbf3[ii]+HIbf4[ii]+HeIbf1[ii]+\
    HeIbf2[ii]+HeIbf3[ii]+HeIbf4[ii]+HeIbf5[ii]+HeIIbf1[ii]+HeIIbf2[ii])

    elec_op.append(oe)

    k_total.append(ff_total[ii]+bf_total[ii]+elec_op[ii])

    


# ==================================================================================
# ==================================== PLOTS =======================================
# ==================================================================================


import matplotlib.pyplot as plt

grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.2)

Tlabel = '5000', '8000'
fig0 = plt.figure()
pos = (0,0),(0,1)


for ii in 0,1:
    fig0.add_subplot(grid[pos[ii]])
    plt.yscale('log',basey=10)
    plt.xscale('log',basex=10)

    plt.plot(wl,elec_op[ii],label=r'$\kappa_e$',ls=':',color='navy')

    plt.plot(wl,Hionff[ii],label=r'$\kappa_{ff}(H^{-})$',ls='--',color='C0')
    plt.plot(wl,HIff[ii],label=r'$\kappa_{ff}(HI)$',ls='--',color='C2')
    plt.plot(wl,HeIff[ii],label=r'$\kappa_{ff}(HeI)$',ls='--',color='C5')
    plt.plot(wl,HeIIff[ii],label=r'$\kappa_{ff}(HeII)$',ls='--',color='C4')
    plt.plot(wl,ff_total[ii],label=r'$\kappa_{ff}$',ls='--',color='C3')
    plt.plot(wl,Hionbf[ii],label=r'$\kappa_{bf}(H^{-})$',ls='-.',color='C0')
    plt.plot(wl,HIbf1[ii]+HIbf2[ii]+HIbf3[ii]+HIbf4[ii],label=r'$\kappa_{bf}(HI)$',ls='-.',color='C2')
    plt.plot(wl,HeIbf1[ii]+HeIbf2[ii]+HeIbf3[ii]+HeIbf4[ii]+HeIbf5[ii],label=r'$\kappa_{bf}(HeI)$',ls='-.',color='C5')
    plt.plot(wl,HeIIbf1[ii]+HeIIbf2[ii],label=r'$\kappa_{bf}(HeII)$',ls='-.',color='C4')
    plt.plot(wl,bf_total[ii],label=r'$\kappa_{bf}$',ls='-.',color='C3')

    plt.plot(wl,k_total[ii],label=r'$\kappa_{TOTAL}$',color='k')
    plt.legend(loc='lower right',ncol = 3)
    plt.xlabel(r'$\lambda[\AA]$')
    plt.ylabel(r'$\kappa[cm^{-1}]$')
    plt.title('T = '+str(Tlabel[ii])+' K')
    plt.tight_layout()

    
fig0.add_subplot(grid[1,:])
plt.yscale('log',basey=10)
plt.xscale('log',basex=10)
plt.plot(wl,k_total[0],label='T = 5000 K',color='C8')
plt.plot(wl,k_total[1],label='T = 8000 K',color='C9')
plt.title(r'$\kappa_{TOTAL}$')
plt.xlabel(r'$\lambda[\AA]$')
plt.ylabel(r'$\kappa[cm^{-1}]$')
plt.legend()
plt.tight_layout()


# =================================================================================
# ========================== PARA EL SEGUNDO PLOT =================================
# =================================================================================


k3640 = []
ff3640 = []
bf3640 = []
Hionff3640 = []
Hionbf3640 = []
HIff3640 = []
HIbf3640 = []
HeIff3640 = []
HeIbf3640 = []
HeIIff3640 = []
HeIIbf3640 = []
e3640 = []

fig1 = plt.figure()
for ii in 0,1:
    kt = np.array([])
    ff = np.array([])
    bf = np.array([])
    Hionf = np.array([])
    Hionb = np.array([])
    HIf = np.array([])
    HIb = np.array([])
    HeIf = np.array([])
    HeIb = np.array([])
    HeIIf = np.array([])
    HeIIb = np.array([])
    e = np.array([])

    for nn in range(len(T[ii])):
        Hionf = np.append(Hionf,Hion_op_ff(Pe[ii][nn],3640,T[ii][nn],NN[ii][nn][1]))
        HIf = np.append(HIf,HI_op_ff(3640,NN[ii][nn][2],T[ii][nn],n_e[ii][nn]))
        HeIf = np.append(HeIf,HeI_op_ff(3640,NN[ii][nn][4],T[ii][nn],n_e[ii][nn]))
        HeIIf = np.append(HeIIf,HeII_op_ff(3640,NN[ii][nn][5],T[ii][nn],n_e[ii][nn]))
        Hionb = np.append(Hionb,Hion_op_bf(Pe[ii][nn],3640,wl_Hion, T[ii][nn],NN[ii][nn][0]))
        HIb = np.append(HIb,[HI_op_bf(3640,wl_HI1,1.,n_ex[ii][nn][0][0],T[ii][nn])+\
                HI_op_bf(3640,wl_HI2,2.,n_ex[ii][nn][0][1],T[ii][nn])+\
                HI_op_bf(3640,wl_HI3,3.,n_ex[ii][nn][0][2],T[ii][nn])+\
                HI_op_bf(3640,wl_HI4,4.,n_ex[ii][nn][0][3],T[ii][nn])])
        HeIb = np.append(HeIb,HeI_op_bf(3640,wl_HeI1,1.,n_ex[ii][nn][1][0],T[ii][nn])\
                +HeI_op_bf(3640,wl_HeI2,2.,n_ex[ii][nn][1][1],T[ii][nn])+\
                HeI_op_bf(3640,wl_HeI3,3.,n_ex[ii][nn][1][2],T[ii][nn])+\
                HeI_op_bf(3640,wl_HeI4,4.,n_ex[ii][nn][1][3],T[ii][nn])+\
                HeI_op_bf(3640,wl_HeI5,5.,n_ex[ii][nn][1][4],T[ii][nn]))
        HeIIb = np.append(HeIIb,HeII_op_bf(3640,wl_HeII1,1.,n_ex[ii][nn][2][0],\
                T[ii][nn])+HeII_op_bf(3640,wl_HeII2,2.,n_ex[ii][nn][2][1],T[ii][nn]))
        e = np.append(e,e_op(n_e[ii][nn]))
    ff = np.append(ff,Hionf+HIf+HeIf+HeIIf)
    bf = np.append(bf,Hionb+HIb+HeIb+HeIIb)
    kt = np.append(kt,ff+bf+e)

    e3640.append(e)
    k3640.append(kt)
    ff3640.append(ff)
    bf3640.append(bf)
    Hionff3640.append(Hionf)
    Hionbf3640.append(Hionb)
    HIff3640.append(HIf)
    HIbf3640.append(HIb)
    HeIff3640.append(HeIf)
    HeIbf3640.append(HeIb)
    HeIIff3640.append(HeIIf)
    HeIIbf3640.append(HeIIb)
    
    fig1.add_subplot(grid[pos[ii]])
    plt.yscale('log',basey=10)
    plt.xscale('log',basex=10)   
    plt.plot(T[ii],e,label=r'$\kappa_e$',ls=':',color='navy')
    plt.plot(T[ii],Hionf,label=r'$\kappa_{ff}(H^{-})$',ls='--',color='C0')
    plt.plot(T[ii],Hionb,label=r'$\kappa_{bf}(H^{-})$',ls='-.',color='C0')
    plt.plot(T[ii],HIf,label=r'$\kappa_{ff}(HI)$',ls='--',color='C2')
    plt.plot(T[ii],HIb,label=r'$\kappa_{bf}(HI)$',ls='-.',color='C2')
    plt.plot(T[ii],HeIf,label=r'$\kappa_{ff}(HeI)$',ls='--',color='C5')
    plt.plot(T[ii],HeIb,label=r'$\kappa_{bf}(HeI)$',ls='-.',color='C5')
    plt.plot(T[ii],HeIIf,label=r'$\kappa_{ff}(HeII)$',ls='--',color='C4')
    plt.plot(T[ii],HeIIb,label=r'$\kappa_{bf}(HeII)$',ls='-.',color='C4')
    plt.plot(T[ii],bf,label=r'$\kappa_{bf}$',ls='-.',color='C3')
    plt.plot(T[ii],ff,label=r'$\kappa_{ff}$',ls='--',color='C3')
    plt.plot(T[ii],kt,label=r'$\kappa_{TOTAL}$',color='k')
    plt.legend(loc='lower right',ncol = 3)
    plt.xlabel('T [K]')
    plt.ylabel(r'$\kappa[cm^{-1}]$')
    plt.title(r'$T_{model}$ = '+str(Tlabel[ii])+' K')
    plt.tight_layout()




fig1.add_subplot(grid[1,:])
plt.yscale('log',basey=10)
plt.xscale('log',basex=10)
plt.plot(T[0],k3640[0],label='$T_{model}$ = 5000 K',color='C8')
plt.plot(T[1],k3640[1],label='$T_{model}$ = 8000 K',color='C9')
plt.title(r'$\kappa_{TOTAL}$')
plt.xlabel('T [K]')
plt.ylabel(r'$\kappa[cm^{-1}]$')
plt.legend()
plt.tight_layout()

