import numpy as np
import matplotlib.pyplot as pp


def U_jk(g,X,T):
    k=1.3806488e-16 #erg/K
    erg = 1./6.242e11 #erg/eV
    g = np.array(g)
    X = np.array(X)
    return np.sum(g*np.exp(-X*erg/(k*T)))


def ne(Pe, T, ideal):
  me=9.11e-28 #electron mass [g]
  hb=1.05e-27 #Planck's constant (h-bar) [ergxs]
  C=((3.*np.pi**2.)**(2./3.))*(hb**2.)/(5.*me) #constant (which defines de electron pressure)
  R = 8.314472e7 #[egr/(molK)]
  if ideal == False:
      ne=(Pe/C)**(3/5)
      return ne*6.022e23 #atoms/cm3
  else:
      return Pe/(R*T)*6.022e23 #atoms/cm3


def Saha(ne,T,U1,U2,X):
    k=1.3806488e-16 #erg/K
    erg = 1./6.242e11 #erg/eV
    pobrel = (2.07e-16)*ne*(U1/U2)*(T**(-3./2.))*np.exp(X*erg/(k*T))
    return pobrel


def system(nH,nHe,ne, HminusHI, HIHII, HeIHeII, HeIIHeIII):
    abund = 10**(nH-nHe)
#       0 = nHminus + nHI + nHII - nH - 0nHe
#       0 = nHeI + nHeII + nHeIII - 0nH - nHe
#       0 = abundHe - nH      
#       ne = (-1.)*nHminus + 0.*nHI + 1.*nHII + 0*nHeI + 1*nHeII + 2*nHeIII
#       0 = nHminus - nHI*HminusHI
#       0 = nHI - nHII*HIHII
#       0 = nHeI - nHeII*HeIHeII
#       0 = nHeII - nHeIII*HeIIHeIII
        # Escribimos el sistema de eq en forma de matriz
    coef = np.array([[1., 1., 1., 0., 0., 0., -1., 0.],\
                [0., 0., 0., 1., 1., 1., 0., -1.],\
                [0., 0., 0., 0., 0., 0., -1., abund],\
                [-1., 0., 1., 0., 1., 2., 0., 0.],\
                [1., -HminusHI, 0., 0., 0., 0., 0., 0.],\
                [0., 1., -HIHII, 0., 0., 0., 0., 0.],\
                [0., 0., 0., 1., -HeIHeII, 0., 0., 0.],\
                [0., 0., 0., 0., 1., -HeIIHeIII, 0., 0.]])
    sol = np.array([0.,0.,0.,ne,0.,0.,0.,0.])
    pob = np.linalg.inv(coef).dot(sol)
    return pob


def Boltzmann(N,g, U, X, T):
     # La J se puede encontrar en el NIST
    k=1.3806488e-16 #erg/K
    erg = 1./6.242e11 #erg/eV
    n = N*(g/U)*np.exp(-X*erg/(k*T))
    return n
