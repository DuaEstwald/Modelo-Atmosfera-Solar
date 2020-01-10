#|************************************************************************************************|
#|--------------------------------------- H & He opacities ---------------------------------------|
#|************************************************************************************************|
#|                            By Elena Arjona Galvez & Andres Laza Ramos                          |
#|************************************************************************************************|

#|0|=================================== 0. Module Importation ==================================|0|

import numpy as np

import matplotlib.pyplot as pp
pp.ion()
pp.close('all')

#|0|============================================================================================|0|

#**************************************************************************************************

#|1|===================================== 1. What is this? =====================================|1|

  #The following .py file is a module containing definition (i.e., functions) that allow working
  #out the opacities of all (at least, those that contribute more) atomic states of Hydrogen and
  #Helium, this is H-, HI, HII, HeI, HeII, HeIII for bound-free (bf) and free-free (ff) cases.
  #Also, as the effective section are needed to work out the opacities, this module also contains
  #definitions that allow them to be worked out.
  
  #The module is organized as follows:
  #	- 0. Module importation (mor basic python modules needed for this module to properly work)
  #	- 1. What is this? (description of the module, containing the index)
  #	- 2. General definitions (mathematicla definitions of effective section and opacity for general cases)
  #	- 3. Hydrogen (effective sections and opacities for H-, HI & HII)
  #		- 2.1. H- (negative hydrogen ion, this is,  one proton with two electros; bf, ff)
  #		- 2.2. HI (neutral hydrogen, this is, one proton and one electron; bf, ff)
  #	- 4. Helium (effective sections and opacities for HeI, HeII & HeIII)
  #		- 3.1. HeI (neutral helium, this is, a helium nucleus with 2 electrons; ff)
  #		- 3.2. HeII (first level excitef helieum, this is, a helium nucleus with one electron; bf, ff)
  #	- 5. Electron scattering

#|1|============================================================================================|1|

#**************************************************************************************************

#|2|================================== 2. General definitions ==================================|2|

  #Mathematically, the effective section has a general expression for atoms depending on their temperature T
  #(i.e. the ffective temperature fo the atmosphere), the absorbed wavelenght and the their atomic
  #level n. Some atoms or configurations might have different expression, that will be specify 
  #when needed.
  
  #For the definitions below it must be noted that we are considering LTE.
  
  #Effective section
  #Bound-free
def es_bf(wl,Z,n,g=1): #Input must be at least wavelength (angstroms), atomic number and atomic level; the Gauntlet factor can be choosen to be 1 or its real definition by g='1' or g='real', respectively. Default is 1.
  if g==1:
    c=3e8 #light speed [m/s]

    fr=c/(wl*1e-10) #frequencies associated to the given wavelengths [s-1]
    
    return (2.815e29)*(Z**4)/((n**5)*(fr**3)) #effective section for this case [cm2]
    
  
  elif g=='real':
    R=1.0968e-3 #Rydberg's constant [angstroms-1]
    c=3e8 #light speed [m/s]

    fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]
  
    g=1-0.3456*((wl*R/n**2)-0.5)/((wl*R)**(1/3)) #Gauntlet factor's true definition 
    return (2.815e29)*(Z**4)*g/((n**5)*(fr**3)) #effective section for this case [cm2]
    

  elif g!=1 or g!='real':
    return print("g must be set to 1 or 'real'")
     
  
  #Free-free
def es_ff(wl,Z,T,g=1): #Input must be at least wavelength (angstroms), atomic number and temperature (IS units); the Gauntlet factor can be choosen to be 1 or its real definition by g='1' or g='real', respectively. Default is 1.
  if g==1:
   c=3.e8 #light speed (IS units)
    
   fr=c/(wl*1e-10) #frecuency associated to the given wavelenght [s-1]
   
   return (3.7e8)*(Z**2.)/(np.sqrt(T)*(fr**3.)) #effective section for this case
   

  elif g=='real':
    R=1.0968e-3 #Rydberg's constant [angstroms]
    c=3e8 #light speed [m/s]
    k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
    h=6.626e-34 #Planck's constante [kgm2/s]
  
    fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]
    
    g=1+0.3456*((k*T/fr*h)+0.5)/(wl*R)**(1/3)#Gauntlet factor's true definition 
    return (3.7e8)*(Z**2)*g/(np.sqrt(T)*(fr**3)) #effective section for this case
    
  
  elif g!=1 or g!='real':
    return print("g must be set to 1 or 'real'")

  
  #Opacity
  #Bound-free
def k_bf(wl,wl0,Z,n,N,T,g=1.): #imput must be wavelength (angstroms), atomic number, atomic level, population, temperature (K)
    if wl0 >= wl:
        h=6.626e-34 #Planck's constante [kgm2/s]
        k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
        c=3e8 #Light speed [m/s]
        fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]
        es=es_bf(wl,Z,n,g)
        return es*N*(1.-np.exp(-h*fr/(k*T)))  #matrix of len(wl) rows and len(n) columns, where line i contains the opacity for the i wavelenght in the wavelength array, and the element j of that row contains the opacity for the level j in the n array.
    else:
        return 0.
  #Free-free
def k_ff(wl,Z,T,N,ne,g=1): #imput must be the effective section (cm2), wavelength (angstroms), temperature (IS units), atomic level and number of electrons. Both wl and n can be arrays
  h=6.626e-34 #Planck's constante [kgm2/s]
  k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
  c=3e8 #Light speed [m/s]
 
  fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]
  
  es=es_ff(wl,Z,T,g)
  return es*N*ne*(1.-np.exp(-h*fr/(k*T)))

  
  

#|2|============================================================================================|2|

#**************************************************************************************************

#|3|======================================== 3. Hydrogen =======================================|3|

#|3.1|---------------------------------------- 3.1. H- ---------------------------------------|3.1|

  #Ref.: Gray (2005), page 155-157
  #Bound-Free
def Hion_es_bf(wl): #imput must be a wavelength (or array of them) in angstroms
  #Constants needed for the expression
  a0=1.99654
  a1=-1.18267e-5
  a2=2.64243e-6
  a3=-4.40524e-10
  a4=3.23992e-14
  a5=-1.39568e-18
  a6=2.78701e-23
  
  return (a0+a1*wl+a2*wl**2+a3*wl**3+a4*wl**4+a5*wl**5+a6*wl**6)*1e-18 #cm2/H-
  

  #Free-Free
def Hion_es_ff(wl,T): #input must be the star's effective wavelength and temperature, in angstroms and Kelvin
  theta=5040/T
  #coefficients and function needed for the definition
  f0=-2.2763-1.6850*np.log10(wl)+0.76661*(np.log10(wl))**2-0.053346*(np.log10(wl))**3
  f1=15.2827-9.2846*np.log10(wl)+1.99381*(np.log10(wl))**2-0.142631*(np.log10(wl))**3
  f2=-197.789+190.266*np.log10(wl)-67.9975*(np.log10(wl))**2+10.6913*(np.log10(wl))**3-0.625151*(np.log10(wl))**4
  F=-26+f0+f1*np.log10(theta)+f2*(np.log10(theta))**2
  
  es=10**F
  
  return es # in unit of cm2 per neutral hydrogen atom
  
  #Opacity
  #Bound-Free
def Hion_op_bf(Pe,wl,wl0, T,N): #input must be the electron pressure, wavelength and temperature, in angstroms and Kelvin
    if wl0>=wl:

        h=6.626e-34 #Planck's constante [kgm2/s]
        k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
        c=3e8 #Light speed [m/s]
        fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]

        es=Hion_es_bf(wl)
        theta=5040/T
  
#        k=4.158e-10*es*Pe*(theta**(5./2.))*10.**(0.754*theta)
  
        return es*N*(1.-np.exp(-h*fr/(k*T))) # in unit of cm2 per neutral hydrogen atom
    else:
        return 0.
  #Free-free
def Hion_op_ff(Pe,wl,T,N):#input must be the electron pressure, wavelength, temperature and H- population,and N the population of HI, [erg/cm3],[angstroms],[K]
  es=Hion_es_ff(wl,T)
  
  k=Pe*es*N 
  return k
#|3.1|----------------------------------------------------------------------------------------|3.1|

#|3.2|---------------------------------------- 3.2. HI ---------------------------------------|3.2|

  #Hydrogenoid atoms use the effective section and opacity explained in Section 2

  #Effective section
  #Bound-free
def HI_es_bf(wl,n,g=1.): #input must be wavelength and atomic level, both in IS units
  Z=1.
  es=es_bf(wl,Z,n,g)
  
  return es

  #Free-free
def HI_es_ff(wl,T,g=1.):#input must be wavelength and temperature, both in IS units
  Z=1.
  es=es_ff(wl,Z,T,g)
  
  return es

  #Opacity
  #Bound-free
def HI_op_bf(wl,wl0,n,N,T,g=1.): #input must be wavelength, atomic level, HI population and temperature, all in IS units
  Z=1.
  k=k_bf(wl,wl0,Z,n,N,T,g)
  
  return k

def HeII_op_bf(wl,wl0,n,N,T,g=1.): #input must be wavelength, atomic level, HeII population and Temperature, all in IS units
  Z=2.
  k=k_bf(wl,wl0,Z,n,N,T,g)

  #Free-free
def HI_op_ff(wl,N,T,ne,g=1.): #input must be wavelength,  atomic level, HI population, temperature and electron density, all in IS units
  Z=1.
  k=k_ff(wl,Z,T,N,ne,g)
  
  return k

#|3.2|----------------------------------------------------------------------------------------|3.2|

#|3.3|--------------------------------------- 3.3. HII ---------------------------------------|3.3|

  #Actually, this one is not needed, asi que pal lobby

#|3.3|----------------------------------------------------------------------------------------|3.3|

#|3|============================================================================================|2|

#**************************************************************************************************

#|4|========================================= 4. Helium ========================================|4|

#|4.1|--------------------------------------- 4.1. HeI ---------------------------------------|4.1|

  #For n>=2, due to the similarity of the energy levels of H and He, the effective sections can
  #be written as proportional to those of HI
  
  #Effective sections
  #Bound-free
def HeI_es_bf(wl,T,n,g=1.): #input must be wavelength, temperature and atomic level, all in IS units
  
    if n >= 2.:
        keV=8.6173324e-5 #Boltzmann's constant [eV/K]
        esbf=HI_es_bf(wl,n,g)
        return 4*esbf*np.exp(-10.92/(keV*T))
    else:
        c = 3e8
        fr=c/(wl*1e-10)
        return 2.951209e14*fr**(-2)
  
  
  #Free-free
def HeI_es_ff(wl,T,g=1.): #input must be wavelength, temperature and atomic level, all in IS units
    keV=8.6173324e-5 #Boltzmann's constant [eV/K]
    es=HI_es_ff(wl,T,g)
    return 4*es*np.exp(-10.92/(keV*T))
  

  #Opacity
def HeI_op_bf(wl,wl0,n,N,T,g=1.):
    if wl0 >= wl:
        h=6.626e-34 #Planck's constante [kgm2/s]
        k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
        c=3e8 #Light speed [m/s]
        fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]
        es=HeI_es_bf(wl,T,n,g)
        return es*N*(1.-np.exp(-h*fr/(k*T)))  #matrix of len(wl) rows and len(n) columns, where line i contains the opacity for the i wavelenght in the wavelength array, and the element j of that row contains the opacity for the level j in the n array.
    else:
        return 0.

def HeI_op_ff(wl,N,T,ne,g=1.):
    Z = 1.
    h=6.626e-34 #Planck's constante [kgm2/s]
    k=1.3864852e-23 #Boltzmann's constant [m2kg/s2K]
    c=3e8 #Light speed [m/s]

    fr=c/(wl*1e-10) #frequency associated to the given wavelenght [s-1]

    es=HeI_es_ff(wl,T,g)
    return es*N*ne*(1.-np.exp(-h*fr/(k*T)))


#|4.1|----------------------------------------------------------------------------------------|4.1|

#|4.2|--------------------------------------- 4.2. HeII --------------------------------------|4.2|

  #Hydrogenoid atoms use the effective section and opacity explained in Section 2
  
  #Effective section
  #Bound-free
def HeII_es_bf(wl,n,g=1.): #input must be wavelength and atomic level, both in IS units
  Z=2.
  es=es_bf(wl,Z,n,g)
  
  return es

  #Free-free
def HeII_es_ff(wl,T,g=1.): #input must be wavelength and temperature, both in IS units
  Z=2.
  es=es_ff(wl,Z,T,g)
  
  return es
  
  #Opacity
  #Bound-free
def HeII_op_bf(wl,wl0,n,N,T,g=1.): #input must be wavelength, atomic level, HeII population and Temperature, all in IS units
  Z=2.
  k=k_bf(wl,wl0,Z,n,N,T,g)
  
  return k

  #Free-free
def HeII_op_ff(wl,N,T,ne,g=1.): #input must be wavelength, atomic level, HeII population, temperature and electtron density, all in IS units
  Z=2.
  k=k_ff(wl,Z,T,N,ne,g)
  
  return k

#|4.2|----------------------------------------------------------------------------------------|4.2|

#|4|============================================================================================|4|

#**************************************************************************************************

#|5|================================== 5. Electron scattering ==================================|5|

  #For electron scattering, we will only consider the non-relativistic Thomsom scattering, which
  #considers a constant effective section

def e_op(ne): #input must be electorn density
  es=6.25e-25 #electron effective section
  
  k=es*ne
  
  return k

#|5|============================================================================================|5|
