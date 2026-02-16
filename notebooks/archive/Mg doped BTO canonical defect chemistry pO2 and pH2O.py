# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 12:48:34 2021

@author: jrmcg
"""

# Translation of Ryu's code into python

# This code is for solving the charge neutrality equation in doped BaTiO3
# polycrystalline ceramics. Note that the results can be a bit changed by
# the equilibirum constants. In the case of single valent acceptor-doped
# BaTiO3, all acceptor was assumed to be ionzied. Thus, their 
# concentration is not a function of oxygen partial pressure. 
# However, in some case, they are considered to be partially ionized.
# Please look at the paper published in JAP 108 064101 (2010) by Clive
# Randall, S.H. Yoon, and K. H. Hur.

#%%

import numpy as np
import pandas as pd

k = 8.61733e-5; # boltzmann constant in eV/K
k_j = 1.3806e-23 #boltzmann constant in J/K
R = 8.314 #J/Mol*K
e = 1.6E-19;
TK = 1000+273; # equilibrium temp in K
TC = TK-273; # equilibrium temp in C
TQK = 450+273;# quenched temp in K
TQC = TQK-273; # quenched temp in C
#a = np.logspace(-20,5,50); # oxygen partial pressure range
a = 0.2
h = np.logspace(-20,4,50) #water vapor partial pressure

#names = ["Cacc","Cacc_1","Cyt_1","CVO2","CVO2_1","CVTi4","CVTi4_1","CVBa2","CVBa2_1","n","nq","p","pq"]
#df = pd.DataFrame(columns = names) #create pandas df for data storage

#%% run if doing oxygen po2 dependent
M = np.zeros(len(a)); # definding vectors to fill in with values later
Cacc = np.zeros(len(a))
Cacc_1 = np.zeros(len(a))
Cyt_1 = np.zeros(len(a))
CVO2 = np.zeros(len(a))
CVO2_1 = np.zeros(len(a))
CVTi4 = np.zeros(len(a))
CVTi4_1 = np.zeros(len(a))
CVBa2 = np.zeros(len(a))
CVBa2_1 = np.zeros(len(a))
n = np.zeros(len(a))
nq = np.zeros(len(a))
p = np.zeros(len(a))
pq = np.zeros(len(a))
#%%
M = np.zeros(len(h)); # definding vectors to fill in with values later
Cacc = np.zeros(len(h))
Cacc_1 = np.zeros(len(h))
Cyt_1 = np.zeros(len(h))
CVO2 = np.zeros(len(h))
CVO2_1 = np.zeros(len(h))
CVTi4 = np.zeros(len(h))
CVTi4_1 = np.zeros(len(h))
CVBa2 = np.zeros(len(h))
CVBa2_1 = np.zeros(len(h))
n = np.zeros(len(h))
nq = np.zeros(len(h))
p = np.zeros(len(h))
pq = np.zeros(len(h))
#%%

un = 8080*(TQK**(-3/2))*np.exp(-0.021/(k*TQK)) # mobility of electron
up = un/2 # monbility of hole
uVO = (52/TQK)*np.exp(-0.84/(k*TQK)) # mobility of oxygen vacancy

Cmgt_1 = 1*7.65E+19 # single-vaelnt Mg concentration = 0.5mol%
Cacc_1 = 5E18 # unintentionally doped single-valent acceptor concentration
#C_Hydrogen_1 = 1.64E18 #unintentionally doped single-valent donor concentration - hydrogen from humidity testing
#Cyt_1 =0 #1*7.65E+19  # single-valent Y concentration = 0.5mol%

KR = 2.40E74*np.exp(-5.88/(k*TK)) # K for the reduction reaction in cm-9atm^1/2 published in J. Nowotny and M. Rekas, “Defect structure, electrical properties and transport in Barium Titanate. VI. General defect model,?Ceram. Int., 20 [4] 257?63 (1994).
#KR = 1.06E71*np.exp(-5.69/(k*TK)) #undoped BaTiO3 
Ki = 6.8E44*np.exp(-2.90/(k*TK)) # K for the intrisic equailibirum in cm-6 published in C.-R. SONG, H.-I. YOO, and J.-Y. KIM, “Mn-doped BaTiO3: Electrical Transport Properties in Equilibrium State,?J. Electroceramics, 1 [1] 27?9 (1997).
#Ki = 2.28E42*np.exp(-2.15/(k*TK))
#KS = (10**116)*np.exp(-10.5/((k*TK))) # K for the Schottky reaction in cm-15 from S. H. Yoon, K. H. Lee, and H. Kim, J. Am. Ceram. Soc. 83,2463 Íì2000? #for donor doped BTO
KS = (3.4*10**105)*np.exp(-2.795/(k*TK)) #acceptor and nominally undoped BaTiO3
#Khydr_mol = 8.1*10**-9 #WASER, R. (1988). Solubility of Hydrogen Defects in Doped and Undoped BaTiO3. Journal of the American Ceramic Society, 71(1), 58–63. https://doi.org/10.1111/j.1151-2916.1988.tb05760.x
delS = -97.5
delH = -21700
K_hydr_pa = (np.exp(delS/R)/np.exp(delH/(R*TK))) #10**5pa-1 #measured at 1atm ph2O #Kreuer, K. D. (2003). PROTON-CONDUCTING OXIDES. Annu. Rev. Mater. Res, 33, 333–359. https://doi.org/10.1146/annurev.matsci.33.022802.091825
K_hydr = K_hydr_pa
K_hydr_pa_q = (np.exp(delS/R)/np.exp(delH/(R*TQK)))
K_hydr_q = K_hydr_pa_q/1E5
#%%

from scipy.optimize import brentq
from scipy.optimize import fsolve

#for changing partial oxygen pressure
for i in range(len(a)): # This part is for obtaining a full equilibirum defect structure at TK(TC)
    y = a[i] #oxygen pressure atm
    #x = Symbol('x', real = True) # a number of electron, n, which is defined as a variable, x
    #CVO2_1 = (KR/(x**2)*(y**(-1/2))) # a converted form of oxygen vacancies
    #CVTi4_1 = (KS/(KR**3))**(1/2)*x**3*(y**(3/4)); # a converted form of cation vacancies
    def f(x): 
        return (x + 2*Cmgt_1 + Cacc_1 + 6*(KS/(KR**3))**(1/2)*x**3*(y**(3/4)) - (Ki/x) - 2*(KR/(x**2)*(y**(-1/2))))
        #return (x + 2*Cmgt_1 + Cacc_1 - (Ki/x) - 2*(KR/(x**2)*(y**(-1/2))))
   # eqn1 = x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4_1 - (Ki/x) - 2*CVO2_1 - Cyt_1 # the total charge nuetrality equation
   # f = implemented_function('f', lambda x: eqn1)
   # fn1 = lambdify(x, f(x))
    sol1 = brentq(f,1,1E30)
    #sol1 = fsolve(f,1E15)
   # sol1 = solveset(eqn1,x,domain=S.Reals)
    n[i] = np.array(sol1).astype(np.float64) # a concentration of electron
    CVO2[i] = KR/(n[i]**2)*(y**(-1/2)); # a concentration of oxygen vacancies (Reduction reaction)
    CVTi4[i] = ((KS/(KR**3))**(1/2))*(n[i]**3)*(y**(3/4)) # a concentration of cation vacancies (Schottky reaction)
    
#%%
#for changing partial water vapor pressure
from scipy.optimize import brentq
Cmgt = np.zeros(len(h))

for i in range(len(h)): # This part is for obtaining a full equilibirum defect structure at TK(TC)
    y = a #oxygen pressure atm
    pH2O = h[i]
    #x = Symbol('x', real = True) # a number of electron, n, which is defined as a variable, x
    #CVO2_1 = (KR/(x**2)*(y**(-1/2))) # a converted form of oxygen vacancies
    #CVTi4_1 = (KS/(KR**3))**(1/2)*x**3*(y**(3/4)); # a converted form of cation vacancies
    def f(x): 
        #return (x + 2*Cmgt_1 + Cacc_1 + 6*(KS/(KR**3))**(1/2)*x**3*(y**(3/4)) - (Ki/x) - 2*(KR/(x**2)*(y**(-1/2))))#(K_hydr*(KR/(x**2)*(y**(-1/2)))*pH2O)**1/2)
        return (x + 2*Cmgt_1 + Cacc_1 + 6*((KS/(KR**3))**(1/2)*x**3*(y**(3/4))) - (Ki/x) - 2*(KR/(x**2)*(y**(-1/2))))
   # eqn1 = x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4_1 - (Ki/x) - 2*CVO2_1 - Cyt_1 # the total charge nuetrality equation
   # f = implemented_function('f', lambda x: eqn1)
   # fn1 = lambdify(x, f(x))
    sol1 = brentq(f,1E-20,1E50)
   # sol1 = solveset(eqn1,x,domain=S.Reals)
    Cmgt[i] = Cmgt_1
    n[i] = np.array(sol1).astype(np.float64) # a concentration of electron
    CVO2[i] = KR/(n[i]**2)*(y**(-1/2)) # a concentration of oxygen vacancies (Reduction reaction)
    CVTi4[i] = ((KS/(KR**3))**(1/2))*(n[i]**3)*(y**(3/4)) # a concentration of cation vacancies (Schottky reaction)
#%%
#for changing partial water vapor pressure hydration reaction

Kiq=6.8E44*np.exp(-2.90/(k*TQK)); # K for the intrisic equailibirum in cm-6 published in C.-R. SONG, H.-I. YOO, and J.-Y. KIM, “Mn-doped BaTiO3: Electrical Transport Properties in Equilibrium State,?J. Electroceramics, 1 [1] 27?9 (1997).
#Kiq = 2.28E42*np.exp(-2.15/(k*TK))

Cmgt = np.zeros(len(h))
Cyt =  np.zeros(len(h))
sigman =  np.zeros(len(h))
sigmap = np.zeros(len(h))
sigmaVO2 =  np.zeros(len(h))
sigmatotal =  np.zeros(len(h))
#C_hydrogen = np.zeros(len(a))
COHO = np.zeros(len(h))

for i in range(len(h)): # This part is for obtaining a quenched defect structure at TQK(TQC) where oxygen vacancies are immobile
    y = a
    pH2O = h[i]
    #x = Symbol('x', real = True)
    CVO2q_1 = CVO2
    CVTi4q_1 = CVTi4
    #eqn2 = Eq(x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i] + Cyt_1,0)
    def g(x): 
        #return (x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i]-(K_hydr_q*(KR/(x**2)*(y**(-1/2)))*pH2O)**1/2)
        return (x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i]-((K_hydr_q*(KR/(x**2)*(y**(-1/2)))*pH2O)**1/2))
    sol2 = brentq(g,1,1E50, maxiter=500)
    nq[i] = np.array(sol2).astype(np.float64)
    pq[i] = Kiq/nq[i]
    COHO[i] = (K_hydr_q*(KR/(nq[i]**2)*(y**(-1/2)))*pH2O)**1/2 #concentration of proton defects
    Cmgt[i] = Cmgt_1
    Cacc[i] = Cacc_1
    #C_hydrogen[i] = C_Hydrogen_1
    #Cyt[i] = Cyt_1
    sigman[i] = e*un*nq[i]
    sigmap[i] = e*up*pq[i]
    sigmaVO2[i] = e*uVO*CVO2[i]
    #sigmaH
    sigmatotal[i] = sigman[i] + sigmap[i] + sigmaVO2[i]
#%%
#for changing partial water vapor pressure hydrogenation reaction

Kiq=6.8E44*np.exp(-2.90/(k*TQK)); # K for the intrisic equailibirum in cm-6 published in C.-R. SONG, H.-I. YOO, and J.-Y. KIM, “Mn-doped BaTiO3: Electrical Transport Properties in Equilibrium State,?J. Electroceramics, 1 [1] 27?9 (1997).
#Kiq = 2.28E42*np.exp(-2.15/(k*TK))

Cmgt = np.zeros(len(h))
Cyt =  np.zeros(len(h))
sigman =  np.zeros(len(h))
sigmap = np.zeros(len(h))
sigmaVO2 =  np.zeros(len(h))
sigmatotal =  np.zeros(len(h))
#C_hydrogen = np.zeros(len(a))
COHO = np.zeros(len(h))

for i in range(len(h)): # This part is for obtaining a quenched defect structure at TQK(TQC) where oxygen vacancies are immobile
    y = a
    pH2O = h[i]
    #x = Symbol('x', real = True)
    CVO2q_1 = CVO2
    CVTi4q_1 = CVTi4
    #eqn2 = Eq(x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i] + Cyt_1,0)
    def g(x): 
        #return (x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i]-(K_hydr_q*(KR/(x**2)*(y**(-1/2)))*pH2O)**1/2)
        return (x + 2*Cmgt_1 + Cacc_1 + 6*CVTi4q_1[i] - (Kiq/x) - 2*CVO2q_1[i]-((K_hydr_q*((Kiq/x)**2)*pH2O*(y**(-1/2)))**1/2))
    sol2 = brentq(g,1,1E50, maxiter=500)
    nq[i] = np.array(sol2).astype(np.float64)
    pq[i] = Kiq/nq[i]
    COHO[i] = (K_hydr_q*(KR/(nq[i]**2)*(y**(-1/2)))*pH2O)**1/2 #concentration of proton defects
    Cmgt[i] = Cmgt_1
    Cacc[i] = Cacc_1
    #C_hydrogen[i] = C_Hydrogen_1
    #Cyt[i] = Cyt_1
    sigman[i] = e*un*nq[i]
    sigmap[i] = e*up*pq[i]
    sigmaVO2[i] = e*uVO*CVO2[i]
    #sigmaH
    sigmatotal[i] = sigman[i] + sigmap[i] + sigmaVO2[i]

#%%
# pO2 Quench
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(np.log10(a),np.log10(nq), marker='s',c =[255/255, 0/255, 0/255], label=r'n')
ax1.plot(np.log10(a),np.log10(pq),marker='o', c =[0/255, 0/255, 255/255], label = r'p')
ax1.plot(np.log10(a),np.log10(CVO2),marker='d',c =[0/255, 255/255, 0/255], label = r'$V_O^{2+}$')
ax1.plot(np.log10(a),np.log10(Cacc),marker='v', c =[255/255, 0/255, 255/255], label = r'$A_{Ti}^{1-}$') 
ax1.plot(np.log10(a),np.log10(Cmgt),marker='+', c =[0/255, 255/255, 255/255], label = r'$Mg_{Ti}^{2-}$')
#ax1.plot(np.log10(a),np.log10(Cyt),marker='*', c =[0/255, 0/255, 0/255], label = r'$Y_{Ba}^{1+}$' )
ax1.plot(np.log10(a),np.log10(CVTi4),marker='h', c =[0/255, 0/255, 128/255], label = r'$V_{Ti}^{4+}$')
ax1.plot(np.log10(a),np.log10(COHO),marker='s', c =[128/255, 0/255, 128/255], label = r'$H^{1+}$')
ax1.set_xlabel(r'$log_{10} \mathit{p}O_2 (atm)$')
ax1.set_ylabel(r'$log_{10} [ ] (cm^{-3})$')
ax1.legend(loc='best',fontsize='small')
#%%
# pO2 High T
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(np.log10(a),np.log10(n), marker='s',c =[255/255, 0/255, 0/255], label=r'n')
ax1.plot(np.log10(a),np.log10(Ki/n),marker='o', c =[0/255, 0/255, 255/255], label = r'p')
ax1.plot(np.log10(a),np.log10(CVO2),marker='d',c =[0/255, 255/255, 0/255], label = r'$V_O^{2+}$')
ax1.plot(np.log10(a),np.log10(Cacc),marker='v', c =[255/255, 0/255, 255/255], label = r'$A_{Ti}^{1-}$') 
ax1.plot(np.log10(a),np.log10(Cmgt),marker='+', c =[0/255, 255/255, 255/255], label = r'$Mg_{Ti}^{2-}$')
#ax1.plot(np.log10(a),np.log10(Cyt),marker='*', c =[0/255, 0/255, 0/255], label = r'$Y_{Ba}^{1+}$' )
ax1.plot(np.log10(a),np.log10(CVTi4),marker='h', c =[0/255, 0/255, 128/255], label = r'$V_{Ti}^{4+}$')
#ax1.plot(np.log10(a),np.log10(COHO),marker='s', c =[128/255, 0/255, 128/255], label = r'$H^{1+}$')
ax1.set_xlabel(r'$log_{10} \mathit{p}O_2 (atm)$')
ax1.set_ylabel(r'$log_{10} [ ] (cm^{-3})$')
ax1.legend(loc='best',fontsize='small')
#%% for hydration

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(np.log10(h),np.log10(nq), marker='s',c =[255/255, 0/255, 0/255], label=r'n')
ax1.plot(np.log10(h),np.log10(pq),marker='o', c =[0/255, 0/255, 255/255], label = r'p')
ax1.plot(np.log10(h),np.log10(CVO2q_1),marker='d',c =[0/255, 255/255, 0/255], label = r'$V_O^{2+}$')
ax1.plot(np.log10(h),np.log10(Cacc),marker='v', c =[255/255, 0/255, 255/255], label = r'$A_{Ti}^{1-}$') 
ax1.plot(np.log10(h),np.log10(Cmgt),marker='+', c =[0/255, 255/255, 255/255], label = r'$Mg_{Ti}^{2-}$')
#ax1.plot(np.log10(a),np.log10(Cyt),marker='*', c =[0/255, 0/255, 0/255], label = r'$Y_{Ba}^{1+}$' )
ax1.plot(np.log10(h),np.log10(CVTi4q_1),marker='h', c =[0/255, 0/255, 128/255], label = r'$V_{Ti}^{4+}$')
ax1.plot(np.log10(h),np.log10(COHO),marker='s', c =[128/255, 0/255, 128/255], label = r'$OH_{O}^{1+}$')
ax1.set_xlabel(r'$log_{10} \mathit{p}H_2O (atm)$')
ax1.set_ylabel(r'$log_{10} [ ] (cm^{-3})$')
ax1.legend(loc='best',fontsize='small')

#%% hydration high T
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(np.log10(h),np.log10(n), marker='s',c =[255/255, 0/255, 0/255], label=r'n')
ax1.plot(np.log10(h),np.log10(Ki/n),marker='o', c =[0/255, 0/255, 255/255], label = r'p')
ax1.plot(np.log10(h),np.log10(CVO2),marker='d',c =[0/255, 255/255, 0/255], label = r'$V_O^{2+}$')
ax1.plot(np.log10(h),np.log10(Cacc),marker='v', c =[255/255, 0/255, 255/255], label = r'$A_{Ti}^{1-}$') 
ax1.plot(np.log10(h),np.log10(Cmgt),marker='+', c =[0/255, 255/255, 255/255], label = r'$Mg_{Ti}^{2-}$')
#ax1.plot(np.log10(a),np.log10(Cyt),marker='*', c =[0/255, 0/255, 0/255], label = r'$Y_{Ba}^{1+}$' )
ax1.plot(np.log10(h),np.log10(CVTi4),marker='h', c =[125/255, 0/255, 125/255], label = r'$V_{Ti}^{4+}$')
#ax1.plot(np.log10(h),np.log10(COHO),marker='s', c =[128/255, 0/255, 128/255], label = r'$OH_{O}^{1+}$')
ax1.set_xlabel(r'$log_{10} \mathit{p}H_2O (atm)$')
ax1.set_ylabel(r'$log_{10} [ ] (cm^{-3})$')
ax1.legend(loc='best',fontsize='small')
#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(np.log10(a),np.log10(sigmatotal), marker='x',c =[0/255, 0/255, 0/255], label=r'log$\sigma_T$')
ax1.plot(np.log10(a),np.log10(sigman),marker='o', c =[0/255, 0/255, 255/255], label = r'log$\sigma_n$')
ax1.plot(np.log10(a),np.log10(sigmap),marker='d',c =[0/255, 255/255, 0/255], label = r'log$\sigma_p$')
ax1.plot(np.log10(a),np.log10(sigmaVO2),marker='v', c =[255/255, 0/255, 255/255], label = r'log$\sigma_{V_O^{2+}}$') 
ax1.set_xlabel(r'$log_{10} \mathit{p}O_2 (atm)$')
ax1.set_ylabel(r'$log_{10} \sigma (\Omega^{-1} cm^{-1})$')
ax1.legend()

#%% for hydration
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

mg_pH2O = np.log10([9.87E-4,4.93E-2,6.908E-2,9.869E-2,0.1381,0.2566])
mg_conductivity = np.log10([5.50078E-5,1.193E-5,1.1729E-5,1.1599E-5,1.160936E-5,1.16146E-5])

fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(mg_pH2O,mg_conductivity, marker='s',c =[255/255, 20/255, 147/255], label=r'Measured log$\sigma_T$')
#ax1.plot(np.log10(9.37E-2),np.log10(1.5E-5), marker='s',c =[0/255, 0/255, 0/255])
ax1.plot(np.log10(h),np.log10(sigmatotal), marker='x',c =[0/255, 0/255, 0/255], label=r'log$\sigma_T$')
ax1.plot(np.log10(h),np.log10(sigman),marker='o', c =[255/255, 0/255, 0/255], label = r'log$\sigma_n$')
ax1.plot(np.log10(h),np.log10(sigmap),marker='d',c =[0/255, 0/255, 255/255], label = r'log$\sigma_p$')
ax1.plot(np.log10(h),np.log10(sigmaVO2),marker='v', c =[0/255, 255/255, 0/255], label = r'log$\sigma_{V_O^{2+}}$') 
ax1.set_xlabel(r'$log_{10} \mathit{p}H_2O (atm)$')
ax1.set_ylabel(r'$log_{10} \sigma (\Omega^{-1} cm^{-1})$')
ax1.legend()


#%% for hydration - total conductivity only
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


fig = plt.figure(figsize=(12, 10), dpi=80)
ax1= plt.subplot(111)
#plt.yscale('symlog',subsy=[2, 3, 4, 5, 6, 7, 8, 9])
ax1.set_title(r"$T_E = $" +str(TC)+"$^\circ C\/\/T_Q = $" +str(TQC) +"$^\circ C $")
#.set_ylim(ymin=0, ymax=1.5E6)
#ax1.set_xlim(xmin=0, xmax=) ##set x limits if needed
matplotlib.rcParams.update({'font.size': 22})
ax1.tick_params(length=6, width=2)
ax1.minorticks_on()
ax1.plot(mg_pH2O,mg_conductivity, marker='s',c =[255/255, 20/255, 147/255], label=r'Measured log$\sigma_T$')
#ax1.plot(np.log10(9.37E-2),np.log10(1.5E-5), marker='s',c =[0/255, 0/255, 0/255])
ax1.plot(np.log10(h),np.log10(sigmatotal), marker='x',c =[0/255, 0/255, 0/255], label=r'Predicted log$\sigma_T$')
#ax1.plot(np.log10(h),np.log10(sigman),marker='o', c =[255/255, 0/255, 0/255], label = r'log$\sigma_n$')
#ax1.plot(np.log10(h),np.log10(sigmap),marker='d',c =[0/255, 0/255, 255/255], label = r'log$\sigma_p$')
#ax1.plot(np.log10(h),np.log10(sigmaVO2),marker='v', c =[0/255, 255/255, 0/255], label = r'log$\sigma_{V_O^{2+}}$') 
ax1.set_xlim(xmin=-5 , xmax=5)
ax1.set_ylim(ymin=-6 , ymax=-4)
ax1.set_xlabel(r'$log_{10} \mathit{p}H_2O (atm)$')
ax1.set_ylabel(r'$log_{10} \sigma (\Omega^{-1} cm^{-1})$')
ax1.legend()

#%%

charge_sum = nq - pq - 2*CVO2 + Cacc + 2*Cmgt_1 + 6*CVTi4 - COHO
#charge_sum_HT =  n - (Ki/n) - 2*CVO2 + Cacc + 2*Cmgt + 6*CVTi4
charge_sum_HT =  n - (Ki/n) - 2*CVO2 + Cacc + 2*Cmgt_1 + 6*CVTi4
