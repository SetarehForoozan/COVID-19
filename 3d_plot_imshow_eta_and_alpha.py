#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 16:48:43 2020

@author: Setare
"""

#fig configuration

from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 12


def pandemic(k,t):
    
    S_H = k[0]
    S_M = k[1]
    S_L = k[2]
    I_H = k[3]
    I_M = k[4]
    I_L = k[5]  
    R_H = k[6]
    R_M = k[7]
    R_L = k[8]
    
    dS_Hdt = -ps['beta_H'] * (I_H + ps['alpha'] * I_M + ps['alpha'] * I_L) * S_H +  ps['eta_MH'] * S_M
    dS_Mdt = -ps['beta_M'] * (ps['alpha'] * I_H + I_M + ps['alpha'] * I_L) * S_M -  ps['eta_MH'] * S_M +  ps['eta_LM'] * S_L
    dS_Ldt = -ps['beta_L'] * (ps['alpha'] * I_H + ps['alpha'] * I_M + I_L) * S_L -  ps['eta_LM'] * S_L
    
    dI_Hdt =  ps['beta_H'] * (I_H + ps['alpha'] * I_M + ps['alpha'] * I_L) * S_H - ps['gamma_H'] * I_H 
    dI_Mdt =  ps['beta_M'] * (ps['alpha'] * I_H + I_M + ps['alpha'] * I_L) * S_M - ps['gamma_M'] * I_M 
    dI_Ldt =  ps['beta_L'] * (ps['alpha'] * I_H + ps['alpha'] * I_M + I_L) * S_L - ps['gamma_L'] * I_L 
    
    dR_Hdt =  ps['gamma_H'] * I_H 
    dR_Mdt =  ps['gamma_M'] * I_M 
    dR_Ldt =  ps['gamma_L'] * I_L 
    
    
    return [dS_Hdt, dS_Mdt, dS_Ldt, dI_Hdt, dI_Mdt, dI_Ldt, dR_Hdt, dR_Mdt, dR_Ldt]


def RDS(x, y, z): #root diff squared
    return np.sqrt((y - x) ** 2 + (z - y) ** 2)



e = 0.04
b = 2
g = 0.3
a = 0.2

pop_high = 2.57
pop_med = 6.55
pop_low = 1.16
pop_tot =  pop_high + pop_med + pop_low

ps = {'beta_H'  : b ,  'beta_M' : 38/60 * b ,'beta_L' : 19/60 * b,       #Infection rate
      'gamma_H' : g, 'gamma_M' : g,'gamma_L' : g, # death and recovery rate
      'alpha'   : a,                                            #isolation parameter
      'eta_MH'  : e, 'eta_LM' : e}                                #interchage rate


alphas = np.arange(0,1,0.01)
etas =  np.arange(0,1,0.01)

stp = 0.1
t = np.arange(0,100, stp)
time_array = np.array([t,t,t]).T

rds = np.zeros((len(alphas), len(etas)))

k0 = [pop_high /pop_tot, pop_med/pop_tot, pop_low/pop_tot, 0.001, 0.001, 0.001, 0, 0, 0]   #[high, med, low] 

j = 0
for e in etas:
    ps['eta_MH'] = e
    ps['eta_LM'] = e 
    i = 0
    for a in alphas:
        ps['alpha'] = a
        k = odeint(pandemic,k0,t)
        S = k[:,:3]
        I = k[:,3:6]
        R = k[:,6:]
        rds[i,j] = RDS(R[-1,0]/ S[0,0] , R[-1,1]/ S[0,1], R[-1,2]/ S[0,2])
        rds[i,j] = RDS(R[-1,0] , R[-1,1], R[-1,2])
        i += 1
    j += 1



fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(rds)
ax.set_xlabel('alpha')
ax.set_ylabel('eta')

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(alphas, etas)  # `plot_surface` expects `x` and `y` data to be 2D
#uncomment for 3d plot
from matplotlib import cm
surf = ax.plot_surface(X, Y, rds, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

fig.colorbar(surf, shrink=0.5, aspect=5)
'''