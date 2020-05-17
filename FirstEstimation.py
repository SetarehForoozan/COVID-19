# -*- coding: utf-8 -*-
"""
Created on Fri May 15 00:40:25 2020

@author: Asus
"""


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#Differential equation

#Recovery Rate
gamma_H = 1
gamma_M = 1
gamma_L = 1

#Ebtela Rate!
beta_H = 1.8
beta_M = 1.5
beta_L = 1

#alpha explains how much three groups are isolated 
alpha = 1.2

#how the M group lose its saving and moves to H group
eta = 1/14

#how they enter the L group!
epsilon = 0.001

def pandemic(k,t):
    
    S_H = k[0]
    S_M = k[1]
    S_L = k[2]
    I_H = k[3]
    I_M = k[4]
    I_L = k[5]
    
    dS_Hdt = -beta_H * (alpha * I_H + I_M + I_L) * (S_H + S_M + S_L)
    dS_Mdt = -beta_M * (alpha * I_M + I_H + I_L) * (S_H + S_M + S_L)
    dS_Ldt = -beta_L * (alpha * I_L + I_H + I_M) * (S_H + S_M + S_L)
    
    dI_Hdt = beta_H * (alpha * I_H + I_M + I_L) * (S_H + S_M + S_L) - gamma_H * I_H + eta * I_M
    dI_Mdt = beta_M * (alpha * I_M + I_H + I_L) * (S_H + S_M + S_L) - gamma_M * I_M - eta * I_M -  epsilon * I_M
    dI_Ldt = beta_L * (alpha * I_L + I_H + I_M) * (S_H + S_M + S_L) - gamma_L * I_L + epsilon * I_M
    
    return [dS_Hdt, dS_Mdt, dS_Ldt, dI_Hdt, dI_Mdt, dI_Ldt]

#Initial conditions and t range
stp = 0.1
t = np.arange(0,15, stp)
k0=[0.6, 0.3, 0.1, 0, 0, 0.00001]

#Solve
k = odeint(pandemic,k0,t)

"""
s = np.size(k[:, 1])
t = []
for i in range(s):
    if i == (s-1):
        break
    if ((k[i, 1] > (k[i-1, 1])) & (k[i, 1] > (k[i+1, 1]))):
        t.append(i * stp)

T = (t[np.size(t)-1]-t[0])/(np.size(t)-1)
print(T)


"""
#Plotting
plt.plot(t,k[:,3],'r-', linewidth=2.0, alpha = 0.6)
plt.plot(t,k[:,4],'b-', linewidth=2.0, alpha = 0.6)
plt.plot(t,k[:,5],'g-', linewidth=2.0, alpha = 0.6)
#plt.plot(t,k[:,0],'-', color = 'orange', linewidth=2.0, alpha = 0.6)



#print(k[20, 0])
plt.title('$Population Dynamics$')
plt.xlabel("$t_{(days)}$")
plt.ylabel("[S, I]")
plt.legend(["$I_H$", "$I_M$", "$I_L$"])
plt.savefig("Plot.png")
plt.show()
