#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:34:00 2019

@author: mariandm
"""
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt


def alone(y0,t,pars):
    
    N1=y0[0]
    R=y0[1]
    
    dN1dt=pars['gamma1']*pars['e1']*N1*R-pars['m1']*N1
    dRdt=pars['w']*pars['R0']-pars['w']*R-pars['gamma1']*N1*R
    
    dydt=[dN1dt,dRdt]
    
    return dydt

pars={}
pars['e1']=0.5
pars['w']=0.1
pars['R0']=100
gammas=[0.1,2,10,4,2,1,0.3]
ms=[1,4,100,10,5,10,1]
N1_0=10
R_0=50

gammas=[0.1,2,10,4,2,1,0.3]
ms=[1,4,100,10,5,10,1]

t=np.linspace(0,20)

y0=[N1_0,R_0]

finaldent=[]
for i in range(0,len(gammas)):
    pars['gamma1']=gammas[i]
    pars['m1']=ms[i]
    if i==2:
        tmp=np.zeros(len(t))
        tmp[0]=10
        plt.plot(t,tmp)
        finaldent.append(0)
    else:
        y=integrate.odeint(alone,y0,t,args=(pars,))
        finaldent.append(y[-1,0])
        plt.plot(t,y[:,0])
    
plt.plot(t,t*0, linestyle="--", color="r")
#plt.title("Bacteria " + str(i+1))
plt.xlabel("Time", fontsize=15)
plt.ylabel("Density", fontsize=15)
plt.ylim([-1,27])
plt.legend(["B1","B2","B3","B4","B5","B6","B7"])
plt.savefig("P1Figure1.pdf")
plt.show()


finaldent
