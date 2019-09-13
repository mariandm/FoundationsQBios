#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 10:10:42 2019

@author: mariandm
"""

#Problem 2

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import *

#Parameters
pars={}
pars['a1']=5
pars['a2']=5
pars['beta']=3
pars['gamma']=3


#Initial conditions
x0=[1,1.1]
t0=0
tf=5
t=np.linspace(t0,tf)

##Writing the code

#Simulate Box1
x=integrate.odeint(ToggleBox1,x0,t,args=(pars,))



#Identifying nullclines

#####TRATAR DE PONERLE COORDENADAS A LA ORILLA DEL VECTOR

# Define the nullclines
unullcline = lambda v: pars['a1']/(1+v**pars['beta'])
vnullcline = lambda v: np.sign((pars['a2']/v-1))*np.abs(pars['a2']/v-1)**(1/pars['gamma'])


###NULLCLINES SOLITAS
# Define a vector of protein and mRNA values
vvec=np.linspace(0.04,5)
#vvec = np.arange(0,60,1)*0.1
# Plot the nullclines
plt.plot(unullcline(vvec), vvec, linewidth=3)
plt.plot(vnullcline(vvec), vvec, linewidth=3)
#for i in range(0,len(x0)):
#    x=integrate.odeint(ToggleBox1,x0[i],t,args=(pars,))
#    plt.plot(x[:,0],x[:,1],linewidth=3,color='lightgray')
#plt.plot(x[:,0],x[:,1],linewidth=3,color='b')
plt.plot(1.4,1.4,marker="o",color="black")
plt.plot(5,0.03,marker="o",color="black")
plt.plot(0.03,5,marker="o",color="black")
# Label the plot
plt.legend(['dudt = 0','dvdt = 0'],fontsize=18)
#plt.xlim(min(mvec),max(mvec))
#plt.ylim(min(pvec),max(pvec))
plt.xlabel('u concentration (nMol)',fontsize=20)
plt.ylabel('v concentration (nMol)',fontsize=20)
plt.savefig("Figure3a_nullclines.pdf")
plt.close()


####NULLCLINES CON DIFERENTES CONDICIONES INICIALES
x0=[[8,1],[1,1],[1,1.1],[1,8],[0,1],[1,0],[4,3],[7,3],[3,7]]

# Define a vector of protein and mRNA values
vvec=np.linspace(0.04,5)
#vvec = np.arange(0,60,1)*0.1
# Plot the nullclines
plt.plot(unullcline(vvec), vvec, linewidth=3, color="lightgray")
plt.plot(vnullcline(vvec), vvec, linewidth=3, color="lightgray")
for i in range(0,len(x0)):
    x=integrate.odeint(ToggleBox1,x0[i],t,args=(pars,))
    plt.plot(x[:,0],x[:,1],linewidth=3)
plt.plot(1.4,1.4,marker="o",color="black")
plt.plot(5,0.03,marker="o",color="black")
plt.plot(0.03,5,marker="o",color="black")
# Label the plot
plt.legend(['nullclines'],fontsize=18)
#plt.xlim(min(mvec),max(mvec))
#plt.ylim(min(pvec),max(pvec))
plt.xlabel('u concentration (nMol)',fontsize=20)
plt.ylabel('v concentration (nMol)',fontsize=20)
plt.savefig("Figure3b_nullexperiment.pdf")
plt.close()


t0=0
tf=10
t=np.linspace(t0,tf)

for i in range(0,len(x0)):
    x=integrate.odeint(ToggleBox1,x0[i],t,args=(pars,))
    plt.plot(t,x[:,0],linewidth=3)
plt.ylabel("u concentration (nMol)")
plt.xlabel("Times (hrs)")
plt.title("u concentration at different time points")
plt.savefig("Figure3c_uconc.pdf")
plt.close()
    
for i in range(0,len(x0)):
    x=integrate.odeint(ToggleBox1,x0[i],t,args=(pars,))
    plt.plot(t,x[:,1],linewidth=3)
plt.ylabel("v concentration (nMol)")
plt.xlabel("Times (hrs)")
plt.title("v concentration at different time points")
plt.savefig("Figure3d_vconc.pdf")
plt.close()

#a2s=np.array(range(1,21))*0.5
#for i in a2s:
#    pars['a2'] = i
#    plt.plot(unullcline(vvec), vvec, linewidth=3)
#    plt.plot(vnullcline(vvec), vvec, linewidth=3)
#    plt.show()