#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 18:17:16 2019

@author: mariandm
"""

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def model_fn(x,t,pars):
    v = x[0]
    w = x[1]
    dxdt = [v*(pars['a']-v)*(v-1)-w+pars['I'], pars['b']*v-pars['gamma']*w]
    return dxdt

pars={}
pars['a'] = 0.1
pars['b']=0.01
pars['gamma']=0.02


y0=[0, 0]
t = np.linspace(0,1000,1000)
#plt.plot(t,y[:,0],linewidth=3)
#plt.xlabel('Time',fontsize=20)
#plt.ylabel('Voltage, v',fontsize=20)
#plt.xlim([0,1000])
#plt.ylim([-0.05,0.1])
#plt.tick_params(labelsize=20)

##########################################
#############3a##########################

for currI in np.arange(0,0.255,0.01):
    pars['I']=currI
    y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)

    plt.plot(y[:,1],y[:,0],linewidth=3)

    #plt.xlim([0,1000])
    #plt.ylim([-0.05,0.1])
    #plt.tick_params(labelsize=20)
plt.xlabel('Recovery Variable, w',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.savefig("Figure3.pdf")
plt.close()

##########################################
#############3b##########################
for currI in np.arange(0,0.04,0.01):
    pars['I']=currI
    y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)

    plt.plot(y[:,1],y[:,0],linewidth=3)

plt.xlabel('Recovery Variable, w',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.savefig("Figure4a.pdf")
plt.close()

for currI in np.arange(0.04,0.22,0.01):
    pars['I']=currI
    y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)

    plt.plot(y[:,1],y[:,0],linewidth=3)
plt.xlabel('Recovery Variable, w',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.savefig("Figure4c.pdf")
plt.close()

for currI in np.arange(0.22,0.251,0.01):
    pars['I']=currI
    y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)
    plt.plot(y[:,1],y[:,0],linewidth=3)
plt.xlabel('Recovery Variable, w',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.savefig("Figure4e.pdf")
plt.close()


pars['I']=0.02
y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)
plt.plot(t,y[:,0],linewidth=3)
plt.xlabel('Time',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.xlim([-10,1000])
plt.title("I=0.02",fontsize=20)
plt.savefig("Figure4b.pdf")
plt.close()

pars['I']=0.2
y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)
plt.plot(t,y[:,0],linewidth=3)
plt.xlabel('Time',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.xlim([-10,1000])
plt.title("I=0.2",fontsize=20)
plt.savefig("Figure4d.pdf")
plt.close()

pars['I']=0.23
y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)
plt.plot(t,y[:,0],linewidth=3)
plt.xlabel('Time',fontsize=15)
plt.ylabel('Voltage, v',fontsize=15)
plt.xlim([-10,1000])
plt.title("I=0.23",fontsize=20)
plt.savefig("Figure4f.pdf")
plt.close()


##########################################
#############3c##########################
pars['I']=0.125
for e in [1,10,15,30]:
    pars['b']=0.01*e
    pars['gamma']=0.02*e
    y = integrate.odeint(model_fn,y0,t,args=(pars,),rtol=10**-6)
    plt.plot(t,y[:,0],linewidth=3)
    plt.xlabel('Time',fontsize=15)
    plt.ylabel('Voltage, v',fontsize=15)
    plt.xlim([-10,1000])
    plt.title("e = %d"%e,fontsize=20)
    plt.tick_params(labelsize=15)
    plt.savefig("Figure5a_%d.pdf"%e)
    plt.show()
    
    plt.plot(y[:,1],y[:,0],linewidth=3,color="purple")
    plt.xlabel('Response, w',fontsize=15)
    plt.ylabel('Voltage, v',fontsize=15)
    #plt.xlim([-10,1000])
    plt.title("e = %d"%e,fontsize=20)
    plt.tick_params(labelsize=15)
    plt.savefig("Figure5b_%d.pdf"%e)
    plt.show()
    
   

    
