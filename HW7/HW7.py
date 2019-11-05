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

t=np.linspace(0,20)

y0=[N1_0,R_0]

finaldent=[]
for i in range(0,len(gammas)):
    pars['gamma1']=gammas[i]
    pars['m1']=ms[i]
    y=integrate.odeint(alone,y0,t,args=(pars,))
    finaldent.append(y[-1,0])
    plt.plot(t,y[:,0])
    
plt.plot(t,t*0, linestyle="--", color="r")
plt.title("Bacteria " + str(i+1))
plt.xlabel("Time")
plt.ylabel("Density")
#plt.ylim([0,max(y[:,0])+1])
plt.show()
        
    #plt.plot(t,y[:,1],color="green")
    #plt.legend(["R"])
    #plt.xlabel("Time")
    #plt.ylabel("Density")
    #plt.show()

####Problematico
t=np.linspace(0,0.3)
pars['gamma1']=10
pars['m1']=100
y=integrate.odeint(alone,[10,50],t,args=(pars,))
y[:,0]
plt.plot(t,y[:,0])
plt.title("Bacteria 3")
plt.xlabel("Time")
plt.ylabel("Density")
plt.ylim([-1,max(y[:,0])+1])
plt.show()

plt.plot(t,y[:,1],color="green")
plt.legend(["R"])
plt.xlabel("Time")
plt.ylabel("Density")
plt.show()

###########################################################################

def compt(y0,t,pars):
    
    N1=y0[0]
    N2=y0[1]
    R=y0[2]
    
    dN1dt=pars['gamma1']*pars['e1']*N1*R-pars['m1']*N1
    dN2dt=pars['gamma2']*pars['e2']*N2*R-pars['m2']*N2
    dRdt=pars['w']*pars['R0']-pars['w']*R-pars['gamma1']*N1*R-pars['gamma2']*N2*R
    
    dydt=[dN1dt,dN2dt,dRdt]
    
    return dydt

pars={}
pars['e1']=0.5
pars['e2']=0.5
pars['w']=0.1
pars['R0']=100

N1_0=10
N2_0=10
R_0=50

gammas=[0.1,2,4,2,1,0.3]
ms=[1,4,10,5,10,1]
ind=[1,2,4,5,6,7]

t=np.linspace(0,20)


y0=[N1_0,N2_0,R_0]

N1_relative=np.zeros([len(gammas),len(gammas)])

for i in range(0,len(gammas)):
    for j in range(0,len(gammas)):
        pars['gamma1']=gammas[i]
        pars['gamma2']=gammas[j]
        pars['m1']=ms[i]
        pars['m2']=ms[j] 
        y=integrate.odeint(compt,y0,t,args=(pars,))
        plt.plot(t,y[:,0])
        plt.plot(t,y[:,1])
        plt.title("N1 = Bacteria " + str(ind[i]) + ", N2 = Bacteria " + str(ind[j]))
        plt.legend(["N1","N2"])
        plt.xlabel("Time")
        plt.ylabel("Density")
        plt.ylim([-1,max(list(y[:,1]) + list(y[:,0]))+1])
        plt.show()
        
        plt.plot(t,y[:,2],color="green")
        plt.legend(["R"])
        plt.xlabel("Time")
        plt.ylabel("Density")
        plt.ylim([-1,max(y[:,2])+1])
        plt.show()
        
        y0tmp=y[-1,0]
        y1tmp=y[-1,1]
        if(y0tmp<0):
            y0tmp=0
        if(y1tmp<0):
            y1tmp=0
        N1_relative[i,j]=y0tmp/(y0tmp+y1tmp)
        
#ticks=[1,2,4,5,6,7]
N1_relative
neg=plt.imshow(N1_relative)
plt.colorbar(neg)
plt.xlabel("N1")
plt.ylabel("N2")
plt.title("Proportion of Bacteria N1 at end of experiment")
plt.xticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.yticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.show()

gammavsgamma=np.zeros([len(gammas),len(gammas)])
for i in range(0,len(gammas)):
    for j in range(0,len(gammas)):
        a=(gammas[i]>gammas[j])*1
        b=(gammas[i]<gammas[j])*-1
        c=(gammas[i]==gammas[j])*0
        tot=a+b+c
        gammavsgamma[i,j]=tot        

neg=plt.imshow(gammavsgamma)
plt.colorbar(neg)
plt.xlabel("N1")
plt.ylabel("N2")
plt.title("Comparison between gamma values")
plt.xticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.yticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.set_cmap('RdBu')
plt.show()

mvsm=np.zeros([len(ms),len(ms)])
for i in range(0,len(ms)):
    for j in range(0,len(ms)):
        a=(ms[i]>ms[j])*1
        b=(ms[i]<ms[j])*-1
        c=(ms[i]==ms[j])*0
        tot=a+b+c
        mvsm[i,j]=tot        

neg=plt.imshow(mvsm)
plt.colorbar(neg)
plt.xlabel("N1")
plt.ylabel("N2")
plt.title("Comparison between m values")
plt.xticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.yticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.set_cmap('RdBu')
plt.show()


ratiovsratio=np.zeros([len(ms),len(ms)])
for i in range(0,len(ms)):
    for j in range(0,len(ms)):
        a=(gammas[i]/ms[i]>gammas[j]/ms[j])*1
        b=(gammas[i]/ms[i]<gammas[j]/ms[j])*-1
        c=(gammas[i]/ms[i]==gammas[j]/ms[j])*0
        tot=a+b+c
        ratiovsratio[i,j]=tot        

neg=plt.imshow(ratiovsratio)
plt.colorbar(neg)
plt.xlabel("N1")
plt.ylabel("N2")
plt.title("Comparison between m values")
plt.xticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.yticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.set_cmap('RdBu')
plt.show()


#B6vsB1 mismo ratio pero gana B6. Ambos valores son mas grandes para B6. En experimentos solos. B6 va mas abajo
#B4vsB5 mismo ratio pero gana B4. Tiene los valores mas grandes. En experimentos solos B4 va mas abajo. Interesaaaante.
#Esto cumplira para todo?

del finaldent[1]

dentvsdent=np.zeros([len(ms),len(ms)])
for i in range(0,len(ms)):
    for j in range(0,len(ms)):
        a=(finaldent[i]>finaldent[j])*1
        b=(finaldent[i]<finaldent[j])*-1
        c=(finaldent[i]==finaldent[j])*0
        tot=a+b+c
        dentvsdent[i,j]=tot        

neg=plt.imshow(dentvsdent)
plt.colorbar(neg)
plt.xlabel("N1")
plt.ylabel("N2")
plt.title("Comparison of final densities in isolated experiments")
plt.xticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.yticks(range(0,6), ('B1', 'B2', 'B4', 'B5', 'B6','B7'))
plt.set_cmap('RdBu')
plt.show()

finaldent #TAMPOCO
#Ratio es el que mejor predice.


