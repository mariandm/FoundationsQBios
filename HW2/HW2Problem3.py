#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 15:33:49 2019

@author: mariandm
"""

#Problem 3
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import *

##Subpregunta i
#Hacer plots como en la subfigura 2 pero ahora sin quitar el iptg 
#Hacer nullclines pero con diferentes valores de iptg
#Ver como cambian los steady states con el uso de iptg

#Parameters
pars={}
pars['a1']=156.25
pars['a2']=15.6
pars['beta']=2.5
pars['gamma']=1
pars['K']=2.9618*10**-5
pars['IPTG']=0.02 #Segun figura 4 en paper
pars['eta']=2.0015

x0=[0,0]
t0=0
tf=15
t=np.linspace(t0,tf)

f, (ax1, ax2) = plt.subplots(2, 1)

x=integrate.odeint(ModelB,x0,t,args=(pars,))

#plt.plot(t,x[:,0],linewidth=3)
ax1.plot(t,x[:,1],linewidth=3,color='r')
ax1.set_xlabel("Time (hr)")
ax1.set_ylabel("Concentration (nM)")
ax1.set_title("Production of protein v with IPTG induction")
ax1.set_ylim([0,16])


pars['IPTG']=0
x=integrate.odeint(ModelB,x0,t,args=(pars,))

#plt.plot(t,x[:,0],linewidth=3)
ax2.plot(t,x[:,1],linewidth=3,color='r')
ax2.set_xlabel("Time (hr)")
ax2.set_ylabel("Concentration (nM)")
ax2.set_title("Production of protein v, control")
ax2.set_ylim([0,16])

f.tight_layout()
f.savefig("fig4.pdf", bbox_inches='tight')
plt.close()



#It stays when running for 15 hours

##################################################

##Subpregunta ii
#Parameters
pars={}
pars['a1']=156.25
pars['a2']=15.6
pars['beta']=2.5
pars['gamma']=1
pars['K']=2.9618*10**-5
pars['IPTG']=0 #Segun figura 4 en paper
pars['eta']=2.0015

#Experiment: IPTG 5 hours, then no IPTG 10 hours

f, (ax1, ax2) = plt.subplots(2, 1)

x0=[1,1]
t0=0
tf=5
t=np.linspace(t0,tf)

pars['IPTG']=0.02
x=integrate.odeint(ModelB,x0,t,args=(pars,))
upre=x[:,0]
vpre=x[:,1]
u0=x[:,0][-1]
v0=x[:,1][-1]
pars['IPTG']=0
t2=np.linspace(tf,tf*3)
x=integrate.odeint(ModelB,[u0,v0],t,args=(pars,))
ucon=list(upre)+list(x[:,0])
vcon=list(vpre)+list(x[:,1])
time=list(t)+list(t2)

#plt.plot(time,ucon,linewidth=3,color='b')
ax1.plot(time,vcon,linewidth=3,color='r')
ax1.set_xlabel("Time (hrs)")
ax1.set_ylabel("Concentration (nM)")
ax1.set_title("IPTG Induction Removal Experiment")
ax1.set_ylim([0,16])
x = np.arange(0.0, 5, 0.1)
y1 = 0 * x +16
ax1.fill_between(x, 0, y1,color="lightgray")
#plt.legend(["u",'v'])


##CONTROL: No IPTG
#ESto pasa porque u tiene una a1 mucho mas grande en comparacion de a2

pars['IPTG']=0
x0=[1,1]
t0=0
tf=5
t=np.linspace(t0,tf)

x=integrate.odeint(ModelB,x0,t,args=(pars,))
upre=x[:,0]
vpre=x[:,1]
u0=x[:,0][-1]
v0=x[:,1][-1]
pars['IPTG']=0
t2=np.linspace(tf,tf*3)
x=integrate.odeint(ModelB,[u0,v0],t,args=(pars,))
ucon=list(upre)+list(x[:,0])
vcon=list(vpre)+list(x[:,1])
time=list(t)+list(t2)

#plt.plot(time,ucon,linewidth=3,color='b')
ax2.plot(time,vcon,linewidth=3,color='r')
ax2.set_xlabel("Time (hrs)")
ax2.set_ylabel("Concentration (nM)")
#plt.legend(["u",'v'])
ax2.set_title("IPTG Induction Removal Control")
ax2.set_ylim([0,16])


f.tight_layout()
f.savefig("fig5.pdf", bbox_inches='tight')
plt.close()

###########################################

#subpregunta iii
pars={}
pars['a1']=156.25
pars['a2']=15.6
pars['beta']=2.5
pars['gamma']=1
pars['K']=2.9618*10**-5
pars['IPTG']=0 #Segun figura 4 en paper
pars['eta']=2.0015

x0=[1,1]
t0=0
tf=17
t=np.linspace(t0,tf)

IPTG_conditions=[0,10**(-6),10**(-5.5), 10**(-5),10**(-4.9),10**(-4.8),10**(-4.7),10**(-4.6),10**(-4.5),10**(-4), 10**(-3), 10**(-2)]
uline=[]
vline=[]

for i in IPTG_conditions:
    pars['IPTG']=i
    x=integrate.odeint(ModelB,x0,t,args=(pars,))
    ufinal=x[:,0][-1]
    vfinal=x[:,1][-1]
    uline.append(ufinal)
    vline.append(vfinal)

#plt.plot(np.log10(IPTG_conditions),uline,linewidth=3,color='b')
plt.plot(np.log10(IPTG_conditions),vline,linewidth=3,color='r',marker="o")
plt.xlabel("log10 IPTG (M)")
plt.ylabel("v concentration (nM)")
#plt.legend(["u",'v'])
plt.savefig("fig6.pdf")
plt.close()
    
#NULLCLINES

vvec=np.linspace(0,15,70)


pars['IPTG']=10**(-4.75)
C = ((1+pars['IPTG']/pars['K'])**(pars['eta']*pars['gamma']))
unullcline = lambda v: pars['a1']/(1+v**pars['beta'])
vnullcline = lambda v: (C*pars['a2']/v-C)**(1/pars['gamma'])
plt.plot(unullcline(vvec), vvec, linewidth=3)
plt.plot(vnullcline(vvec), vvec, linewidth=3)
plt.plot(151.04,0.26,marker="o",color="black")
plt.plot(0.1949,14.5,marker="o",color="black")
plt.plot(12.56,2.65,marker="o",color="black")
plt.xlabel('u concentration (nMol)',fontsize=20)
plt.ylabel('v concentration (nMol)',fontsize=20)
plt.legend(["u",'v'])
plt.savefig("Fig7a.pdf")
plt.close()


x0=[[50,14],[0,2],[2,0],[100,10],[150,14],[160,6]]

# Define a vector of protein and mRNA values
#vvec = np.arange(0,60,1)*0.1
# Plot the nullclines
plt.plot(unullcline(vvec), vvec, linewidth=3, color="lightgray")
plt.plot(vnullcline(vvec), vvec, linewidth=3, color="lightgray")
for i in range(0,len(x0)):
    x=integrate.odeint(ModelB,x0[i],t,args=(pars,))
    plt.plot(x[:,0],x[:,1],linewidth=2)
# Label the plot
plt.plot(151.04,0.26,marker="o",color="black")
plt.plot(0.1949,14.5,marker="o",color="black")
plt.plot(12.56,2.65,marker="o",color="black")
plt.legend(['nullclines'],fontsize=18)
#plt.xlim(min(mvec),max(mvec))
#plt.ylim(min(pvec),max(pvec))
plt.xlabel('u concentration (nMol)',fontsize=20)
plt.ylabel('v concentration (nMol)',fontsize=20)
plt.savefig("Fig7b.pdf")
plt.close()

t0=0
tf=20
t=np.linspace(t0,tf)

for i in range(0,len(x0)):
    x=integrate.odeint(ModelB,x0[i],t,args=(pars,))
    plt.plot(t,x[:,0],linewidth=3)
plt.ylabel("u concentration (nMol)")
plt.xlabel("Times (hrs)")
plt.title("u concentration at different time points")
plt.savefig("Figure7c_uconc.pdf")
plt.close()
    
for i in range(0,len(x0)):
    x=integrate.odeint(ModelB,x0[i],t,args=(pars,))
    plt.plot(t,x[:,1],linewidth=3)
plt.ylabel("v concentration (nMol)")
plt.xlabel("Times (hrs)")
plt.title("v concentration at different time points")
plt.savefig("Figure7d_vconc.pdf")
plt.close()
