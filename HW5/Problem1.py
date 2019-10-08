#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:33:23 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import *

    
    
Xtot =1000
Xm=350
Xmp=200
x=[Xm,Xmp]
t=np.linspace(0,1000)

vr = 1
R = 20
Kr = 2
vd = 100

vb=1
B=40
Kb = 200
A=2
vma=100
va = vma/(1+A)

tchanges=[0,200,400,600]
Achanges=[1,0.5,2,0.25]


a = integrate.odeint(chemotaxis,x,t,args=(vr,R,Kr,vd,vb,B,Kb,vma,Xtot,tchanges,Achanges))
Xzero=np.array([Xtot]*len(a[:,0]))
Xzero = Xzero-a[:,0]-a[:,1]
plt.plot(t,a[:,0])
plt.plot(t,a[:,1])
plt.plot(t,Xzero)
plt.plot(t,[200]*len(Xzero),color="red", linestyle="--",linewidth=0.8)
plt.scatter(tchanges,np.zeros(len(tchanges)),color="k")
plt.xlabel("Time (min)")
plt.ylabel("Complex Concentration (nM)")
plt.ylim([0,700])
plt.xlim([0,1000])
plt.legend(["Xm","Xm,p","X0","Equilibrium","A pulses"],fontsize=8)
plt.savefig("Figure1.pdf")
plt.close()

t=np.linspace(0,500)
tchanges=[0,200]

Aat=[0.25,0.5,1.5,2]
fig, axs = plt.subplots(2,2,figsize=(6, 8))

for i in range(0,len(Aat)):
    ax0 = int(i/2)
    ax1 = i - ax0*2
    Aa=Aat[i]
    Achanges=[1,Aa]
    a = integrate.odeint(chemotaxis,x,t,args=(vr,R,Kr,vd,vb,B,Kb,vma,Xtot,tchanges,Achanges))
    Xzero=np.array([Xtot]*len(a[:,0]))
    Xzero = Xzero-a[:,0]-a[:,1]
    axs[ax0][ax1].plot(t,a[:,0])
    axs[ax0][ax1].plot(t,a[:,1])
    axs[ax0][ax1].plot(t,Xzero)
    axs[ax0][ax1].plot(t,[200]*len(Xzero),color="red", linestyle="--",linewidth=0.8)
    axs[ax0][ax1].scatter(tchanges,np.zeros(len(tchanges)),color="k")
    y = range(0,int(max(np.ndarray.flatten(a))))
    #axs[ax0][ax1].plot(np.ones([len(y),len(tchanges)])*tchanges,y,color="red", linestyle="--")
    axs[ax0][ax1].set_xlabel("Time (min)")
    axs[ax0][ax1].set_ylabel("Complex Concentration (nM)")
    axs[ax0][ax1].legend(["Xm","Xm,p","X0","Equilibrium","A pulses"],fontsize=8)
    axs[ax0][ax1].set_title("A attract = %.2f"%Aa)
    axs[ax0][ax1].set_ylim([0,600])
    axs[ax0][ax1].set_xlim([0,500])
fig.tight_layout()
fig.savefig("Figure2.pdf")
plt.close()


fig, axs = plt.subplots(2,2,figsize=(6, 7))
tchanges=[0,200]
Achanges=[1,2]
Rs=[5,10,15,20]
for i in range(0,len(Rs)):
    ax0 = int(i/2)
    ax1 = i - ax0*2
    R=Rs[i]
    a = integrate.odeint(chemotaxis,x,t,args=(vr,R,Kr,vd,vb,B,Kb,vma,Xtot,tchanges,Achanges))
    Xzero=np.array([Xtot]*len(a[:,0]))
    Xzero = Xzero-a[:,0]-a[:,1]
    #axs[ax0][ax1].plot(t,a[:,0])
    axs[ax0][ax1].plot(t,a[:,1],color="darkorange")
    axs[ax0][ax1].plot(t,[a[:,1][-1]]*len(Xzero),color="red", linestyle="--",linewidth=0.8)
    axs[ax0][ax1].scatter(tchanges,np.zeros(len(tchanges)),color="k")
    axs[ax0][ax1].set_xlabel("Time (min)")
    axs[ax0][ax1].set_ylabel("Complex Concentration (nM)")
    axs[ax0][ax1].legend(["Xm,p","Equilibrium","A pulses"],fontsize=8)
    axs[ax0][ax1].set_title("R = %d"%R)
    axs[ax0][ax1].set_xlim([0,500])
    axs[ax0][ax1].set_ylim([0,250])
    print(t[np.where(a[:,1]==min(a[:,1]))[0]])
    print(t[max(np.where(a<((a[:,1][-1] - min(a[:,1]))/2 + min(a[:,1])))[0])] - t[np.where(a[:,1]==min(a[:,1]))[0]])
fig.tight_layout()
fig.savefig("Figure3.pdf")
plt.close()


t=np.linspace(0,500)
tchanges=[0,200]
R=20
Aat=[4,8]
fig, axs = plt.subplots(1,2,figsize=(8, 3))

for i in range(0,len(Aat)):
    Aa=Aat[i]
    Achanges=[1,Aa]
    a = integrate.odeint(chemotaxis,x,t,args=(vr,R,Kr,vd,vb,B,Kb,vma,Xtot,tchanges,Achanges))
    Xzero=np.array([Xtot]*len(a[:,0]))
    Xzero = Xzero-a[:,0]-a[:,1]
    axs[i].plot(t,a[:,0])
    axs[i].plot(t,a[:,1])
    axs[i].plot(t,Xzero)
    axs[i].plot(t,[200]*len(Xzero),color="red", linestyle="--",linewidth=0.8)
    axs[i].scatter(tchanges,np.zeros(len(tchanges)),color="k")
    axs[i].set_xlabel("Time (min)")
    axs[i].set_ylabel("Complex Concentration (nM)")
    axs[i].legend(["Xm","Xm,p","X0","Equilibrium","A pulses"],fontsize=8)
    axs[i].set_title("A attract = %d"%Aa)
    axs[i].set_xlim([0,500])
    axs[i].set_ylim([0,900])
fig.tight_layout()
fig.savefig("Figure4.pdf")
plt.close()