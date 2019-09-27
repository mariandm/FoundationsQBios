#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 21:58:17 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import *

A00=np.arange(0,50,10)
A01=np.arange(50,200,50)
A02=np.arange(200,400,100)
A0s=np.array(list(A00)+list(A01)+list(A02))
A0s[0]=1
fixes=np.zeros(len(A0s))

N = 500 # Number of individuals
s=0.05
nsamples=100

for j in range(0,len(A0s)):
    cfix=0
    A_0=A0s[j]
    for i in range(0,nsamples):  
        A,fracA,fixation,currg=moran(N,A_0,s)
        cfix+=fixation
    fixes[j]=cfix

pis=(1-(1+s)**(-A0s))/(1-(1+s)**(-N))
fixes=np.array(fixes)/nsamples


plt.plot(A0s,fixes,marker="o",color="rebeccapurple")
plt.plot(A0s,pis,marker="x",color="thistle")
plt.plot(np.ones(len(np.linspace(0,1.05)))*1/s,np.linspace(0,1.05),color="red",linestyle="--")
plt.xlabel("Initial number of A individuals")
plt.ylabel("Probability of Fixation")
plt.legend(["Experiments","Theory","Critical value: 1/s"])
plt.xlim([0,320])
plt.ylim([0,1.05])
plt.savefig("Figure10b.pdf")

#If the mutant makes the leap and survives the genetic drift, the probability of fixation is very high
