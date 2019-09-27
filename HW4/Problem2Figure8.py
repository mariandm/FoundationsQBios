#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 21:20:09 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import *

ss=np.arange(0.05,0.55,0.05)
fixes=np.zeros(len(ss))

N = 500 # Number of individuals
A_0 = 1 #Fraction that are A
nsamples=100

for j in range(0,len(ss)):
    cfix=0
    s=ss[j]
    for i in range(0,nsamples):  
        A,fracA,fixation,currg=moran(N,A_0,s)
        cfix+=fixation
    fixes[j]=cfix

pis=ss/(ss+1)
fixes=np.array(fixes)/nsamples

plt.plot(ss,fixes,marker="o",color="rebeccapurple")
plt.plot(ss,pis,marker="x",color="thistle")
plt.xlabel("S")
plt.ylabel("Probability of Fixation")
plt.legend(["Experiments","Theory"])
plt.xlim([0,0.5])
plt.savefig("Figure9.pdf")