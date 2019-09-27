#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 19:51:35 2019

@author: mariandm
"""
import numpy as np
import matplotlib.pyplot as plt
from Functions import *

#plt.ylim([0,1])


N = 500 # Number of individuals
A_0 = 1 #Fraction that are A
s=0.05 #A reproduction advantage
nsamples=10


fixationTimes=[[],[]] # fixation is fixationTimes[1], elimination is fixationTimes[0]
cfix=0
while(cfix<3):  
    A,fracA,fixation,currg=moran(N,A_0,s)
    cfix+=fixation
    if(fixation):
        plt.plot(range(0,currg+1),fracA)
    fixationTimes[fixation].append(currg)
plt.xlabel("Generations")
plt.ylabel("Fraction of A individuals")
plt.legend(["Run1","Run2","Run3"])
#plt.savefig("Problem2Figure8.pdf")

print(len(fixationTimes[0])+len(fixationTimes[1]))

