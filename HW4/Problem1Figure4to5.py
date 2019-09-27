#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 19:24:13 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt

def fixtimes(nsamples,N,A0,generations):
    fixation = 0
    elimination = 0
    fixationTimes=[]
    eliminationTimes=[]
    for rep in range(0,nsamples):
        numA_0 = A0 #Fraction that are A
        g = generations # Generations
        fracA = np.zeros(g+1)
        fracA[0] = numA_0/N
        for i in range(1,g+1):
            nextA = np.random.binomial(N,fracA[i-1])
            fracA[i]=nextA/N
            if fracA[i]==1:
                fixationTimes.append(i)
                fixation+=1
                break
            if fracA[i]==0:
                eliminationTimes.append(i)
                elimination+=1
                break
    return(fixationTimes,eliminationTimes,fixation,elimination)
    
N=100
A0=50
nsamples=1000
generations=3000

fixationTimes,eliminationTimes,fixation,elimination=fixtimes(nsamples,N,A0,generations)



#FIGURE 5

fig, axs = plt.subplots(1,3,figsize=(15, 4))

#FIXATION
fixtimemean = np.mean(fixationTimes)
fixtimemedian = np.median(fixationTimes)
a=axs[0].hist(fixationTimes,bins=20,color="lightblue")
#Mean
axs[0].plot(np.ones(len(a[0]))*fixtimemean,a[0],color="red", linestyle="--")
#Median
axs[0].plot(np.ones(len(a[0]))*fixtimemedian,a[0],color="blue", linestyle="--")
axs[0].set_xlim([0,600])
axs[0].legend(["Mean","Median","Experiments"])
axs[0].set_ylabel("Number of experiments")
axs[0].set_xlabel("Generation at fixation")
axs[0].set_title("Fixation")
#plt.savefig("HW4Fig2a.pdf")

#ELIMINATION
elitimemean = np.mean(eliminationTimes)
elitimemedian = np.median(eliminationTimes)

a=axs[1].hist(eliminationTimes,bins=20,color="lightblue")
#Mean
axs[1].plot(np.ones(len(a[0]))*elitimemean,a[0],color="red", linestyle="--")
#Median
axs[1].plot(np.ones(len(a[0]))*elitimemedian,a[0],color="blue", linestyle="--")
axs[1].set_xlim([0,600])
axs[1].legend(["Mean","Median","Experiments"])
axs[1].set_ylabel("Number of experiments")
axs[1].set_xlabel("Generation at elimination")
axs[1].set_title("Elimination")
#plt.savefig("HW4Fig2a.pdf")

#BOTH 
bothTimes=fixationTimes+eliminationTimes
bothtimemean = np.mean(bothTimes)
bothtimemedian = np.median(bothTimes)

a=axs[2].hist(bothTimes,bins=20,color="lightblue")
#Mean
axs[2].plot(np.ones(len(a[0]))*bothtimemean,a[0],color="red", linestyle="--")
#Median
axs[2].plot(np.ones(len(a[0]))*bothtimemedian,a[0],color="blue", linestyle="--")
axs[2].set_xlim([0,600])
axs[2].legend(["Mean","Median","Experiments"])
axs[2].set_ylabel("Number of experiments")
axs[2].set_xlabel("Generation at fixation or elimination")
axs[2].set_xlim([0,600])
axs[2].set_title("Either Fixation or Elimination")

fig.savefig("Figure5.pdf")
plt.close(fig)

print("FIXATION")
print(fixation)
print(fixtimemean)
print("ELIMINATION")
print(elimination)
print(elitimemean)
print("BOTH")
print(bothtimemean)


################################################################

