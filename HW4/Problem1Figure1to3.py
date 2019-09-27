#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 11:29:44 2019

@author: mariandm
"""
import numpy as np
import matplotlib.pyplot as plt

#Problem 1

#Figure 1: 3 experiments with fixation
fixation=0
for rep in range(0,30):
    N = 200 # Number of individuals
    numA_0 = N/2 #Fraction that are A
    g = 100 # Generations
    fracA = np.zeros(g+1)
    fracA[0] = numA_0/N
    for i in range(1,g+1):
        nextA = np.random.binomial(N,fracA[i-1])
        fracA[i]=nextA/N
    if fracA[g]==1 or fracA[g]==0:
        fixation+=1
        plt.plot(np.arange(g+1),fracA,marker='.',markersize=8)
        if fixation==3:
            break
plt.ylim([0,1])
plt.ylabel("Fraction of type A")
plt.xlabel("Generations")
#plt.savefig("HW4Fig1.pdf")
plt.close()

#Figure 2: Approximate number of time steps
fixationTimes=[]
for rep in range(0,1000):
    N = 200 # Number of individuals
    numA_0 = N/2 #Fraction that are A
    g = 1500 # Generations
    fracA = np.zeros(g+1)
    fracA[0] = numA_0/N
    for i in range(1,g+1):
        nextA = np.random.binomial(N,fracA[i-1])
        fracA[i]=nextA/N
        if fracA[i]==1 or fracA[i]==0:
            fixationTimes.append(i)
            break

fixtimemean = np.mean(fixationTimes)
fixtimemedian = np.median(fixationTimes)

a=plt.hist(fixationTimes,bins=20,color="lightblue")
#Mean
plt.plot(np.ones(len(a[0]))*fixtimemean,a[0],color="red", linestyle="--")
#Median
plt.plot(np.ones(len(a[0]))*fixtimemedian,a[0],color="blue", linestyle="--")

plt.legend(["Mean","Median","Experiments"])
plt.ylabel("Number of experiments")
plt.xlabel("Generations")
#plt.savefig("HW4Fig2.pdf")
plt.close()

#Exponential
expo = np.random.exponential(fixtimemean,size=1000)
plt.hist([fixationTimes,expo],bins=20, rwidth=0.9,color=["lightblue","gray"])
plt.xlim([0,1000])
plt.ylabel("Number of experiments")
plt.xlabel("Generation at fixation")
plt.legend(["Experiments","Exponential Distribution"])
#plt.savefig("HW4Fig2b.pdf")
plt.close()

#Figure3

def fixtimes(nsamples,N,A0,generations):
    fixation = 0
    elimination = 0
    fixationTimes=[]
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
                fixationTimes.append(i)
                elimination+=1
                break
    return(fixationTimes,fixation,elimination)

As=[25,50,100,150,175]
nsamples=1000
N=200
generations=3000

fixis=[]
colors=["darkorange","darkcyan","lightblue","mediumpurple","palevioletred"]

fig, axs = plt.subplots(1,5,figsize=(15, 4))
prob_fix=np.zeros(len(As))
for i in range(0,len(As)):
    A0=As[i]
    percentage=A0/N
    fixis,fixation,elimination=fixtimes(nsamples,N,A0,generations)
    axs[i].hist(fixis,bins=20,color=colors[i])
    axs[i].set_xlim([0,1000])
    axs[i].set_ylabel("Number of experiments")
    axs[i].set_xlabel("Generation at fixation")
    axs[i].set_title('Fraction of Type A =%.2f'%percentage)
    prob_fix[i]=fixation/nsamples
fig.tight_layout()
#fig.savefig("Figure3.pdf", bbox_inches='tight')
plt.show()


plt.scatter(np.array(As)/N,prob_fix,s=40,color="rebeccapurple")
plt.scatter(np.array(As)/N,np.array(As)/N, marker="x",s=40,color="thistle")
plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel("Initial fraction of Type A")
plt.ylabel("Probabiity of Fixation")
plt.legend(["Experiments","Theoretical"])
#plt.savefig("Figure4.pdf")
plt.show()

