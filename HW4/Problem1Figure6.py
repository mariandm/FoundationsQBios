#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 23:21:17 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
numgen=1000
N = 100
Tmat = np.zeros([N+1, N+1])
for jj in range(N + 1):
    currp = jj/N
    xvec = np.array(range(N + 1))
    Tmat[:, jj] = stats.binom.pmf(xvec, N, currp)

nvec = np.arange(N+1)
X0vec = np.zeros((N+1,1))
X0vec[int(N/2+1),0]=1
theormean = np.zeros(numgen+1)
theorvar = np.zeros(numgen+1)
theorp = np.zeros((numgen+1,N+1))
Xcur = np.matrix(X0vec) #note these are numpy matrices, not arrays
TMAT = np.matrix(Tmat)
theormean[0]=np.dot(X0vec.flatten(),nvec)
theorp[0,:]=X0vec.flatten()
for i in range(1,numgen+1):
    Xcur = TMAT*Xcur
    probvec = np.array(Xcur).flatten()
    theormean[i]=np.dot(np.array(probvec).flatten(),nvec)
    secondmoment = np.dot(np.array(probvec).flatten(),nvec**2)
    theorvar[i] = secondmoment-theormean[i]**2
    theorp[i] = probvec
  
x=np.ones([numgen+1,N+1])*range(N+1)
plt.plot(x[5],theorp[5])
plt.plot(x[10],theorp[10])
plt.plot(x[20],theorp[20])
plt.plot(x[50],theorp[50])
plt.xlabel("Number of A individuals")
plt.ylabel("Probability")
plt.legend(["Generation 5","Generation 10","Generation 20","Generation 50"])
#plt.savefig("Figure7a.pdf")
plt.close()

plt.plot(range(0,numgen+1),theorp[:,-1])


nvec = np.arange(N+1)
X0vec = np.zeros((N+1,1))
X0vec[80,0]=1
theormean = np.zeros(numgen+1)
theorvar = np.zeros(numgen+1)
theorp = np.zeros((numgen+1,N+1))
Xcur = np.matrix(X0vec) #note these are numpy matrices, not arrays
TMAT = np.matrix(Tmat)
theormean[0]=np.dot(X0vec.flatten(),nvec)
theorp[0,:]=X0vec.flatten()
for i in range(1,numgen+1):
    Xcur = TMAT*Xcur
    probvec = np.array(Xcur).flatten()
    theormean[i]=np.dot(np.array(probvec).flatten(),nvec)
    secondmoment = np.dot(np.array(probvec).flatten(),nvec**2)
    theorvar[i] = secondmoment-theormean[i]**2
    theorp[i] = probvec
    
plt.plot(range(0,numgen+1),theorp[:,-1])



nvec = np.arange(N+1)
X0vec = np.zeros((N+1,1))
X0vec[20,0]=1
theormean = np.zeros(numgen+1)
theorvar = np.zeros(numgen+1)
theorp = np.zeros((numgen+1,N+1))
Xcur = np.matrix(X0vec) #note these are numpy matrices, not arrays
TMAT = np.matrix(Tmat)
theormean[0]=np.dot(X0vec.flatten(),nvec)
theorp[0,:]=X0vec.flatten()
for i in range(1,numgen+1):
    Xcur = TMAT*Xcur
    probvec = np.array(Xcur).flatten()
    theormean[i]=np.dot(np.array(probvec).flatten(),nvec)
    secondmoment = np.dot(np.array(probvec).flatten(),nvec**2)
    theorvar[i] = secondmoment-theormean[i]**2
    theorp[i] = probvec
    
plt.plot(range(0,numgen+1),theorp[:,-1])

plt.xlabel("Generations")
plt.ylabel("Probability of Fixation")
plt.title("Probability of Fixation given different initial densities of A, N=100")
plt.legend(["A0=50","A0=80","A0=20"])
plt.savefig("Figure7b.pdf")