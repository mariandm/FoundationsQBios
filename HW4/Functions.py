#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 21:18:49 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt


def moran(N,A_0,s):
    A=[]
    A.append(A_0)
    fracA = []
    fracA.append(A_0/N)
    currg=0
    fixation=0
    elimination=0
    while(fixation==0 and elimination==0):
        currg+=1
        currA=A[currg-1]
        for j in range(0,N):
            ReproProbA = currA*(1+s)/(currA*s+N)
            DeathProbA = currA/N
            currA+=np.random.binomial(1,ReproProbA) - np.random.binomial(1,DeathProbA)
        A.append(currA)
        fracA.append(currA/N)
        if(fracA[currg]==1):
            fixation=1
        if(fracA[currg]==0):
            elimination=1
    return(A,fracA,fixation,currg)