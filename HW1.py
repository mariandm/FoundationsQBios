#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 17:54:43 2019

@author: mariandm
"""

#Problem 1. Simulating the Luria-Delbruck Experiment, One Generation
#Write a program to simulate just one generation of the LD experiment-stochastically.  
#Simulate C=500 cultures each of which has N=1,000 cells and μ= 10-3, i.e., a very high mutation rate.   
#What is the distribution of resistant mutants that you observe across all the cultures?  
#Are they similar or dissimilar to each other?  Specify your measurement of c(m), i.e., the number of cultures with m resistant mutants.  
#Is this distribution well fit by a Poisson distribution?  
#If so, what is the best fitshape parameter of the Poisson density function and how does that relate to the microscopic value of mutation you used to generate the output?  
#Finally, to what extent are the fluctuations “large” or “small”?

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

C=500 #Number of cultures
N=1000 #Initial bacteria concentration
mu=10**-3 #mutation rate
g=1 #Number of generations


def LDreplicates(C,N,mu,g):
    mTotal=[]
    for rep in range(0,C):
        S=N
        m=0
        for gen in range(0,g):
            #Replication
            S=S*2
            m=m*2
            #Generating mutations
            newm=np.sum(np.random.rand(1, S) <= mu)
            S=S-newm
            m=m+newm
        mTotal.append(m)
    return(mTotal)
    
###LO SILENCIAMOS PORQUE QUEREMOS QUEDARNOS CON ESTE SET
#Test1=LDreplicates(C,N,mu,g)

#####What is the distribution of resistant mutants that you observe across all the cultures? 
a = plt.hist(Test1, bins=max(Test1)+1)
plt.xlabel("Number of resistance mutants [m]")
plt.ylabel("Number of cultures with m resistance mutants [C(m)]")
plt.savefig('Fig1_E1_Cm.pdf', bbox_inches='tight')
plt.close()
#Como contesto esto?
#Sacar la media y la varianza

#Make a table of this

Test1_mutant = list(set(Test1))
Test1_mean=np.mean(Test1)
Test1_var=np.var(Test1)

#Mean and variance are similar
print("mutants")
print(Test1_mutant)
print("cmutant")
for i in Test1_mutant:
    print(Test1.count(i))
print("mean")
print(Test1_mean)
print("variance")
print(Test1_var)


#How to fit a poisson dist
#N*(2**g)

#the term in the prob is the expected number of mutants per culture
#En las notas de Joshua está que Nmu es la lambda
#mu por 2000 es dos entonces podemos probar valores similares 0.5,1,1.5,2,2.5,3

#Poisson_mu0=np.random.poisson(mu*N*(2**g),C)
Poisson_mumean=np.random.poisson(Test1_mean,C)

plt.hist([Test1,Poisson_mumean,Poisson_mu0],rwidth=1,bins=max(Test1)+1)
plt.legend(["Experimental data","Poisson lambda= data's mean","Poisson lambda = theory"])
plt.xlabel("Number of resistance mutants [m]")
plt.ylabel("Number of cultures with m resistance mutants [C(m)]")
plt.savefig('Fig2_E1_Cm.pdf', bbox_inches='tight')
plt.close()



#####Are they similar or dissimilar to each other?  Specify your measurement of c(m), i.e., the number of cultures with m resistant mutants.  
#Maybe do a count table
#print(Test1.count(0))

#####Is this distribution well fit by a Poisson distribution? 
#Ponle una poisson con la lambda mu y has un histograma intercalado

#####If so, what is the best fitshape parameter of the Poisson density function and how does that relate to the microscopic value of mutation you used to generate the output? 
#Probar varios valores de lambda cercanos a mi valor y hacer una prueba estadística. O enseñar algunos histogramas empalmados (tampoco hay que volverse locos y mostrar 10, con 9 xd)

#####Finally, to what extent are the fluctuations “large” or “small”?


