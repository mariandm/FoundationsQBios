#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 12:26:45 2019

@author: mariandm
"""

#Problem 1

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from Functions import *

#Parameters
pars={}
pars['beta_minus']=5
pars['beta_plus']=50
pars['K']= 25
pars['n']=150
pars['alpha']=1


####FIGURA 1
xaxis=np.array([-pars['K'],-0.0001, 0, pars['K']])
yaxis= np.heaviside(xaxis,1)*(pars['beta_plus']-pars['beta_minus'])+pars['beta_minus']
plt.plot(xaxis+pars['K'],yaxis,color="lightgray")
#highest x
plt.plot(xaxis+pars['K'],(xaxis+pars['K'])*(pars['beta_plus']/pars['K']))
#medium x
plt.plot(xaxis+pars['K'],((xaxis+pars['K'])*(pars['beta_minus']/pars['K'])+(xaxis+pars['K'])*(pars['beta_plus']/pars['K']))/2)
#lowest x
plt.plot(xaxis+pars['K'],(xaxis+pars['K'])*(pars['beta_minus']/pars['K']))
#biology
plt.plot(xaxis+pars['K'],(xaxis+pars['K'])*0.5)
plt.legend(["Production",r'$\alpha ≥ \beta^+/K$',r'$\beta_-/K < \alpha < \beta^+/K$',r'$\alpha ≤ \beta_-/K$',r'$\alpha = 0.5$'])
plt.xlabel("Time (hrs)")
plt.ylabel("Change, dx/dt (nM/hr)")

plt.savefig('Figure1.pdf', bbox_inches='tight')

plt.show()


####FIGURA 2 #Hacer muchos otros numeros para n
f, (ax1s, ax2s) = plt.subplots(2, 2)

alpha1=2
alpha2=1.25
alpha3=0.2
pars['n']=100

ax1=ax1s[0]

x=np.linspace(0,60)
y=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n']) + pars['beta_minus']
ax1.plot(x,y, color="lightgray")
otroy=alpha1*x
ax1.plot(x,otroy)
otroy=alpha2*x
ax1.plot(x,otroy)
otroy=alpha3*x
ax1.plot(x,otroy)
ax1.set_xlabel("Time (hrs)")
ax1.set_ylabel("Change, dx/dt (nM/hr)")
ax1.set_title("n=100")
ax1.legend(["Production",r'$\alpha = 2$',r'$\alpha = 1.25$',r'$\alpha = 0.2$'],prop={'size': 6})

ax2=ax1s[1]

pars['n']=10
x=np.linspace(0,60)
y=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n']) + pars['beta_minus']
ax2.plot(x,y, color="lightgray")
otroy=alpha1*x
ax2.plot(x,otroy)
otroy=alpha2*x
ax2.plot(x,otroy)
otroy=alpha3*x
ax2.plot(x,otroy)
ax2.set_xlabel("Time (hrs)")
ax2.set_ylabel("Change, dx/dt (nM/hr)")
ax2.set_title("n=10")
ax2.legend(["Production",r'$\alpha = 2$',r'$\alpha = 1.25$',r'$\alpha = 0.2$'],prop={'size': 6})

ax3=ax2s[0]
pars['n']=5
x=np.linspace(0,60)
y=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n']) + pars['beta_minus']
ax3.plot(x,y, color="lightgray")
otroy=alpha1*x
ax3.plot(x,otroy)
otroy=alpha2*x
ax3.plot(x,otroy)
otroy=alpha3*x
ax3.plot(x,otroy)
ax3.set_xlabel("Time (hrs)")
ax3.set_ylabel("Change, dx/dt (nM/hr)")
ax3.set_title("n=5")
ax3.legend(["Production",r'$\alpha = 2$',r'$\alpha = 1.25$',r'$\alpha = 0.2$'],prop={'size': 6})

ax4=ax2s[1]
pars['n']=3
x=np.linspace(0,60)
y=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n']) + pars['beta_minus']
ax4.plot(x,y, color="lightgray")
otroy=alpha1*x
ax4.plot(x,otroy)
otroy=alpha2*x
ax4.plot(x,otroy)
otroy=alpha3*x
ax4.plot(x,otroy)
ax4.set_xlabel("Time (hrs)")
ax4.set_ylabel("Change, dx/dt (nM/hr)")
ax4.set_title("n=3")
ax4.legend(["Production",r'$\alpha = 2$',r'$\alpha = 1.25$',r'$\alpha = 0.2$'],prop={'size': 6})

f.tight_layout()
f.savefig("fig2.pdf", bbox_inches='tight')


