#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:51:40 2019

@author: mariandm
"""

import numpy as np
from scipy import integrate

def selfregulation(x,t,pars):
    dxdt=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n'])-pars['alpha']*x + pars['beta_minus']
    return(dxdt)
    
def selfproduction(x,t,pars):
    dxdt=(pars['beta_plus']*x**pars['n'])/(pars['K']**pars['n']+x**pars['n']) + pars['beta_minus']
    return(dxdt)
    
def selfdegradation(x,t,pars):
    dxdt=-pars['alpha']*x
    return(dxdt)
    

def ToggleBox1(x,t,pars):
    u=x[0]
    v=x[1]
    dudt = pars['a1'] / (1 + v**pars['beta']) - u
    dvdt = pars['a2'] / (1 + u**pars['gamma']) - v
    dxdt=np.array([dudt,dvdt])
    return(dxdt)
    
def ModelB(x,t,pars):
    u=x[0]
    v=x[1]
    uterm=u/(1+pars['IPTG']/pars['K'])**pars['eta']
    dudt = pars['a1'] / (1 + v**pars['beta']) - u
    dvdt = pars['a2'] / (1 + uterm**pars['gamma']) - v
    dxdt=np.array([dudt,dvdt])
    return(dxdt)

