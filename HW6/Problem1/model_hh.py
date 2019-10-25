# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:52:49 2018

@author: Alexander Bo Lee
"""

from impulse_t import impulse_t
import numpy as np

# model_hh
def model_hh(y,t,pars):
    '''
    function dydt = model_hh(y,t,pars)
    This simulates the HH model using parameters embedded in pars
    a stimulus current in pars['Idrive']
    '''

    #Pulse
    #units of muA/cm^2
    #Equivalent to mSmV/cm^2
    
    # Variables
    V = y[0]
    n = y[1]
    m = y[2]
    h = y[3]
    
    # Impulses
    I = impulse_t(t,pars['tryI']) # Specified in a function
    
    # Dynamics
    Vdot = (1/pars['C'])*(I-pars['gKbar']*n**4*(V-pars['EK'])-pars['gNabar']*m**3*h*(V-pars['ENa'])-pars['gL']*(V-pars['EL']))
    # To avoid some issues, do the following:
    if V != 10:
        ndot = pars['alphan'](V)*(1-n)-pars['betan'](V)*n
    else:
        ndot = 0.1*(1-n)-pars['betan'](V)*n
    if V != 25:
        mdot = pars['alpham'](V)*(1-m)-pars['betam'](V)*m
    else:
        mdot = 1*(1-m)-pars['betam'](V)*m
    hdot = pars['alphah'](V)*(1-h)-pars['betah'](V)*h
    
    # Fill in
    dydt = np.array([Vdot, ndot, mdot, hdot])
    return dydt