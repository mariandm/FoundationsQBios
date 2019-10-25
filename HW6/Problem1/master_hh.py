# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:39:28 2018

@author: Alexander Bo Lee
Based off of lab and code by Joshua Weitz
"""

from __future__ import division
import numpy as np
from scipy import integrate
from plot_hh import plot_hh
from model_hh_refrac import model_hh
from impulse_t_refrac import impulse_t

# master_hh
# A script for HH models
#
# This simulates the HH model using parameters embedded in pars
# a stimulus current, and a time range over which the simulation should run

pars={}
# Parameters, including on/off functions
pars['gKbar'] = 36    # mS/cm^2
pars['gNabar'] = 120  # mS/cm^2
pars['gL'] = 0.3      # mS/cm^2
pars['EK'] = -12       # mV
pars['ENa'] = 120     # mV
pars['EL'] = 10.6     # mV
pars['C'] = 1         # muF/cm^2
pars['alphan'] = lambda V: 0.01*(10-V)/(np.exp(1-V/10)-1)
pars['betan'] = lambda V: 0.125*np.exp(-V/80)
pars['alpham'] = lambda V: 0.1*(25-V)/(np.exp(2.5-V/10)-1)
pars['betam'] = lambda V: 4*np.exp(-V/18)
pars['alphah'] = lambda V: 0.07*np.exp(-V/20)
pars['betah'] = lambda V: (np.exp(3-V/10)+1)**-1

# Initial conditions
#pars['V0'] = 0
#pars['n0'] = 0
#pars['m0'] = 0
#pars['h0'] = 0
pars['V0']=0
pars['n0']=pars['alphan'](pars['V0'])/(pars['alphan'](pars['V0'])+pars['betan'](pars['V0']))
pars['m0']=pars['alpham'](pars['V0'])/(pars['alpham'](pars['V0'])+pars['betam'](pars['V0']))
pars['h0']=pars['alphah'](pars['V0'])/(pars['alphah'](pars['V0'])+pars['alphah'](pars['V0']))

pars['tryI']=25
pars['RefracTime']=12.07

# Run the model
t0 = 0
tf = 40
tstep = 0.02
t = np.arange(t0,tf+tstep,tstep)
y0 = np.array([pars['V0'], pars['n0'], pars['m0'], pars['h0']])
y = integrate.odeint(model_hh,y0,t,args=(pars,))

# Store the results
# FILL IN WHEREVER YOU SEE ...
dyn = {}
dyn['t'] = t
dyn['V'] = y[:,0]
dyn['n'] = y[:,1]
dyn['m'] = y[:,2]
dyn['h'] = y[:,3]
dyn['gK'] = pars['gKbar']*dyn['n']**4
dyn['gNa'] = pars['gNabar']*dyn['m']**3*dyn['h']
dyn['gL'] = pars['gL']*np.ones(len(dyn['t']))
dyn['IK'] = pars['gKbar'] * dyn['n']**4 * (dyn['V']-pars['EK'])
dyn['INa'] = pars['gNabar'] * dyn['m']**3 * dyn['h'] * (dyn['V']-pars['ENa'])
dyn['IL'] = pars['gL'] * (dyn['V']-pars['EL'])
dyn['appliedI'] = impulse_t(dyn['t'],pars['tryI'],pars['RefracTime'])

plot_hh(dyn,pars,savePDF=True,figName='Refractory1207')