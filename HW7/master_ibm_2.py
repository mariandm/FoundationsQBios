# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 18:47:19 2019

@author: Alex Bo Lee

An individual-based model of predators eating prey
Translated from Joshua Weitz
"""

import numpy as np
import ibm_functions

#Basic information
info={}
info['prey_density']=10 # density/cm^2
info['maxX']=10 # cm
info['maxY']=10 # cm
info['tf']=50 # sec
info['dt']=0.1 #time step
info['replenish_prey']=1 #should prey regenerate?
info['viz_dyn']=1 #1 for animation, 0 for no animation

# Define the prey
prey = {}
prey['num'] = info['prey_density']*info['maxX']*info['maxY']
prey['pos']=np.random.uniform(size=(prey['num'],2))
prey['pos'][:,0]=prey['pos'][:,0]*info['maxX']
prey['pos'][:,1]=prey['pos'][:,1]*info['maxY']
prey['diffusion']=0.005 #cm^2/sec

#Place the predator
predator={}
predator['pos'] = np.array( [[ info['maxX']/2, info['maxY']/2 ]] )
predator['theta']=np.random.uniform()*2*np.pi #Angle of movement
predator['r']=1.25 #Radius
predator['k']=0.3 #Detection
predator['f']=0.05 #successful capture per time
predator['vel']=.1 #cm/sec
predator['handling_time']=0.0 #sec
predator['tau']=6 #run time length, sec
predator['trun'] = 0 # Initialize each run

#simulate eating
reps=1
totend = np.zeros(reps)

for ens in np.arange(reps):
    print(ens)
    [t,numeaten] = ibm_functions.ibm_predation(info,predator,prey)
    total_eaten = np.cumsum(numeaten)
    totend[ens] = total_eaten[-1]