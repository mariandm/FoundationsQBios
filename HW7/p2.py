#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 20:48:21 2019

@author: mariandm
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 18:47:19 2019

@author: Alex Bo Lee

An individual-based model of predators eating prey
Translated from Joshua Weitz
"""

import numpy as np
import ibm_functions
import matplotlib.pyplot as plt

#Basic information
info={}
info['prey_density']=10 # density/cm^2
info['maxX']=10 # cm
info['maxY']=10 # cm
info['tf']=25 # sec
info['dt']=0.1 #time step
info['replenish_prey']=1 #should prey regenerate?
info['viz_dyn']=0 #1 for animation, 0 for no animation

# Define the prey
prey = {}
prey['num'] = info['prey_density']*info['maxX']*info['maxY']
prey['pos']=np.random.uniform(size=(prey['num'],2))
prey['pos'][:,0]=prey['pos'][:,0]*info['maxX']
prey['pos'][:,1]=prey['pos'][:,1]*info['maxY']
prey['diffusion']=0.5 #cm^2/sec

#Place the predator
predator={}
predator['pos'] = np.array( [[ info['maxX']/2, info['maxY']/2 ]] )
predator['theta']=np.random.uniform()*2*np.pi #Angle of movement
predator['r']=1.25 #Radius
predator['k']=0.3 #Detection
predator['f']=0.3 #successful capture per time
predator['vel']=.1 #cm/sec
predator['handling_time']=0.0 #sec
predator['tau']=6 #run time length, sec
predator['trun'] = 0 # Initialize each run



reps=20
totend1 = []

for k in np.arange(0.1,1.1,0.1):
    predator['k']=k
    b=np.pi*(predator['r']**2)*info['prey_density']*predator['f']
    
    for ens in np.arange(reps):
        print(ens)
        [t,numeaten] = ibm_functions.ibm_predation(info,predator,prey)
        if ens==0:
            total_eaten=np.zeros(len(t))
        total_eaten += np.cumsum(numeaten)
    total_eaten=total_eaten/reps
    #plt.plot(t,total_eaten)
    #plt.plot(t,b*p*t)

    #plt.close()
    totend1.append(total_eaten[-1])
    

#plt.plot(np.arange(0.3,1,0.1),b*np.arange(0.3,1,0.1)*info['tf'])



predator['k']=0.3

totend = []

for f in np.arange(0.1,1.1,0.1):
    predator['f']=f
    b=np.pi*(predator['r']**2)*info['prey_density']*predator['k']
    
    for ens in np.arange(reps):
        print(ens)
        [t,numeaten] = ibm_functions.ibm_predation(info,predator,prey)
        if ens==0:
            total_eaten=np.zeros(len(t))
        total_eaten += np.cumsum(numeaten)
    total_eaten=total_eaten/reps
    #plt.plot(t,total_eaten)
    #plt.plot(t,b*p*t)

    #plt.close()
    totend.append(total_eaten[-1])
    
b=np.arange(0.1,1.1,0.1)*info['prey_density']*0.3*np.pi*(predator['r']**2)*info['prey_density']
plt.plot(np.arange(0.1,1.1,0.1),totend1,color="palevioletred")
plt.plot(np.arange(0.1,1.1,0.1),totend,color="lightseagreen")
plt.plot(np.arange(0.1,1.1,0.1),b)
plt.xlabel("Variable",fontsize=15)
plt.ylabel("Total eaten at time 25 s",fontsize=15)
plt.legend(["Detection", "Capture","Theory"])
plt.title("Effects of detection and capture", fontsize=20)
plt.savefig("P2Figure3a.pdf")
plt.show()

