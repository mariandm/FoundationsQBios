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
from scipy.stats import linregress

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



reps=50
totend = []

for p in np.arange(2,41,2):
    info['prey_density']=p
    b=np.pi*(predator['r']**2)*predator['k']*predator['f']
    
    for ens in np.arange(reps):
        print(ens)
        [t,numeaten] = ibm_functions.ibm_predation(info,predator,prey)
        if ens==0:
            total_eaten=np.zeros(len(t))
        total_eaten += np.cumsum(numeaten)
    total_eaten=total_eaten/reps
    plt.plot(t,total_eaten)
    plt.plot(t,b*p*t)

    plt.close()
    totend.append(total_eaten[-1])
    
plt.plot(np.arange(2,41,2),totend)
plt.plot(np.arange(2,41,2),b*np.arange(2,41,2)*info['tf'])
plt.xlabel("Prey density",fontsize=15)
plt.ylabel("Total eaten at time 25",fontsize=15)
plt.legend(["IBM data", "Theory"])
plt.title("Handling time = 0 seconds", fontsize=20)
plt.savefig("P2Figure1a.pdf")
plt.show()

plt.plot(totend, b*np.arange(2,41,2)*info['tf'])
plt.xlabel("Data",fontsize=15)
plt.ylabel("Theory",fontsize=15)
ln=linregress(totend, b*np.arange(2,41,2)*info['tf'])
plt.legend(["Slope = "+ str(ln[0])])
plt.savefig("P2Figure1b.pdf")
plt.show()


predator['handling_time']=1.0 #sec


reps=50
totend = []

b=np.pi*(predator['r']**2)*predator['k']*predator['f']
for p in np.arange(2,41,2):
    info['prey_density']=p
    
    
    for ens in np.arange(reps):
        print(ens)
        [t,numeaten] = ibm_functions.ibm_predation(info,predator,prey)
        if ens==0:
            total_eaten=np.zeros(len(t))
        total_eaten += np.cumsum(numeaten)
    total_eaten=total_eaten/reps
    plt.plot(t,total_eaten)
    plt.plot(t,b*p*t/(1+a*p))

    plt.close()
    totend.append(total_eaten[-1])
    
a=b*predator['handling_time']
plt.plot(np.arange(2,41,2),totend)
plt.plot(np.arange(2,41,2),b*np.arange(2,41,2)*info['tf']/(1+a*np.arange(2,41,2)))
plt.xlabel("Prey density",fontsize=15)
plt.ylabel("Total eaten at time 25",fontsize=15)
plt.legend(["IBM data", "Theory"])
plt.title("Handling time = 1 second",fontsize=20)
plt.savefig("P2Figure2a.pdf")
plt.show()


#Let's assume that the fitness of a thing is det

plt.plot(totend, b*np.arange(2,41,2)*info['tf']/(1+a*np.arange(2,41,2)))
plt.xlabel("Data",fontsize=15)
plt.ylabel("Theory",fontsize=15)
ln=linregress(totend, b*np.arange(2,41,2)*info['tf']/(1+a*np.arange(2,41,2)))
plt.legend(["Slope = "+ str(ln[0])])
plt.savefig("P2Figure2b.pdf")
plt.show()
