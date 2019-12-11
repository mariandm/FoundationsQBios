#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:18:27 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

ESS = 1e4*np.array([[1.1111,0.5555,0.0000,3.3333],[1.1111,0.5556,0.0000,3.3333],[1.1157,0.4545,0.1240,3.3884],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.1111,0.5556,0.0000,3.3333],[1.1111,0.5556,0.0000,3.3333],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750]])
matrixCoefficients = np.array([[0.7,0.9,0.4,0.6,0.5,0.8],[0.7,0.8,0.4,0.6,0.5,0.9],[0.6,0.9,0.4,0.8,0.5,0.7],[0.6,0.9,0.4,0.7,0.5,0.8],[0.6,0.8,0.4,0.7,0.5,0.9],[0.7,0.9,0.4,0.5,0.6,0.8],[0.7,0.8,0.4,0.5,0.6,0.9],[0.6,0.9,0.4,0.5,0.7,0.8],[0.6,0.8,0.4,0.5,0.7,0.9],[0.6,0.7,0.4,0.5,0.8,0.9],[0.5,0.9,0.4,0.8,0.6,0.7],[0.5,0.9,0.4,0.7,0.6,0.8],[0.5,0.8,0.4,0.7,0.6,0.9],[0.5,0.9,0.4,0.6,0.7,0.8],[0.5,0.8,0.4,0.6,0.7,0.9],[0.5,0.7,0.4,0.6,0.8,0.9],[0.4,0.9,0.5,0.8,0.6,0.7],[0.4,0.9,0.5,0.7,0.6,0.8],[0.4,0.8,0.5,0.7,0.6,0.9],[0.4,0.9,0.5,0.6,0.7,0.8],[0.4,0.8,0.5,0.6,0.7,0.9],[0.4,0.7,0.5,0.6,0.8,0.9]])

# Set matrixIndex = 7 for Representative patient #1
# Set matrixInded = 5 for Representative patient #2
matrixIndex = 6


treatmentType = 'Adaptive'
subtreatmentType=2
alphas = matrixCoefficients[matrixIndex,:]

#Why is there a "scale factor"
scale = .01
r = np.array([0.27726, 0.34657, 0.66542])
r = r*scale

# PSA dynamics
sigmaPSA = 0.5;

#Initial tumor densities set at 40% of ESS values
y0 = ESS[matrixIndex, :]* 0.4;
y0[y0<1e-9] = 1e-9

# Set simulation time.
maxSimulationTime = 10000

#Give abiraterone at what % of ESS PSA?
maxPSAPercent = 0.8;
PSA_zenith = ESS[matrixIndex,3] * maxPSAPercent
PSA_nadir = PSA_zenith * maxPSAPercent/2


#AbiOnOffFlag = np.zeros([maxSimulationTime])


# Find first treatment time to build treatments.
y = copy.deepcopy(y0)
firstTreatmentTime = 0
while y[3] < PSA_zenith:
    
    firstTreatmentTime = firstTreatmentTime + 1;
    
    # Update carrying capacities with current symbiotic T+.
    k = [15000, 10000, 10000]
    
    # T+, TP, T-, and PSA ODE's
    dydt = np.zeros([4]);
    
    dydt[0] = y[0] * r[0] * (1 - ( ( y[0] + alphas[0] * y[1] + alphas[1] * y[2] ) / k[0] ) )
    dydt[1] = y[1] * r[1] * (1 - ( ( alphas[2] * y[0] + y[1] + alphas[3] * y[2] ) / k[1] ) )
    dydt[2] = y[2] * r[2] * (1 - ( ( alphas[4] * y[0] + alphas[5] * y[1] + y[2] ) / k[2] ) )
    dydt[3] = sum(y[0:3]) - sigmaPSA * y[3]

    y = y + dydt;


# No abiraterone treatment sets Abi flag to 0 for entire simulation (and this is used to zero out for all other treatment types).
AbiOnOffFlag = np.zeros([maxSimulationTime])
ADTOnOffFlag = np.zeros([maxSimulationTime])


if treatmentType=='Adaptive':
    if subtreatmentType==2:
        ADTOnOffFlag[firstTreatmentTime:firstTreatmentTime + 800] = 1
        AbiOnOffFlag[firstTreatmentTime:firstTreatmentTime + 800] = 0

# SOC full dose abiraterone set Abi flag to 1 after first treatment time
if treatmentType=='MTD':
    if subtreatmentType==0:
        AbiOnOffFlag[firstTreatmentTime:maxSimulationTime] = 0
        ADTOnOffFlag[firstTreatmentTime:maxSimulationTime] = 1
    elif subtreatmentType==1:
        AbiOnOffFlag[firstTreatmentTime:maxSimulationTime] = 1
        ADTOnOffFlag[firstTreatmentTime:maxSimulationTime] = 1



# Metronomic gives 800 long induction with cycles thereafter.
if treatmentType=='Metronomic': 
    
    #Induction period
    ADTOnOffFlag[firstTreatmentTime:firstTreatmentTime + 800] = 1
    
    if subtreatmentType==0: #ADT+Abi, none
    # Cycles start 2000 time units after the first treatment time, are 1000
    # time units apart, and administer Abiraterone for 200 time units.
        for cycles in range(0,10):
            AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1

    elif subtreatmentType==1: #ADT+Abi, none, ADT, none
    # Cycles start 2000 time units after the first treatment time, are 1000
    # time units apart, and administer Abiraterone for 200 time units.
        for cycles in range(0,10,2):
            AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1
        for cycles in range(1,10,2):
            AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 0
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1

    elif subtreatmentType==2: #ADT+abi, ADT, none
    # Cycles start 2000 time units after the first treatment time, are 1000
    # time units apart, and administer Abiraterone for 200 time units.
        for cycles in range(0,10):
            AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200):firstTreatmentTime + 2000 + (cycles * 1200) + 200] = 1
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200):firstTreatmentTime + 2000 + (cycles * 1200) + 200] = 1

            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200) + 200:firstTreatmentTime + 2000 + (cycles * 1200) + 400] = 1

    elif subtreatmentType==3: #ADT, ADT+abi, none
    # Cycles start 2000 time units after the first treatment time, are 1000
    # time units apart, and administer Abiraterone for 200 time units.
        for cycles in range(0,10):
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200):firstTreatmentTime + 2000 + (cycles * 1200) + 200] = 1

            AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200) + 200:firstTreatmentTime + 2000 + (cycles * 1200) + 400] = 1
            ADTOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1200) + 200:firstTreatmentTime + 2000 + (cycles * 1200) + 400] = 1



# Run the ODE
# -----------

# Set initial state.
c=0
y = copy.deepcopy(y0)

# Create and initialize matrix for ODE solution
allSolution = []
allSolution.append(list(y))
count=True
warning=False
#Run ODE with treatment selection.
for time in range(1,maxSimulationTime):
    # Adaptive Abi is built during the simulation. Turns Abi on once the 
    # PSA zenith value is reached and turns it off once the nadir is reached.
    if treatmentType=='Adaptive':
        if subtreatmentType==0:
            if y[3] > PSA_zenith:
                AbiOnOffFlag[time] = 1
                ADTOnOffFlag[time] = 1
            elif (y[3] < PSA_nadir):
                AbiOnOffFlag[time] = 0
                ADTOnOffFlag[time] = 0
            else:
                AbiOnOffFlag[time] = AbiOnOffFlag[time - 1] 
                ADTOnOffFlag[time] = AbiOnOffFlag[time - 1] 
        elif subtreatmentType==1:
            if y[3] > PSA_zenith:
                if(count):
                    AbiOnOffFlag[time] = 1
                    ADTOnOffFlag[time] = 1
                else:
                    AbiOnOffFlag[time] = 0
                    ADTOnOffFlag[time] = 1
            elif (y[3] < PSA_nadir):
                AbiOnOffFlag[time] = 0
                ADTOnOffFlag[time] = 0
            else:
                AbiOnOffFlag[time] = AbiOnOffFlag[time - 1] 
                ADTOnOffFlag[time] = AbiOnOffFlag[time - 1] 

            if ((ADTOnOffFlag[-2] - ADTOnOffFlag[-1])>0):
                count=not(count)
###############################################                
        elif subtreatmentType==2 and time>(firstTreatmentTime+2000):
            if y[3] > PSA_zenith:
                AbiOnOffFlag[time] = 1
                ADTOnOffFlag[time] = 1
            elif (y[3] < PSA_nadir):
                AbiOnOffFlag[time] = 0
                ADTOnOffFlag[time] = 0
            else:
                AbiOnOffFlag[time] = AbiOnOffFlag[time - 1] 
                ADTOnOffFlag[time] = AbiOnOffFlag[time - 1] 
#################################################                

    if treatmentType=='MTD' and subtreatmentType==2 and time>firstTreatmentTime:
        if y[3] < PSA_zenith:
            warning=True
        if(y[3]>PSA_zenith and warning):
            count=False
        if(count):
            AbiOnOffFlag[time]=0
            ADTOnOffFlag[time]=1
        else:
            AbiOnOffFlag[time]=1
            ADTOnOffFlag[time]=1
        
    
    # If Abi is being given, then use Abi parameters.
    if (ADTOnOffFlag[time] == 0):
        k = [15000, 10000, 10000]
    elif (AbiOnOffFlag[time] == 1):
        k = [y[1] * 0.5, 100, 10000]
    
    # If Abi is not being given, use naive parameters.
    elif AbiOnOffFlag[time] == 0:
        k = [y[1] * 1.5, 10000, 10000]
    
    # T+, TP, T-, and PSA ODE's
    dydt = np.zeros([4])
    
    dydt[0] = y[0] * r[0] * (1 - ( ( y[0] + alphas[0] * y[1] + alphas[1] * y[2] ) / k[0] ) )
    dydt[1] = y[1] * r[1] * (1 - ( ( alphas[2] * y[0] + y[1] + alphas[3] * y[2] ) / k[1] ) )
    dydt[2] = y[2] * r[2] * (1 - ( ( alphas[4] * y[0] + alphas[5] * y[1] + y[2] ) / k[2] ) )
    dydt[3] = sum(y[0:3]) - sigmaPSA * y[3]
    
    y = y + dydt;
        
    # Lower bound check to keep a small density of any cell type.
    y[np.where(y<1e-9)] = 1e-9;
    
    y[y<1e-9] = 1e-9

    
    # Append to solutions vector. 
    allSolution.append(list(y))
    


    
allSolution=np.array(allSolution)
plt.plot(range(maxSimulationTime),allSolution[:,3])
plt.xlim([0,10000])
plt.show()

plt.plot(range(maxSimulationTime),allSolution[:,0],color="b",label="T+")
plt.plot(range(maxSimulationTime),allSolution[:,1],color="r",label="TP")
plt.plot(range(maxSimulationTime),allSolution[:,2],color="g",label="T-")
plt.ylim([-100,15000])
plt.xlim([0,10000])
plt.legend()



#Weno, ya está la weonada esta en piton. Ahora vamo a convertirlo a estocastico vale chico.