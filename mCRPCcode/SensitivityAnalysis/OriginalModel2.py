#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:18:27 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

ESS = np.array([[3196.34703196346,7077.62557077625,2739.72602739727,26027.3972602740],[5607.47663551401,7476.63551401869,467.289719626177,27102.8037383178],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[2786.37770897833,6501.54798761610,3405.57275541796,25386.9969040248],[6270.62706270627,7260.72607260726,330.033003300335,27722.7722772277],[5140.18691588785,7476.63551401869,934.579439252341,27102.8037383178],[6060.60606060606,7575.75757575758,1.e-09,27272.7272727273],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[7141.90425687459,7142.65465924746,0.833764523678801,28570.7853484204],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[4310.34482758620,6034.48275862069,2586.20689655173,25862.0689655172],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484]])
matrixCoefficients = np.array([[0.7,0.9,0.4,0.6,0.5,0.8],[0.7,0.8,0.4,0.6,0.5,0.9],[0.6,0.9,0.4,0.8,0.5,0.7],[0.6,0.9,0.4,0.7,0.5,0.8],[0.6,0.8,0.4,0.7,0.5,0.9],[0.7,0.9,0.4,0.5,0.6,0.8],[0.7,0.8,0.4,0.5,0.6,0.9],[0.6,0.9,0.4,0.5,0.7,0.8],[0.6,0.8,0.4,0.5,0.7,0.9],[0.6,0.7,0.4,0.5,0.8,0.9],[0.5,0.9,0.4,0.8,0.6,0.7],[0.5,0.9,0.4,0.7,0.6,0.8],[0.5,0.8,0.4,0.7,0.6,0.9],[0.5,0.9,0.4,0.6,0.7,0.8],[0.5,0.8,0.4,0.6,0.7,0.9],[0.5,0.7,0.4,0.6,0.8,0.9],[0.4,0.9,0.5,0.8,0.6,0.7],[0.4,0.9,0.5,0.7,0.6,0.8],[0.4,0.8,0.5,0.7,0.6,0.9],[0.4,0.9,0.5,0.6,0.7,0.8],[0.4,0.8,0.5,0.6,0.7,0.9],[0.4,0.7,0.5,0.6,0.8,0.9]])

# Set matrixIndex = 7 for Representative patient #1
# Set matrixInded = 5 for Representative patient #2
treatmentType = 'MTD'
setNumber=500

for curr_matrix in range(17,22):
    matrixIndex = curr_matrix


    alphas = matrixCoefficients[matrixIndex,:]
    
    
    
    finalDensity=np.zeros(setNumber)
    TminusProportion=np.zeros(setNumber)
    TminusColonizationTime=np.zeros(setNumber)
    
    
    
    for curr_set in range(setNumber):
        print(curr_set)
    
        #Creating sets
        
        #Original Parameters
        scale = .01
        r = np.array([0.27726, 0.34657, 0.66542])
        r = r*scale
        # PSA dynamics
        sigmaPSA = 0.5;
        
        
        k1 = [1.5, 10000, 10000]
        k2 = [0.5, 100, 10000]
        
        #randomizing the parameters
        random=np.random.uniform(-.1,.1,10)
        r=[r[0]+r[0]*random[0],r[1]+r[1]*random[1],r[2]+r[2]*random[2]]
        sigmaPSA=sigmaPSA+sigmaPSA*random[3]
        k1=[k1[0]+k1[0]*random[4],k1[1]+k1[1]*random[5],k1[2]+k1[2]*random[6]]
        k2=[k2[0]+k2[0]*random[7],k2[1]+k2[1]*random[8],k2[2]+k2[2]*random[9]]
        
        ########################################
        
        
        #Initial tumor densities set at 40% of ESS values
        y0 = ESS[matrixIndex, :]* 0.4;
        
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
            k = [y[1] * k1[0], k1[1], k1[2]]
            
            # T+, TP, T-, and PSA ODE's
            dydt = np.zeros([4]);
            
            dydt[0] = y[0] * r[0] * (1 - ( ( y[0] + alphas[0] * y[1] + alphas[1] * y[2] ) / k[0] ) )
            dydt[1] = y[1] * r[1] * (1 - ( ( alphas[2] * y[0] + y[1] + alphas[3] * y[2] ) / k[1] ) )
            dydt[2] = y[2] * r[2] * (1 - ( ( alphas[4] * y[0] + alphas[5] * y[1] + y[2] ) / k[2] ) )
            dydt[3] = sum(y[0:3]) - sigmaPSA * y[3]
        
            y = y + dydt;
        
        
        # No abiraterone treatment sets Abi flag to 0 for entire simulation (and this is used to zero out for all other treatment types).
        AbiOnOffFlag = np.zeros([maxSimulationTime])
        
        
        # SOC full dose abiraterone set Abi flag to 1 after first treatment time
        if treatmentType=='MTD':
            AbiOnOffFlag[firstTreatmentTime:maxSimulationTime] = 1;
        
        
        
        # Metronomic gives 800 long induction with cycles thereafter.
        if treatmentType=='Metronomic':
            
            #Induction period
            AbiOnOffFlag[firstTreatmentTime:firstTreatmentTime + 800] = 1
            
            # Cycles start 2000 time units after the first treatment time, are 1000
            # time units apart, and administer Abiraterone for 200 time units.
            for cycles in range(0,10):
                AbiOnOffFlag[firstTreatmentTime + 2000 + (cycles * 1000):firstTreatmentTime + 2000 + (cycles * 1000) + 200] = 1
        
        
        
        # Run the ODE
        # -----------
        
        # Set initial state.
        c=0
        y = copy.deepcopy(y0)
        
        # Create and initialize matrix for ODE solution
        allSolution = []
        allSolution.append(list(y))
        count=1
        #Run ODE with treatment selection.
        for time in range(1,maxSimulationTime):
            count=count+1
            # Adaptive Abi is built during the simulation. Turns Abi on once the 
            # PSA zenith value is reached and turns it off once the nadir is reached.
            if treatmentType=='Adaptive':
                if y[3] > PSA_zenith:
                    AbiOnOffFlag[time] = 1
                elif (y[3] < PSA_nadir):
                    AbiOnOffFlag[time] = 0
                else:
                    AbiOnOffFlag[time] = AbiOnOffFlag[time - 1] 
            
            # If Abi is being given, then use Abi parameters.
            if (AbiOnOffFlag[time] == 1):
                k = [y[1] * k2[0], k2[1], k2[2]]
            
            # If Abi is not being given, use naive parameters.
            elif AbiOnOffFlag[time] == 0:
                k = [y[1] * k1[0], k1[1], k1[2]]
            
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
            if(c==0):
                a=count
                c=1
            
            # Append to solutions vector. 
            allSolution.append(list(y))
            
        
        
        allSolution=np.array(allSolution)
        
        finalDensity[curr_set]=np.sum(y[0:3])
        TminusProportion[curr_set]=y[2]/finalDensity[curr_set]
        print(TminusProportion[curr_set])
        #if(TminusProportion[curr_set]<0.2):
         #   print(r)
         #   print(k1)
         #   print(k2)
         #   print(sigmaPSA)
        #TminusColonizationTime.append()
    
    #plt.plot(range(maxSimulationTime),allSolution[:,3])
    #plt.xlim([0,10000])
    #plt.show()
    
    #plt.plot(range(maxSimulationTime),allSolution[:,0],color="b",label="T+")
    #plt.plot(range(maxSimulationTime),allSolution[:,1],color="r",label="TP")
    #plt.plot(range(maxSimulationTime),allSolution[:,2],color="g",label="T-")
    #plt.ylim([-100,15000])
    #plt.xlim([0,10000])
    #plt.legend()
    
    
    #plt.hist(finalDensity,color="dimgrey")
    #plt.xlabel("Population Density")
    #plt.ylabel("Frequency")
    #plt.xlim([8000,12000])
    #plt.title("Final Population Density, Treatment: "+treatmentType+", Patient "+str(curr_matrix+1))
    #plt.savefig("Sensitivity_FinalPopulationDensity_Treatment"+treatmentType+"_matrixIndex"+str(curr_matrix)+".pdf")
    #plt.close()
    
    plt.hist(TminusProportion,color="silver",bins=1,width=0.1)
    plt.xlabel("Proportion of T- cells")
    plt.ylabel("Frequency")
    plt.xlim([-0.1,1.1])
    plt.title("Final Proportion of T- cells, Treatment: "+treatmentType+", Patient "+str(curr_matrix+1))
    plt.savefig("Sensitivity_FinalProportion_Treatment"+treatmentType+"_matrixIndex"+str(curr_matrix)+".pdf")
    plt.close()    
    