#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 18:07:15 2019

@author: mariandm
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import copy

start_time = time.time()

ESS = np.array([[3196.34703196346,7077.62557077625,2739.72602739727,26027.3972602740],[5607.47663551401,7476.63551401869,467.289719626177,27102.8037383178],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[2786.37770897833,6501.54798761610,3405.57275541796,25386.9969040248],[6270.62706270627,7260.72607260726,330.033003300335,27722.7722772277],[5140.18691588785,7476.63551401869,934.579439252341,27102.8037383178],[6060.60606060606,7575.75757575758,1.e-09,27272.7272727273],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[6617.64705882353,7352.94117647059,1.e-09,27941.1764705882],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[7141.90425687459,7142.65465924746,0.833764523678801,28570.7853484204],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[7142.85714285714,7142.85714285714,1.e-09,28571.4285714286],[1.e-09,4545.45454545455,6818.18181818182,22727.2727272727],[4310.34482758620,6034.48275862069,2586.20689655173,25862.0689655172],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484],[7096.77419354839,6451.61290322581,1.e-09,27096.7741935484]])
matrixCoefficients = np.array([[0.7,0.9,0.4,0.6,0.5,0.8],[0.7,0.8,0.4,0.6,0.5,0.9],[0.6,0.9,0.4,0.8,0.5,0.7],[0.6,0.9,0.4,0.7,0.5,0.8],[0.6,0.8,0.4,0.7,0.5,0.9],[0.7,0.9,0.4,0.5,0.6,0.8],[0.7,0.8,0.4,0.5,0.6,0.9],[0.6,0.9,0.4,0.5,0.7,0.8],[0.6,0.8,0.4,0.5,0.7,0.9],[0.6,0.7,0.4,0.5,0.8,0.9],[0.5,0.9,0.4,0.8,0.6,0.7],[0.5,0.9,0.4,0.7,0.6,0.8],[0.5,0.8,0.4,0.7,0.6,0.9],[0.5,0.9,0.4,0.6,0.7,0.8],[0.5,0.8,0.4,0.6,0.7,0.9],[0.5,0.7,0.4,0.6,0.8,0.9],[0.4,0.9,0.5,0.8,0.6,0.7],[0.4,0.9,0.5,0.7,0.6,0.8],[0.4,0.8,0.5,0.7,0.6,0.9],[0.4,0.9,0.5,0.6,0.7,0.8],[0.4,0.8,0.5,0.6,0.7,0.9],[0.4,0.7,0.5,0.6,0.8,0.9]])

# Set matrixIndex = 7 for Representative patient #1
# Set matrixInded = 5 for Representative patient #2

treatmentType = 'MTD'

#Why is there a "scale factor"
scale = .01
r = np.array([0.27726, 0.34657, 0.66542])
r = r*scale

# PSA dynamics
sigmaPSA = 0.5;



# Set simulation time.
maxSimulationTime = 10000
replicateNumber = 50

dicti=dict.fromkeys(range(22))
for k in range(22):
    dicti[k]=dict.fromkeys(range(3))
    for l in range(3):
        dicti[k][l]=dict.fromkeys(["count","Time","PSA"],0)                
        dicti[k][l]["densities"]=np.array([0.0,0.0,0.0])



for curr_matrixIndex in range(len(ESS)):
    matrixIndex = curr_matrixIndex
    alphas = matrixCoefficients[matrixIndex,:]
    
    #Initial tumor densities set at 40% of ESS values
    y0 = ESS[matrixIndex, :]* 0.4
    y0 = np.ceil(y0)
            
    if (y0[2]<1):y0[2]=1

    #Give abiraterone at what % of ESS PSA?
    maxPSAPercent = 0.8;
    PSA_zenith = ESS[matrixIndex,3] * maxPSAPercent
    PSA_nadir = PSA_zenith * maxPSAPercent/2

    for curr_replicate in range(1,replicateNumber+1):
        
        AbiOnOffFlag = [0]

        # Set initial state.
        y = copy.deepcopy(y0)

        # Create and initialize matrix for ODE solution
        allSolution = []
        allSolution.append(list(y))
        
        time=[0]
        firstTreatment=1
        while time[-1] < maxSimulationTime:
        
            if treatmentType=="MTD":
                if firstTreatment==1:
                    if y[3]<PSA_zenith:
                        AbiOnOffFlag.append(0)
                    else:
                        firstTreatment=0
                else:
                    AbiOnOffFlag.append(1)
            
            if treatmentType=="NoTreatment":
                AbiOnOffFlag.append(0)
            
            # Adaptive Abi is built during the simulation. Turns Abi on once the 
            # PSA zenith value is reached and turns it off once the nadir is reached.
            if treatmentType=="Metronomic":
                if firstTreatment==1:
                    if y[3]<PSA_zenith:
                        AbiOnOffFlag.append(0)
                    else:
                        firstTreatment=0
                        firstTreatmentTime=time[-1]
                        tot=firstTreatmentTime+2000
                        abo=[]
                        cutTimes=[]
                        while tot < maxSimulationTime:
                            tot+=200
                            abo.append(1)
                            cutTimes.append(tot)
                            if tot<maxSimulationTime:
                                tot+=1000
                                cutTimes.append(tot)
                                abo.append(0)
                            cutTimes=[firstTreatmentTime, firstTreatmentTime+800, firstTreatmentTime+2000]+cutTimes
                            abo=[0,1,0]+abo
                else:
                    AbiOnOffFlag.append(abo[np.where(time[-1]<np.array(cutTimes))[0][0]])
            
            
            if treatmentType=='Adaptive':
                if y[3] > PSA_zenith:
                    AbiOnOffFlag.append(1)
                elif (y[3] < PSA_nadir):
                    AbiOnOffFlag.append(0)
                else:
                    AbiOnOffFlag.append(AbiOnOffFlag[-1])
                
        
            # If Abi is being given, then use Abi parameters.
            if (AbiOnOffFlag[-1] == 1):
                k = [y[1] * 0.5, 100, 10000]
            
            # If Abi is not being given, use naive parameters.
            elif AbiOnOffFlag[-1] == 0:
                k = [y[1] * 1.5, 10000, 10000]
            
                
            y[3]= y[3] + sum(y[0:3]) - sigmaPSA * y[3]
            
            dydt = np.zeros([6])
            
            #T+ growth
            dydt[0] = y[0] * r[0] 
            #T+ death
            if k[0]== 0:
                dydt[1]=0
            else:
                dydt[1] = y[0] * r[0] * ( ( y[0] + alphas[0] * y[1] + alphas[1] * y[2] ) / k[0] )
            #TP growth
            dydt[2] = y[1] * r[1] 
            #TP death
            dydt[3] = y[1] * r[1] * ( ( alphas[2] * y[0] + y[1] + alphas[3] * y[2] ) / k[1] )
            #T- growth
            dydt[4] = y[2] * r[2] 
            #T- death
            dydt[5] = y[2] * r[2] * ( ( alphas[4] * y[0] + alphas[5] * y[1] + y[2] ) / k[2] )
        
            
            dt= -(1/sum(dydt))*np.log(np.random.uniform())
            time.append(time[-1] + dt)
            
            #reaction = np.argmin(dydt/sum(dydt)<np.random.uniform())
            reaction = np.where(np.random.uniform()<=np.cumsum(dydt)/sum(dydt))[0][0]
            
            if reaction==0:
                y[0]=y[0]+1
            elif reaction==1:
                y[0]=y[0]-1
            elif reaction==2:
                y[1]=y[1]+1
            elif reaction==3:
                y[1]=y[1]-1
            elif reaction==4:
                y[2]=y[2]+1
            else:
                y[2]=y[2]-1
        
########################################################################            
            #if (y[2]<1):y[2]=1
########################################################################
            allSolution.append(list(y))
            
        allSolution=np.array(allSolution)
        
        
        if curr_replicate%5==0:
        
            start=[]
            end=[]
            a = np.array(AbiOnOffFlag) 
                                
            for i in range(1,len(a)):
                if(not(a[i-1]) and a[i]):
                    start.append(time[i])
                if(a[i-1] and not(a[i])):
                    end.append(time[i])
            if len(start)>len(end):
                end.append(time[-1])
            elif len(end)>len(start):
                start=[0]+list(start)
                                    
            if treatmentType=="Adaptive":
                start=np.array(start)-20
                end=np.array(end)+20
            
            #print("--- %s seconds ---" % (time.time() - start_time))
            plt.plot(time,allSolution[:,3]/y0[3])
            plt.plot(time,np.ones(len(time))*PSA_zenith/y0[3],color="r",linestyle="--")
            plt.plot(time,np.ones(len(time))*PSA_nadir/y0[3],color="r",linestyle="--")
            
            for i in range(len(start)):
                plt.axvspan(start[i],end[i], facecolor='lightgray', alpha=0.5)
            
            plt.savefig("./LongRunOriginalStochastic/PSA_"+treatmentType+"_matrixIndex_"+str(matrixIndex)+"_replicate_"+str(curr_replicate)+".pdf")
            plt.close()
            #plt.show()
            
            plt.plot(time,allSolution[:,0],label="T+",color="b")
            plt.plot(time,allSolution[:,1],label="TP",color="r")
            plt.plot(time,allSolution[:,2],label="T-",color="g")
            
            for i in range(len(start)):
                plt.axvspan(start[i],end[i], facecolor='lightgray', alpha=0.5)
            
            #x=np.linspace(time[np.argmax(np.array(AbiOnOffFlag)>0)],time[np.argmin(np.array(AbiOnOffFlag)>0)])
            #y=np.zeros(len(x))+max(allSolution[:,0])
            #plt.fill_between(x, 0, y,color="lightgray")
            plt.ylim([-100,15000])
            plt.xlim([0,maxSimulationTime])
            #plt.legend(loc="upper right")
            plt.xlabel("Simulated Time")
            plt.ylabel("Population Density")
            #plt.show()
            plt.savefig("./LongRunOriginalStochastic/Densities_"+treatmentType+"_matrixIndex_"+str(matrixIndex)+"_replicate_"+str(curr_replicate)+".pdf")
            plt.close()
        
        dicti[matrixIndex][np.argmax(y[0:3])]["count"]+=1
        dicti[matrixIndex][np.argmax(y[0:3])]["PSA"]+=y[3]
        dicti[matrixIndex][np.argmax(y[0:3])]["densities"]+=y[0:3]
        if np.sum(np.logical_not(np.argmax(allSolution[:,0:3],1)==np.argmax(y[0:3])))>0:               
            dicti[matrixIndex][np.argmax(y[0:3])]["Time"]+=time[np.where(np.logical_not(np.argmax(allSolution[:,0:3],1)==np.argmax(y[0:3])))[0][-1]+1]
        
    for cellType in range(3):
        divide= dicti[matrixIndex][cellType]["count"]
        if divide>0:
            dicti[matrixIndex][cellType]["PSA"]/=divide
            dicti[matrixIndex][cellType]["densities"]/=divide
            dicti[matrixIndex][cellType]["Time"]/=divide

f = open("./LongRunOriginalStochastic/Run_"+treatmentType+".txt", "w")
f.write(str(dicti))
f.close()