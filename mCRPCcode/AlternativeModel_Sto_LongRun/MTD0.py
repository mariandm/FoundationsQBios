#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 16:15:13 2019

@author: mariandm
"""


import numpy as np
import matplotlib.pyplot as plt
import time
import copy


start_time = time.time()

ESS = 1e4*np.array([[1.1111,0.5555,0.0000,3.3333],[1.1111,0.5556,0.0000,3.3333],[1.1157,0.4545,0.1240,3.3884],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.1111,0.5556,0.0000,3.3333],[1.1111,0.5556,0.0000,3.3333],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.1842,0.5263,0.0000,3.4211],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.2500,0.5000,0.0000,3.5000],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750],[1.3750,0.3125,0.0000,3.3750]])
matrixCoefficients = np.array([[0.7,0.9,0.4,0.6,0.5,0.8],[0.7,0.8,0.4,0.6,0.5,0.9],[0.6,0.9,0.4,0.8,0.5,0.7],[0.6,0.9,0.4,0.7,0.5,0.8],[0.6,0.8,0.4,0.7,0.5,0.9],[0.7,0.9,0.4,0.5,0.6,0.8],[0.7,0.8,0.4,0.5,0.6,0.9],[0.6,0.9,0.4,0.5,0.7,0.8],[0.6,0.8,0.4,0.5,0.7,0.9],[0.6,0.7,0.4,0.5,0.8,0.9],[0.5,0.9,0.4,0.8,0.6,0.7],[0.5,0.9,0.4,0.7,0.6,0.8],[0.5,0.8,0.4,0.7,0.6,0.9],[0.5,0.9,0.4,0.6,0.7,0.8],[0.5,0.8,0.4,0.6,0.7,0.9],[0.5,0.7,0.4,0.6,0.8,0.9],[0.4,0.9,0.5,0.8,0.6,0.7],[0.4,0.9,0.5,0.7,0.6,0.8],[0.4,0.8,0.5,0.7,0.6,0.9],[0.4,0.9,0.5,0.6,0.7,0.8],[0.4,0.8,0.5,0.6,0.7,0.9],[0.4,0.7,0.5,0.6,0.8,0.9]])

# Set matrixIndex = 7 for Representative patient #1
# Set matrixInded = 5 for Representative patient #2

TreatmentsVector=["MTD"]
subtreatmentsVector=[[0]]


#Why is there a "scale factor"
scale = .01
r = np.array([0.27726, 0.34657, 0.66542])
r = r*scale

# PSA dynamics
sigmaPSA = 0.5;

# Set simulation time.
maxSimulationTime = 10000
replicateNumber = 50


dicti=dict.fromkeys(TreatmentsVector)
for i in range(len(TreatmentsVector)):
    dicti[TreatmentsVector[i]]=dict.fromkeys(subtreatmentsVector[i])
    for j in range(len(subtreatmentsVector[i])):
        dicti[TreatmentsVector[i]][subtreatmentsVector[i][j]]=dict.fromkeys(range(22))
        for k in range(22):
            dicti[TreatmentsVector[i]][subtreatmentsVector[i][j]][k]=dict.fromkeys(range(3))
            for l in range(3):
                dicti[TreatmentsVector[i]][subtreatmentsVector[i][j]][k][l]=dict.fromkeys(["count","Time","PSA"],0)                
                dicti[TreatmentsVector[i]][subtreatmentsVector[i][j]][k][l]["densities"]=np.array([0.0,0.0,0.0])



#############################we start all the fors here########################

for curr_treatmentType in range(len(TreatmentsVector)):
    treatmentType = TreatmentsVector[curr_treatmentType]
    
    for curr_subtreatmentType in subtreatmentsVector[curr_treatmentType]:
        subtreatmentType = curr_subtreatmentType

        for curr_matrixIndex in range(len(ESS)):
            matrixIndex = curr_matrixIndex
            
            alphas = matrixCoefficients[matrixIndex,:]
        
            #Initial tumor densities set at 40% of ESS values
            y0 = ESS[matrixIndex, :]* 0.4
            y0 = np.ceil(y0)
            
            if (y0[2]<1):y0[2]=1
            
            
            #Give abiraterone at what % of ESS PSA?
            maxPSAPercent = 0.8
            PSA_zenith = ESS[matrixIndex,3] * maxPSAPercent
            PSA_nadir = PSA_zenith * maxPSAPercent/2
        
            for curr_replicate in range(1,replicateNumber+1):
        
                AbiOnOffFlag = [0]
                ADTOnOffFlag = [0]
                
                
                # Set initial state.
                y = copy.deepcopy(y0)
                
                # Create and initialize matrix for ODE solution
                allSolution = []
                allSolution.append(list(y))
                
                time=[0]
                firstTreatment=1
                count=True
                warning=False
                
                while time[-1] < maxSimulationTime:
                
                    if treatmentType=="MTD":
                        if firstTreatment==1:
                            if y[3]<PSA_zenith:
                                AbiOnOffFlag.append(0)
                                ADTOnOffFlag.append(0)
                            else:
                                firstTreatment=0
                        else:
                            if subtreatmentType==0:
                                AbiOnOffFlag.append(0)
                                ADTOnOffFlag.append(1)
                                
                            elif subtreatmentType==1:
                                AbiOnOffFlag.append(1)
                                ADTOnOffFlag.append(1)
                                
                            elif subtreatmentType==2:
                                if y[3] < PSA_zenith:
                                    warning=True
                                if(y[3]>PSA_zenith and warning):
                                    count=False
                                if(count):
                                    AbiOnOffFlag.append(0)
                                    ADTOnOffFlag.append(1)
                                else:
                                    AbiOnOffFlag.append(1)
                                    ADTOnOffFlag.append(1)
                    
                    if treatmentType=="NoTreatment":
                        AbiOnOffFlag.append(0)
                        ADTOnOffFlag.append(0)
                    
                    # Adaptive Abi is built during the simulation. Turns Abi on once the 
                    # PSA zenith value is reached and turns it off once the nadir is reached.
                    if treatmentType=="Metronomic":
                        
                        if firstTreatment==1:
                            if y[3]<PSA_zenith:
                                AbiOnOffFlag.append(0)
                                ADTOnOffFlag.append(0)
                            else:
                                firstTreatment=0
                                firstTreatmentTime=time[-1]
                                tot=firstTreatmentTime
                                abo=[]
                                adto=[]
                                cutTimes=[]
                                
                                if subtreatmentType==0:#ADT+Abi, none
                                    while tot < maxSimulationTime:
                                        #ADT+abi
                                        tot+=200
                                        abo.append(1)
                                        adto.append(1)
                                        cutTimes.append(tot)
                                        #none
                                        tot+=1000
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(0)
                                        cutTimes=[firstTreatmentTime, firstTreatmentTime+800, firstTreatmentTime+2000]+cutTimes
                                        abo=[0,0,0]+abo
                                        adto=[0,1,0]+adto
                                        
                                elif subtreatmentType==1: #ADT+Abi, none, ADT, none
                                    while tot < maxSimulationTime:
                                        tot+=200
                                        abo.append(1)
                                        adto.append(1)
                                        cutTimes.append(tot)
                                        tot+=1000
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(0)
                                        tot+=200
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(1)
                                        tot+=1000
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(0)
                                        cutTimes=[firstTreatmentTime, firstTreatmentTime+800, firstTreatmentTime+2000]+cutTimes
                                        abo=[0,0,0]+abo
                                        adto=[0,1,0]+adto
                                        
                                elif subtreatmentType==2: #ADT+abi, ADT, none
                                    while tot < maxSimulationTime:
                                        tot+=200
                                        abo.append(1)
                                        adto.append(1)
                                        cutTimes.append(tot)
                                        tot+=200
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(1)
                                        tot+=1000
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(0)
                                        cutTimes=[firstTreatmentTime, firstTreatmentTime+800, firstTreatmentTime+2000]+cutTimes
                                        abo=[0,0,0]+abo
                                        adto=[0,1,0]+adto
                                        
                                elif subtreatmentType==3: #ADT, ADT+abi, none
                                    while tot < maxSimulationTime:
                                        tot+=200
                                        abo.append(0)
                                        adto.append(1)
                                        cutTimes.append(tot)
                                        tot+=200
                                        cutTimes.append(tot)
                                        abo.append(1)
                                        adto.append(1)
                                        tot+=1000
                                        cutTimes.append(tot)
                                        abo.append(0)
                                        adto.append(0)
                                        cutTimes=[firstTreatmentTime, firstTreatmentTime+800, firstTreatmentTime+2000]+cutTimes
                                        abo=[0,0,0]+abo
                                        adto=[0,1,0]+adto
                                    
                        else:
                            AbiOnOffFlag.append(abo[np.where(time[-1]<np.array(cutTimes))[0][0]])
                            ADTOnOffFlag.append(abo[np.where(time[-1]<np.array(cutTimes))[0][0]])
                    
                    if treatmentType=='Adaptive':
                        if subtreatmentType==0:
                            if y[3] > PSA_zenith:
                                AbiOnOffFlag.append(1)
                                ADTOnOffFlag.append(1)
                            elif (y[3] < PSA_nadir):
                                AbiOnOffFlag.append(0)
                                ADTOnOffFlag.append(0)
                            else:
                                AbiOnOffFlag.append(AbiOnOffFlag[-1])
                        elif subtreatmentType==1:
                            if y[3] > PSA_zenith:
                                if(count):
                                    AbiOnOffFlag.append(1)
                                    ADTOnOffFlag.append(1)
                                else:
                                    AbiOnOffFlag.append(0)
                                    ADTOnOffFlag.append(1)
                            elif (y[3] < PSA_nadir):
                                AbiOnOffFlag.append(0)
                                ADTOnOffFlag.append(0)
                            else:
                                AbiOnOffFlag.append(AbiOnOffFlag[-1])
                                ADTOnOffFlag.append(AbiOnOffFlag[-1])
                            
                            if ((ADTOnOffFlag[-2] - ADTOnOffFlag[-1])>0):
                                count=not(count)
                        
                
                    if (ADTOnOffFlag[-1] == 0):
                        k = [15000, 10000, 10000]
                    # If Abi is being given, then use Abi parameters.
                    elif (AbiOnOffFlag[-1] == 1):
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
                
                    
                    
                    allSolution.append(list(y))
                    
                allSolution=np.array(allSolution)
                
                if curr_replicate%10==0:
                
                    #print("--- %s seconds ---" % (time.time() - start_time))
                    plt.plot(time,allSolution[:,3]/y0[3])
                    plt.plot(time,np.ones(len(time))*PSA_zenith/y0[3],color="r",linestyle="--")
                    plt.plot(time,np.ones(len(time))*PSA_nadir/y0[3],color="r",linestyle="--")
                    plt.savefig("./tara/PSA_"+treatmentType+"_subtreat_"+str(subtreatmentType)+"_matrixIndex_"+str(matrixIndex)+"_replicate_"+str(curr_replicate)+".pdf")
                    plt.close()
                    
                    
                    plt.plot(time,allSolution[:,0],label="T+",color="b")
                    plt.plot(time,allSolution[:,1],label="TP",color="r")
                    plt.plot(time,allSolution[:,2],label="T-",color="g")
                    #x=np.linspace(time[np.argmax(np.array(AbiOnOffFlag)>0)],time[np.argmin(np.array(AbiOnOffFlag)>0)])
                    #y=np.zeros(len(x))+max(allSolution[:,0])
                    #plt.fill_between(x, 0, y,color="lightgray")
                    plt.ylim([-100,15000])
                    plt.xlim([0,maxSimulationTime])
                    plt.legend(loc="upper right")
                    plt.xlabel("Simulated Time")
                    plt.ylabel("Population Density")
                    plt.savefig("./tara/Densities_"+treatmentType+"_subtreat_"+str(subtreatmentType)+"_matrixIndex_"+str(matrixIndex)+"_replicate_"+str(curr_replicate)+".pdf")
                    plt.close()
                
                dicti[treatmentType][subtreatmentType][matrixIndex][np.argmax(y[0:3])]["count"]+=1
                dicti[treatmentType][subtreatmentType][matrixIndex][np.argmax(y[0:3])]["PSA"]+=y[3]
                dicti[treatmentType][subtreatmentType][matrixIndex][np.argmax(y[0:3])]["densities"]+=y[0:3]
                if np.sum(np.logical_not(np.argmax(allSolution[:,0:3],1)==np.argmax(y[0:3])))>0:               
                    dicti[treatmentType][subtreatmentType][matrixIndex][np.argmax(y[0:3])]["Time"]+=time[np.where(np.logical_not(np.argmax(allSolution[:,0:3],1)==np.argmax(y[0:3])))[0][-1]+1]
        
            for cellType in range(3):
                divide= dicti[treatmentType][subtreatmentType][matrixIndex][cellType]["count"]
                if divide>0:
                    dicti[treatmentType][subtreatmentType][matrixIndex][cellType]["PSA"]/=divide
                    dicti[treatmentType][subtreatmentType][matrixIndex][cellType]["densities"]/=divide
                    dicti[treatmentType][subtreatmentType][matrixIndex][cellType]["Time"]/=divide
     

f = open("./tara/Run_"+treatmentType+"_subtreat_"+str(subtreatmentType)+".txt", "w")
f.write(str(dicti))
f.close()
        
#f = open("demofile2.txt", "a")
#f.write(str(allSolution))
#f.close()

#f = open("demofile3.txt", "a")
#f.write(str(AbiOnOffFlag))
#f.close()


