#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:18:44 2019

@author: mariandm
"""

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter 

s = open('Run_originalDet_2.txt', 'r').read()
lala=eval(s)

Adaptive=lala["Adaptive"]
Metronomic=lala["Metronomic"]
MTD=lala["MTD"]
NoTreatment=lala["NoTreatment"]


alldics=[MTD,Metronomic, Adaptive, NoTreatment]

labels=["MTD","Metronomic", "Adaptive", "NoTreatment"]

treatmentColors=["firebrick","pink","yellowgreen","teal","lightblue","mediumpurple","gold","palevioletred", "slategray","orange"]


#3 bar graphs for final population density
for i in [0,6,11]:#range(22):
    a=[]
    for j in range(len(alldics)):
        a.append(np.sum(alldics[j][i]["densities"]))
    x = np.arange(len(labels))
    plt.bar(x, a,color="gray")#, width)
    plt.xticks(range(4),labels)
    plt.xlabel("Treatment")
    plt.ylabel("Final Population Density")
    plt.ylim([0,15000])
    plt.show()
    
#Heatmap for T- percentage at final time
allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(alldics[j][i]["densities"][2]/np.sum(alldics[j][i]["densities"]))
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="Greens")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label="T- proportion")
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Proportion of T- cells at final time")
plt.tight_layout()
plt.savefig("TminusProportionFinalTime.pdf")
plt.show()

#Heatmap for T- percentage at final time
allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(np.sum(alldics[j][i]["densities"]))
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="Greys")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label="T- proportion")
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Population density at final time")
plt.tight_layout()
plt.savefig("PopulationDensityatFinalTime.pdf")
plt.show()



allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(np.sum(alldics[j][i]["Time"]))
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="bone")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label="T- proportion")
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Time to T- Colonization")
plt.tight_layout()
plt.savefig("TimetoTminusColonization.pdf")
plt.show()