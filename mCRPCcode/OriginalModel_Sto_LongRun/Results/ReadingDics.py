#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:18:44 2019

@author: mariandm
"""

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter 

s = open('Run_forced_Adaptive_2.txt', 'r').read()
lala=eval(s)
Adaptive=lala

s = open('Run_forced_Metronomic_2.txt', 'r').read()
lala=eval(s)
Metronomic=lala

s = open('Run_forced_MTD_2.txt', 'r').read()
lala=eval(s)
MTD=lala

s = open('Run_forced_NoTreatment_2.txt', 'r').read()
lala=eval(s)
NoTreatment=lala


alldics=[MTD,Metronomic, Adaptive, NoTreatment]

labels=["MTD","Metronomic", "Adaptive", "NoTreatment"]

treatmentColors=["firebrick","pink","yellowgreen","teal","lightblue","mediumpurple","gold","palevioletred", "slategray","orange"]



allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(alldics[j][i][2]["count"]/50)
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="Greens")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label=)
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Proportion of experiments in which\n T- dominates the final population")
plt.tight_layout()
plt.savefig("ProportionExperimentsTminusColonization.pdf")
plt.show()

######################

allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(alldics[j][i][2]["Time"])
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="bone")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label=)
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Mean Time to T- colonization")
plt.tight_layout()
plt.savefig("MeanTimeToTminusColonization.pdf")
plt.show()

################################

allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        a=np.sum(alldics[j][i][0]["densities"])*alldics[j][i][0]["count"]/50
        b=np.sum(alldics[j][i][1]["densities"])*alldics[j][i][1]["count"]/50
        c=np.sum(alldics[j][i][2]["densities"])*alldics[j][i][2]["count"]/50
        tmp.append(a+b+c)
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="Greys")
plt.yticks(range(0,4),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label=)
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Mean Final Population Density")
plt.tight_layout()
plt.savefig("MeanFinalPopulationDensity.pdf")
plt.show()

