#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:18:44 2019

@author: mariandm
"""

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter 

s = open('Run_forced_Adaptive_subtreat_0_2.txt', 'r').read()
lala=eval(s)
Adaptive0=lala["Adaptive"][0]

s = open('Run_forced_Adaptive_subtreat_1_2.txt', 'r').read()
lala=eval(s)
Adaptive1=lala["Adaptive"][1]

s = open('Run_forced_Adaptive_subtreat_2_2.txt', 'r').read()
lala=eval(s)
Adaptive2=lala["Adaptive"][2]

s = open('Run_forced_Metronomic_subtreat_0_2.txt', 'r').read()
lala=eval(s)
Metronomic0=lala["Metronomic"][0]

s = open('Run_forced_Metronomic_subtreat_1_2.txt', 'r').read()
lala=eval(s)
Metronomic1=lala["Metronomic"][1]

s = open('Run_forced_Metronomic_subtreat_2_2.txt', 'r').read()
lala=eval(s)
Metronomic2=lala["Metronomic"][2]

s = open('Run_forced_Metronomic_subtreat_3_2.txt', 'r').read()
lala=eval(s)
Metronomic3=lala["Metronomic"][3]

s = open('Run_forced_MTD_subtreat_0_2.txt', 'r').read()
lala=eval(s)
MTD0=lala["MTD"][0]

s = open('Run_forced_MTD_subtreat_1_2.txt', 'r').read()
lala=eval(s)
MTD1=lala["MTD"][1]

s = open('Run_forced_MTD_subtreat_2_2.txt', 'r').read()
lala=eval(s)
MTD2=lala["MTD"][2]

s = open('Run_forced_NoTreatment_subtreat_0_2.txt', 'r').read()
lala=eval(s)
NoTreatment0=lala["NoTreatment"][0]

alldics=[MTD0,MTD1,MTD2,Metronomic0,Metronomic1,Metronomic2,Metronomic3, Adaptive0,Adaptive1,Adaptive2, NoTreatment0]

labels=["MTD 1","MTD 2","MTD 3","Metronomic 1", "Metronomic 2", "Metronomic 3", "Metronomic 4", "Adaptive 1","Adaptive 2","Adaptive 3", "No Treatment"]

treatmentColors=["firebrick","pink","yellowgreen","teal","lightblue","mediumpurple","gold","palevioletred", "slategray","orange"]



allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(alldics[j][i][2]["count"]/50)
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues,cmap="Greens")
plt.yticks(range(0,11),tuple(labels))
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
plt.yticks(range(0,11),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label=)
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Mean Time for T- colonization")
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
plt.yticks(range(0,11),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,fraction=0.02, pad=0.04)#,label=)
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.title("Mean Final Population Density")
plt.tight_layout()
plt.savefig("MeanFinalPopulationDensity.pdf")
plt.show()

