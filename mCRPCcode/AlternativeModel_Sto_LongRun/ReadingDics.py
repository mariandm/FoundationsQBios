#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 13:18:44 2019

@author: mariandm
"""

import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter 

s = open('Run_Adaptive_subtreat_0_2.txt', 'r').read()
lala=eval(s)
Adaptive0=lala["Adaptive"][0]
s = open('Run_Adaptive_subtreat_1_2.txt', 'r').read()
lala=eval(s)
Adaptive1=lala["Adaptive"][1]
s = open('Run_MTD_subtreat_0_2.txt', 'r').read()
lala=eval(s)
MTD0=lala["MTD"][0]
s = open('Run_MTD_subtreat_1_2.txt', 'r').read()
lala=eval(s)
MTD1=lala["MTD"][1]
s = open('Run_MTD_subtreat_2_2.txt', 'r').read()
lala=eval(s)
MTD2=lala["MTD"][2]
s = open('Run_NoTreatment_subtreat_0_2.txt', 'r').read()
lala=eval(s)
NoTreatment0=lala["NoTreatment"][0]
s = open('Run_Metronomic_subtreat_0_2.txt', 'r').read()
lala=eval(s)
Metronomic0=lala["Metronomic"][0]
s = open('Run_Metronomic_subtreat_1_2.txt', 'r').read()
lala=eval(s)
Metronomic1=lala["Metronomic"][1]
s = open('Run_Metronomic_subtreat_2_2.txt', 'r').read()
lala=eval(s)
Metronomic2=lala["Metronomic"][2]
s = open('Run_Metronomic_subtreat_3_2.txt', 'r').read()
lala=eval(s)
Metronomic3=lala["Metronomic"][3]

alldics=[NoTreatment0, MTD0, MTD1, MTD2, Metronomic0, Metronomic1, Metronomic2, Metronomic3, Adaptive0, Adaptive1]
labels=["NoTreatment", "MTD 0", "MTD 1", "MTD 2", "Metronomic 0", "Metronomic 1", "Metronomic 2", "Metronomic 3", "Adaptive 0", "Adaptive 1"]
treatmentColors=["firebrick","pink","yellowgreen","teal","lightblue","mediumpurple","gold","palevioletred", "slategray","orange"]



allvalues=[]
for j in range(len(alldics)):
    tmp=[]
    for i in range(22):
        tmp.append(alldics[j][i][2]["count"]/50)
    allvalues.append(tmp)
allvalues=np.array(allvalues)
neg=plt.imshow(allvalues)
plt.yticks(range(0,10),tuple(labels))
plt.xticks(range(0,22),range(1,23))
plt.colorbar(neg,label="Proportion of runs in which T- dominates the popultaion at final time")
plt.xlabel("Patient")
plt.ylabel("Treatment")
plt.show()




width=0.25
for i in range(22):
    a=[]
    b=[]
    c=[]
    for j in range(len(alldics)):
        a.append(np.sum(alldics[j][i][0]["densities"]))
        b.append(np.sum(alldics[j][i][1]["densities"]))
        c.append(np.sum(alldics[j][i][2]["densities"]))
        
    x = np.arange(len(labels))  # the label locations
    plt.bar(x - width, a, width, label='T+',color="b")
    plt.bar(x, b, width, label='TP',color="r")
    plt.bar(x + width, c, width, label='T-',color="g")
    plt.legend(["T+","TP","T-"])
    plt.xticks(range(10),labels,rotation=90)
    plt.title("Patient #%d"%(i+1))
    for h in range(len(np.arange(0.5,10.5,1))):
        plt.plot([np.arange(0.5,10.5,1)[h]]*2,[0,700],color="k")
    plt.show()
    

width=0.7
for i in range(22):
    tot=[]
    for j in range(len(alldics)):
        a=np.sum(alldics[j][i][0]["densities"])*alldics[j][i][0]["count"]/50
        b=np.sum(alldics[j][i][1]["densities"])*alldics[j][i][1]["count"]/50
        c=np.sum(alldics[j][i][2]["densities"])*alldics[j][i][2]["count"]/50
        tot.append(a+b+c)
        
    x = np.arange(len(labels))  # the label locations
    plt.bar(x, tot, width,color="gray")#,color=treatmentColors)#, label='TP',color="r")
    plt.xticks(range(10),labels,rotation=90)
    plt.title("Patient #%d"%(i+1))
    for h in range(len(np.arange(0.5,10.5,1))):
        plt.plot([np.arange(0.5,10.5,1)[h]]*2,[0,700],color="k")
    plt.show()
    
    
a=[3,4,5] #Each one of these is a treatment #Each number is a patient I think
b=[4,3,6]
colors=["b","r"]
tara=np.array([a,b])
sorted_ind=np.argsort(tara,axis=0)
for h in range(len(a)):
    orden=[tara[sorted_ind[0,h],h]]+list(np.diff(tara[sorted_ind[:,h],h]))
    color_orden=itemgetter(*sorted_ind[:,h])(colors)
    plt.bar([h],[orden[0]],color=color_orden[0])
    previous=[orden[0]]
    for j in range(1,len(orden)):
        plt.bar([h],[orden[j]],bottom=previous,color=color_orden[j])
        previous=orden[j]
plt.xticks(range(3),labels[0:3],rotation=90)
        




lolo=[]
yvalues=list(range(22))
for j in range(len(alldics)):
    xvalues=[]
    for i in range(22):
        xvalues.append(alldics[j][i][2]["Time"])
    lolo.append(xvalues)
    plt.scatter(xvalues,yvalues,color=treatmentColors[j])
    
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width*0.77, box.height])
legend_x = 1
legend_y = 0.5
plt.legend(labels, loc='center left', bbox_to_anchor=(legend_x, legend_y),title="Treatment",title_fontsize=12)
plt.xlabel("Time")
plt.ylabel("Patient")
plt.xlim([110,10000])
plt.yticks([])
plt.show()

lolo=np.array(lolo)
sorted_ind=np.argsort(lolo,axis=0)
for h in range(22):
    orden=[lolo[sorted_ind[0,h],h]]+list(np.diff(lolo[sorted_ind[:,h],h]))
    #orden=lolo[sorted_ind[:,h],h]
    color_orden=itemgetter(*sorted_ind[:,h])(treatmentColors)
    plt.barh([h],[orden[0]],color=color_orden[0])
    previous=[orden[0]]
    for j in range(1,len(orden)):
        plt.barh([h],orden[j],left=previous,color=color_orden[j])
        previous=np.cumsum(orden)[j]
        
axo = plt.subplot(111)
box = axo.get_position()
axo.set_position([box.x0, box.y0, box.width*0.77, box.height])
legend_x = 1
legend_y = 0.5
        
plt.yticks(range(22),range(1,23))
plt.xlabel("Time")
plt.ylabel("Patient Number")
plt.legend(labels, loc='center left', bbox_to_anchor=(legend_x, legend_y),title="Treatment",title_fontsize=12)
plt.show()


