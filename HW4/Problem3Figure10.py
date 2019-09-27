#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 22:52:09 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from Functions import *

from numpy import genfromtxt
my_data = genfromtxt('lang_data.csv', delimiter=',')


for i in my_data:
    x=(1000/11)*np.array(range(0,len(i)))
    plt.plot(x,i,marker=".", alpha=0.5)
plt.legend(["Mutant 1","Mutant 2","Mutant 3","Mutant 4","Mutant 5","Mutant 6","Mutant 7","Mutant 8","Mutant 9","Mutant 10","Mutant 11","Mutant 12","Mutant 13","Mutant 14","Mutant 15"],prop={'size': 7})
plt.xlim([0,1020])
plt.ylim([0,1.05])
plt.xlabel("Generations")
plt.ylabel("Frequency")
plt.title("Population G7 of Lang et al.")
plt.savefig("Figure7.pdf")