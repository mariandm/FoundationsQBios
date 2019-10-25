# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:47:42 2018

@author: Alexander Bo Lee
"""
import numpy as np

# impulse_t 
def impulse_t(t,tryI,RefracTime):
    '''
    function I = impulse_t(t)
    
    specifies the appplied time-varying current
    works if t is a single value or many values
    '''
    secondTime1 = 2.5+RefracTime
    secondTime2 = 2.5+RefracTime+0.5
    
    if isinstance(t,float):
        if t>2 and t<2.5:
            I=tryI
        elif t>secondTime1 and t<secondTime2:
            I=tryI
        else:
            I=0
    else:
        I=np.zeros(len(t))
        for i in range(len(t)):
            if t[i]>2 and t[i]<2.5:
                I[i] = tryI
            elif t[i]>secondTime1 and t[i]<secondTime2:
                I[i] = tryI
            else:
                I[i]=0
    return I