#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:13:55 2019

@author: mariandm
"""

import numpy as np
from scipy import integrate

def pulseChanges(t,tchanges,Achanges,vma):
    ind=max(np.where(t >=np.array(tchanges))[0])
    A=Achanges[ind]
    va = vma/(1+A)
    return(va)


def chemotaxis(x,t,vr,R,Kr,vd,vb,B,Kb,vma,Xtot,tchanges,Achanges):
    Xm=x[0]
    Xmp=x[1]
    X0= Xtot-Xm-Xmp

    va=pulseChanges(t,tchanges,Achanges,vma)
    
    dXmdt = vr*R*X0/(Kr+X0) + vd*Xmp - va*Xm
    dXmpdt = - vb*B*Xmp/(Kb+Xmp) - vd*Xmp + va*Xm
    
    dxdt=np.array([dXmdt,dXmpdt])
    
    return(dxdt)