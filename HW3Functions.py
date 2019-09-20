#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:51:40 2019

@author: mariandm
"""

import numpy as np
from scipy import integrate

def selfregulation(x,t,pars):
    dxdt=np.heaviside(x-pars['K'],1)*(pars['beta_plus']-pars['beta_minus'])+pars['beta_minus']-x*pars['alpha']
    return(dxdt)
    
def stochSelfRegulation(beta_plus,beta_minus,alpha,K,p0,tmax):
    #Initial conditions
    currp = p0 #current proteins
    currt = 0 #current time
    tvals = [currt] #recording of time
    pvals = [currp] #recording of proteins
    #Stochastic process
    while currt<tmax:
        # Calculate rates
        decayrate = alpha*currp
        prodrate = (beta_plus-beta_minus)*((currp-K>0)*1)+beta_minus
        #if(prodrate==beta_plus):
        #break
        total_rate = decayrate+prodrate
        # Find waiting time
        deltat = np.random.exponential(1/total_rate) # Find next event
        currt = currt+deltat # Move forward in time
        prob_produce = prodrate/total_rate # Calculate probability of production
        # Update the state
        if (np.random.uniform()<prob_produce): # Choose production at random
            currp = currp+1 # Increment number of proteins by 1
        else:
            currp = currp-1 # Decrement number of proteins by 1
        # Record the state
        tvals.append(currt) #Store the event time
        pvals.append(currp) #Store the protein number
    
    return(tvals,pvals)

def stochStateTime(beta_plus,beta_minus,alpha,K,p0,tmax):
    #Initial conditions
    currp = p0 #current proteins
    currt = 0 #current time
    #Stochastic process
    while currt<tmax:
        # Calculate rates
        decayrate = alpha*currp
        prodrate = (beta_plus-beta_minus)*((currp-K>0)*1)+beta_minus
        if(p0==beta_minus and prodrate==beta_plus): #starts in off
            return(currt)
        if(p0==beta_plus and prodrate==beta_minus): #starts in on
            return(currt)
        total_rate = decayrate+prodrate
        # Find waiting time
        deltat = np.random.exponential(1/total_rate) # Find next event
        currt = currt+deltat # Move forward in time
        prob_produce = prodrate/total_rate # Calculate probability of production
        # Update the state
        if (np.random.uniform()<prob_produce): # Choose production at random
            currp = currp+1 # Increment number of proteins by 1
        else:
            currp = currp-1 # Decrement number of proteins by 1
    
    return()


def sample_traj(t,y,trange):
    """"Samples a trajectory t,y at discrete intervals"""
    t=np.array(t)
    y=np.array(y)
    ts = np.zeros(np.shape(trange))
    ys = np.zeros(np.shape(trange))
    #initialize
    curt = trange[0]
    ind=np.where(t<=curt)[0]
    curind = ind[-1]
    ts[0]=t[curind]
    ys[0]=y[curind]
    #Scans across the sample interval and then moves the
    # actual dynamics forward until we cross it
    for i in range(1,len(trange)):
        nextt=trange[i]
        while t[curind]<nextt:
            curind=curind+1
        ts[i]=nextt
        ys[i]=y[curind-1]
    return ts, ys

def stochMRNAdetPRO(alpha,beta,m0,tmax,r):
    currt=0
    currm = m0 #current proteins
    currp = m0*r
    tvals = [currt] #recording of time
    mvals = [currm] #recording of proteins
    pvals = [currp] #recording of proteins
    while currt<tmax:
        # Calculate rates
        decayrate = alpha*currm
        prodrate = beta
        total_rate = decayrate+prodrate
        # Find waiting time
        deltat = np.random.exponential(1/total_rate) # Find next event
        currt = currt+deltat # Move forward in time
        prob_produce = prodrate/total_rate # Calculate probability of production
        # Update the state
        if (np.random.uniform()<prob_produce): # Choose production at random
            currm = currm+1 # Increment number of proteins by 1
        else:
            currm = currm-1 # Decrement number of proteins by 1
        # Record the state
        currp=currm*r
        tvals.append(currt) #Store the event time
        mvals.append(currm) #Store the mRNA number
        pvals.append(currp) #Store the protein number
    return(tvals,mvals,pvals)
    
def stochMRNAstochPRO(alpha,beta,m0,p0,tmax,r):
    currm = m0 #current proteins
    currp = p0
    currt = 0 #current time
    tvals = [currt] #recording of time
    mvals = [currm] #recording of proteins
    pvals = [currp] #recording of proteins
    #end of part 1 - initial conditions
    while currt<tmax:
        # Calculate rates
        mdecayrate = alpha*currm
        mprodrate = beta
        pdecayrate = alpha*currp
        pprodrate = r*currm

        total_rate = mdecayrate+mprodrate+pdecayrate+pprodrate
        # Find waiting time
        deltat = np.random.exponential(1/total_rate) # Find next event
        currt = currt+deltat # Move forward in time
        prob_m = (mdecayrate+mprodrate)/total_rate # Calculate probability of production
        # Update the state
        ran=np.random.uniform()
        if (ran<prob_m): # Choose production at random
            prob_prodm = mprodrate/total_rate
            if (ran<prob_prodm):
                currm = currm+1 # Increment number of proteins by 1
            else:
                currm = currm-1 # Decrement number of proteins by 1
        else:
            ran=ran-prob_m
            prob_prodp = pprodrate/total_rate
            if (ran<prob_prodp):
                currp = currp+1 # Increment number of proteins by 1
            else:
                currp = currp-1 # Decrement number of proteins by 1
        # Record the state
        tvals.append(currt) #Store the event time
        mvals.append(currm) #Store the mRNA number
        pvals.append(currp) #Store the protein number
    
    return(tvals,mvals,pvals)
    
def stochPaper(alpha,beta,offon,onoff,m0,tmax,r,env):
    currm = m0 #current proteins
    currp = m0*r
    currt = 0 #current time
    tvals = [currt] #recording of time
    mvals = [currm] #recording of proteins
    pvals = [currp] #recording of proteins
    ss=[env]
    
    states=[offon,onoff]
    betas=[0,beta]
    while currt<tmax:
        # Calculate rates
        mdecayrate = alpha*currm
        mprodrate = betas[env]
        staterate=states[env]
        total_rate = mdecayrate+mprodrate+staterate
        # Find waiting time
        deltat = np.random.exponential(1/total_rate) # Find next event
        currt = currt+deltat # Move forward in time
        #calculating probabilities
        mprodprob=mprodrate/total_rate
        mdecprob=mdecayrate/total_rate
        staprob=staterate/total_rate
        # Update the state
        ran=np.random.uniform()
        if (ran<mprodprob): # Choose production at random
            #production of mrna
            currm=currm+1
            currp=currm*r
        elif (mprodprob<ran<mprodprob+mdecprob):
            #decay
            currm=currm-1
            currp=currm*r
        else:
            #change state
            env=abs(env-1)
        # Record the state
        tvals.append(currt) #Store the event time
        mvals.append(currm) #Store the mRNA number
        pvals.append(currp) #Store the protein number
        ss.append(env) #Store on off
    
    return(tvals,mvals,pvals,ss)