#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 22:16:43 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

a = range(5,180,10)
b = np.array([18,55,71,120,140,138,135,93,89,58,52,53,39,41,18,11,10,3])
#plt.plot(a,b)

chales=[]
for h in range(0,len(a)):
    chales=chales+([a[h]]*b[h])

#e = np.mean(chales)
#v = np.var(chales)
#alpha = (-(e**2)*(e-1)-v*e)/v
#beta = -(e*alpha-alpha)/e


h1,h2,h3,h4=stats.beta.fit(chales)
fitb=stats.beta.pdf(a,h1,h2,h3,h4)
por=max(b)/max(fitb)

#plt.hist(np.random.beta(h1,h2,1000))
plt.scatter(a,b)
plt.plot(x,fitb*por)
plt.legend(["Fitted Beta","Data"])
plt.xlabel("Degrees")
plt.ylabel("Frequency")
print(h1)
print(h2)
#plt.savefig("Figure5.pdf")
plt.close()



######Proof change in angles
##ORIGINAL VERSION
ori_diffangle=[]

tmax = 100 # seconds
v = 1 # microns/sec
currangle = np.random.uniform()*2*np.pi
tumblerate = 2 # per sec
# Initialize the system
currt = 0
currx = 0
curry = 0
Xoft = [[currx, curry]]
t = [0]
ind = 0
# Run and tumble

# Parameters and constraints
tmax = 10000 # seconds
v = 1 # microns/sec
currangle = np.random.uniform()*2*np.pi
tumblerate = 2 # per sec
# Initialize the system
currt = 0
currx = 0
curry = 0
Xoft = [[currx, curry]]
t = [0]
ind = 0
while currt<tmax:
    ind = ind+1
    dt = -(1/tumblerate)*np.log(np.random.uniform()) # Find the next tumble time
    dx = v*np.cos(currangle)*dt #Advance in the x-direction
    dy = v*np.sin(currangle)*dt #Advance in the y-direction
    currt = currt + dt
    currx = currx + dx
    curry = curry + dy
    t.append(currt)
    Xoft.append([currx,curry])
    newangle = np.random.uniform()*2*np.pi
    ori_diffangle.append(abs(currangle-newangle))
    currangle = newangle
    # convert Xoft into a numpy array
Xoft = np.array(Xoft)
t = np.array(t)
    
plt.hist(ori_diffangle,alpha=0.6, label='Original',normed=True)


##MODIFIED VERSION
mod_diffangle=[]
# Parameters and constraints
tmax = 10000 # seconds
v = 1 # microns/sec
currangle = np.random.uniform()*2*np.pi
tumblerate = 2 # per sec
# Initialize the system
currt = 0
currx = 0
curry = 0
Xoft = [[currx, curry]]
t = [0]
ind = 0
while currt<tmax:
    ind = ind+1
    dt = -(1/tumblerate)*np.log(np.random.uniform()) # Find the next tumble time
    dx = v*np.cos(currangle)*dt #Advance in the x-direction
    dy = v*np.sin(currangle)*dt #Advance in the y-direction
    currt = currt + dt
    currx = currx + dx
    curry = curry + dy
    t.append(currt)
    Xoft.append([currx,curry])
    #calculating \beta from a beta distr
    diffangle = (np.random.beta(h1,h2)*h4+h3)*np.random.choice([-1,1])
    diffangle=np.radians(diffangle)
    mod_diffangle.append(abs(diffangle))
    currangle = currangle + diffangle
    # convert Xoft into a numpy array
Xoft = np.array(Xoft)
t = np.array(t)


plt.hist(mod_diffangle,alpha=0.6, label='Modified',color="darkgreen",normed=True)
plt.legend()
plt.xlabel("Angle change in radians")
plt.ylabel("Frequency")
#plt.savefig("Figure6.pdf")
plt.close()


##### 3 examples of each 
fig, axs = plt.subplots(1,2,figsize=(8, 3))

tmax = 100 # seconds
v = 25 # microns/sec
currangle = np.random.uniform()*2*np.pi
tumblerate = 1 # per sec
# Initialize the system
currt = 0
currx = 0
curry = 0
Xoft = [[currx, curry]]
t = [0]
ind = 0
# Run and tumble
for i in range(0,3):
    # Parameters and constraints
    tmax = 100 # seconds
    v = 1 # microns/sec
    currangle = np.random.uniform()*2*np.pi
    tumblerate = 2 # per sec
    # Initialize the system
    currt = 0
    currx = 0
    curry = 0
    Xoft = [[currx, curry]]
    t = [0]
    ind = 0
    while currt<tmax:
        ind = ind+1
        dt = -(1/tumblerate)*np.log(np.random.uniform()) # Find the next tumble time
        dx = v*np.cos(currangle)*dt #Advance in the x-direction
        dy = v*np.sin(currangle)*dt #Advance in the y-direction
        currt = currt + dt
        currx = currx + dx
        curry = curry + dy
        t.append(currt)
        Xoft.append([currx,curry])
        currangle = np.random.uniform()*2*np.pi
    # convert Xoft into a numpy array
    Xoft = np.array(Xoft)
    t = np.array(t)
    
    axs[0].plot(Xoft[:,0],Xoft[:,1])
axs[0].legend(["Run 1","Run 2","Run 3"])
axs[0].set_title("Original Model")

print(len(t))

tmax = 100 # seconds
v = 25 # microns/sec
currangle = np.random.uniform()*2*np.pi
tumblerate = 1 # per sec
# Initialize the system
currt = 0
currx = 0
curry = 0
Xoft = [[currx, curry]]
t = [0]
ind = 0
# Run and tumble
for i in range(0,3):
    # Parameters and constraints
    tmax = 100 # seconds
    v = 1 # microns/sec
    currangle = np.random.uniform()*2*np.pi
    tumblerate = 2 # per sec
    # Initialize the system
    currt = 0
    currx = 0
    curry = 0
    Xoft = [[currx, curry]]
    t = [0]
    ind = 0
    while currt<tmax:
        ind = ind+1
        dt = -(1/tumblerate)*np.log(np.random.uniform()) # Find the next tumble time
        dx = v*np.cos(currangle)*dt #Advance in the x-direction
        dy = v*np.sin(currangle)*dt #Advance in the y-direction
        currt = currt + dt
        currx = currx + dx
        curry = curry + dy
        t.append(currt)
        Xoft.append([currx,curry])
        diffangle = (np.random.beta(h1,h2)*h4+h3)*np.random.choice([-1,1])
        diffangle=np.radians(diffangle)
        currangle = currangle + diffangle
    # convert Xoft into a numpy array
    Xoft = np.array(Xoft)
    t = np.array(t)
    
    axs[1].plot(Xoft[:,0],Xoft[:,1])
axs[1].legend(["Run 1","Run 2","Run 3"])
axs[1].set_title("Modified Model")
#fig.savefig("Figure7.pdf")
plt.close()
print(len(t))


########
#### MSD
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




###Ensemble of samples for MSD and shit
nsamples=200

MSD=np.zeros([1000,2])
#We are gonna have to sample as in no me acuerdo que tarea para poder tener 
for i in range(0,nsamples):
    tmax = 1000 # seconds
    v = 1 # microns/sec
    currangle = np.random.uniform()*2*np.pi
    tumblerate = 0.4 # per sec
    # Initialize the system
    currt = 0
    currx = 0
    curry = 0
    Xoft = [[currx, curry]]
    t = [0]
    ind = 0
    while currt<tmax:
        ind = ind+1
        dt = -(1/tumblerate)*np.log(np.random.uniform()) # Find the next tumble time
        dx = v*np.cos(currangle)*dt #Advance in the x-direction
        dy = v*np.sin(currangle)*dt #Advance in the y-direction
        currt = currt + dt
        currx = currx + dx
        curry = curry + dy
        t.append(currt)
        Xoft.append([currx,curry])
        currangle = np.random.uniform()*2*np.pi
        # convert Xoft into a numpy array
    Xoft = np.array(Xoft)
    t = np.array(t)
    #Aqui sampleamos y eso lo sumamos a MSD
    #We can do just a square cause we ALWAYS start at 0 c:
    a=sample_traj(t,Xoft[:,0],range(0,1000))
    MSD[:,0]+=a[1]**2
    a=sample_traj(t,Xoft[:,1],range(0,1000))
    MSD[:,1]+=a[1]**2
#Aqui dividmos MSD/nsamples
MSD=MSD/nsamples
MSD=MSD[:,0]+MSD[:,1]
MSD=MSD**(1/2)
time=np.array(range(0,1000))
root_time=time**(1/2)

fig, axs = plt.subplots(1,2,figsize=(8, 4))
axs[0].plot(time,MSD)
axs[0].set_ylabel("Root MSD")
axs[0].set_xlabel("Time")
axs[1].plot(root_time,MSD)
axs[1].set_ylabel("Root MSD")
axs[1].set_xlabel("Root Time")

fig.savefig("Figure8.pdf")
plt.close()