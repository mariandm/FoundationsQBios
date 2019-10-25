#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 10:00:25 2019

@author: mariandm
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML


pars={}
# Length of cable
Lx = 100
# Number of points for approximation
Nx = 100
dx = Lx/Nx
# Total cells including ghost cells
totN = Nx+2
#Total time of dynamics
T=200
# Parameters from ODE
pars['a'] = 0.1
pars['b'] = 0.01
pars['gamma'] = 0.02
pars['I'] = 0#.125
pars['D'] = 1 # Diffusion constant
dt = 0.025
totruns = int(np.round(T/dt))
# Preallocate
V = np.zeros(totN)
w = np.zeros(totN)
data = {}
data['allV'] = np.zeros((Nx,totruns+1))
data['allW'] = np.zeros((Nx,totruns+1))
data['tvec'] = np.arange(0,T+dt,dt)



V[11]=2.25

t = 0
cnt = 0
while t<=(T+dt):
    # calculate terms in PDE
    #reaction term in V in FN
    rxnV = V*(pars['a']-V)*(V-1)-w+pars['I']
    # Voltage diffusion
    diffV = pars['D']*(V[:-2]-2*V[1:-1]+V[2:])/dx**2
    #reaction term of w in FN
    rxnw = pars['b']*V-pars['gamma']*w
    # Evolve over time
    V[1:-1] = dt*(rxnV[1:-1]+diffV) + V[1:-1]######
    w[1:-1] = dt*(rxnw[1:-1]) + w[1:-1]
    data['allV'][:,cnt] = V[1:-1]
    data['allW'][:,cnt] = w[1:-1]
    # Fix boundary conditions for no-flux
    V[0] = V[2]
    V[-1] = V[-3]
    # not necessary for w, but we'll do it anyways for consistency
    w[0] = w[2]
    w[-1] = w[-3]
    # update time
    t = t+dt
    cnt = cnt+1
    

    #x = np.arange(100)
    #y = data['allV'][:,cnt-1]
    #plt.plot(x,y)
    #plt.ylim([0,2.25])
    #plt.show()    

# First set up the figure, the axis, and the plot elements we
# want to animate
fig, ax = plt.subplots()
ax.set_xlim([0,100])
ax.set_ylim(np.amin(data['allV']),np.amax(data['allV']))
ax.set_ylabel('Voltage, V',fontsize=20)
ax.set_xlabel('Cell in Cable',fontsize=20)
plt.tight_layout()
line, = ax.plot([],[],linewidth=2)
# initialization function: plot the background of each frame
def init():
    line.set_data([],[])
    return (line,)
# animation function. This is called sequentially
def animate(i):
    x = np.arange(100)
    y = data['allV'][:,i]
    line.set_data(x,y)
    return (line,)
# call the animator. blit=True means only re-draw the parts that
# have changed.
#anim = animation.FuncAnimation(fig,animate, init_func=init,frames=4000,interval=20,blit=True)
#anim.save('im.mp4', writer='ffmpeg',bitrate=-1)
plt.close()

######################################
#######4a SPEED

stime=5
sruns=int(stime/dt)
lala=[]
for i in range(sruns,totruns,sruns):
    future=np.argmax(data['allV'][range(10,100),i])
    present=np.argmax(data['allV'][range(10,100),i-sruns])
    lala.append((future-present)/stime)
    if (future-present)<0:
        print(future,present)
        print(i)

plt.scatter(range(stime,T,stime),lala,color="lightcoral")
plt.plot(range(stime,T,stime),lala,color="lightcoral")
plt.ylim([-0.05,1])
plt.xlabel("Time",fontsize=15)
plt.ylabel("Speed,cell/s",fontsize=15)
plt.savefig("Figure6.pdf")
plt.show()

#####################################
#########4B D affecting traveling

#Vary D and plot different timepoints
Ds=[.5,.75,1.25,1.5]
fig, axs = plt.subplots(4,4,figsize=(10,6))
for de in range(0,len(Ds)):
    pars={}
    # Length of cable
    Lx = 100
    # Number of points for approximation
    Nx = 100
    dx = Lx/Nx
    # Total cells including ghost cells
    totN = Nx+2
    #Total time of dynamics
    T=100
    # Parameters from ODE
    pars['a'] = 0.1
    pars['b'] = 0.01
    pars['gamma'] = 0.02
    pars['I'] = 0#.125
    pars['D'] = Ds[de] # Diffusion constant
    dt = 0.025
    totruns = int(np.round(T/dt))
    # Preallocate
    V = np.zeros(totN)
    w = np.zeros(totN)
    data = {}
    data['allV'] = np.zeros((Nx,totruns+1))
    data['allW'] = np.zeros((Nx,totruns+1))
    data['tvec'] = np.arange(0,T+dt,dt)
    
    
    
    V[11]=2.25
    
    t = 0
    cnt = 0
    plotid=0
    while t<=(T+dt):
        # calculate terms in PDE
        #reaction term in V in FN
        rxnV = V*(pars['a']-V)*(V-1)-w+pars['I']
        # Voltage diffusion
        diffV = pars['D']*(V[:-2]-2*V[1:-1]+V[2:])/dx**2
        #reaction term of w in FN
        rxnw = pars['b']*V-pars['gamma']*w
        # Evolve over time
        V[1:-1] = dt*(rxnV[1:-1]+diffV) + V[1:-1]######
        w[1:-1] = dt*(rxnw[1:-1]) + w[1:-1]
        data['allV'][:,cnt] = V[1:-1]
        data['allW'][:,cnt] = w[1:-1]
        # Fix boundary conditions for no-flux
        V[0] = V[2]
        V[-1] = V[-3]
        # not necessary for w, but we'll do it anyways for consistency
        w[0] = w[2]
        w[-1] = w[-3]
        # update time
        t = t+dt
        cnt = cnt+1
        
        if(cnt==int(1/dt) or cnt==int(25/dt) or cnt==int(50/dt) or cnt==int(90/dt)):
            x = np.arange(100)
            y = data['allV'][:,cnt-1]
            axs[plotid][de].plot(x,y)
            axs[plotid][de].set_ylim([-0.1,2.25])
            axs[plotid][de].set_ylabel("V, Voltage")
            axs[plotid][de].set_xlabel("Cell in cable")
            
            #plt.show() 
            plotid=plotid+1
            

times=[1,25,50,90]
for de in range(0,len(Ds)):          
    lala=Ds[de]
    axs[0][de].set_title("D=%.2f" %lala,fontsize=20)
    secax = axs[de][3].secondary_yaxis('right')  
    secax.set_ylabel("Time=%d" %times[de])       
for ax in fig.get_axes():
    ax.label_outer()
fig.tight_layout()            
fig.savefig("Figure7.pdf")
plt.show()
    

###YA NADA MAS HACER UNA MEGA FIGURA DE 4*4 CON TIME EN EL EJE DE LAS Y Y D EN EL DE LAS X