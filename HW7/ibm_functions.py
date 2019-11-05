# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 19:17:23 2019

@author: Alex Bo Lee
Key to visualization is found on Stack Overflow
https://stackoverflow.com/questions/11874767/how-do-i-plot-in-real-time-in-a-while-loop-using-matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt

def ibm_predation(info,predator,prey):
    """
    [t,numeaten] = ibm_predation(info,predator,prey)
    An individual-based movement model of predators eating prey
    numeaten returns the time series of eating
    """
    
    #Python inputs are not local copies, but pointers to the real thing. If we
    # want to run multiple times, we'll have to reset certain parts
    prey['num'] = info['prey_density']*info['maxX']*info['maxY']
    prey['pos']=np.random.uniform(size=(prey['num'],2))
    prey['pos'][:,0]=prey['pos'][:,0]*info['maxX']
    prey['pos'][:,1]=prey['pos'][:,1]*info['maxY']
    predator['pos'] = np.array( [[ info['maxX']/2, info['maxY']/2 ]] )
    predator['theta']=np.random.uniform()*2*np.pi #Angle of movement
    
    t=np.arange(0,info['tf']+info['dt'],step=info['dt'])
    numt = len(t)
    numeaten = np.zeros(shape=t.shape)
    
    ispeating = 0 # whether predators are eating
    pdoneeating = 0 # The time when predators are done eating
    numprey = len(prey['pos'][:,0])
    
    # Setup the plot
    if info['viz_dyn']:
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim([0,info['maxX']])
        ax.set_ylim([0,info['maxY']])
        
        
        # Draw the predator as it views things
        tmppred_scan = plt.Circle(predator['pos'][0,:],predator['r'],color='g')
        ax.add_patch(tmppred_scan)
        
        # Draw the predator
        tmppred_body = plt.plot(predator['pos'][0,0],predator['pos'][0,1],'ro')
        plt.setp(tmppred_body,markerfacecolor='r')
        
        # Draw the prey
        tmpprey = plt.plot(prey['pos'][:,0],prey['pos'][:,1],'k.')
        
        carcassFlag = False #Should we show a carcass?
        # Keep the same axes
        ax.set_aspect('equal')
        plt.pause(0.001)
        
        # Update the time
        tmps = 't = {time}'.format(time=0)
        tmptime_handle = plt.text(info['maxX']*1.05,info['maxY'],tmps)
        plt.setp(tmptime_handle,fontsize=24)
    
    for i in np.arange(numt):
        if (ispeating and t[i]>pdoneeating): # Is the predator done eating?
            ispeating = 0
        if not ispeating:
            #move predator
            predator['pos'] = (predator['pos']+
                    predator['vel']*
                    np.array([[np.cos(predator['theta']),
                              np.sin(predator['theta'])]])*
                    info['dt'])
            predator['pos'] = restore_pos(predator['pos'],
                    info['maxX'],info['maxY'])
            predator['trun'] = predator['trun']+info['dt']
            if predator['trun']>predator['tau']:
                predator['theta']=np.random.uniform()*2*np.pi
                predator['trun']=0
        
        # Move prey
        if numprey>0:
            rdxy = np.random.normal(size=(prey['num'],2))*\
                np.sqrt(prey['diffusion'])*np.sqrt(info['dt'])
            prey['pos'] = prey['pos']+rdxy
            prey['pos'] = restore_pos(prey['pos'],info['maxX'],info['maxY'])

            
            # Find who is eaten
            if not ispeating:
                #Find prey in range, this is tricky due to the periodic
                # boundary condition
                distX = np.abs(predator['pos'][:,0]-prey['pos'][:,0])
                distY = np.abs(predator['pos'][:,1]-prey['pos'][:,1])
                distX_optimal = np.amin(np.array([distX,info['maxX']-distX]),0)
                distY_optimal = np.amin(np.array([distY,info['maxY']-distY]),0)
                tmpi = np.where(distX_optimal**2+distY_optimal**2<
                                (predator['r']+predator['vel']*info['dt'])**2)
#                tmpi = np.where((np.square(predator['pos'][:,0]-
#                                           prey['pos'][:,0])+
#                                np.square(predator['pos'][:,1]-
#                                          prey['pos'][:,1])) <
#                                (predator['r']+predator['vel']*info['dt'])**2)
                tmpi = tmpi[0] #lose the metadata
                if len(tmpi)>0: # Someone is detected
                    tmpj = (np.random.uniform(size=len(tmpi))<
                            predator['k']*info['dt']*predator['f'])
                    if predator['handling_time'] == 0:
                        numeaten[i]=sum(tmpj) # Eaten
                    elif sum(tmpj) > 0:
                        numeaten[i]=1;
                        #tmpj=[1]; Why is this here?
                        pdoneeating=t[i]+predator['handling_time']
                        ispeating=1
                    else:
                        numeaten[i] = 0
        else:
            numeaten[i] = 0
    
        if numeaten[i] > 0:
            # Remove the old visualization
            if 'tmph_eaten' in locals():
                del tmph_eaten
            if info['viz_dyn']:
                carcassFlag = True
                carcassX = prey['pos'][tmpi[tmpj][0],0]
                carcassY = prey['pos'][tmpi[tmpj][0],1]
                
            if info['replenish_prey']:
                prey['pos'][tmpi[tmpj],0]=(np.random.uniform(size=sum(tmpj))*
                    info['maxX'])
                prey['pos'][tmpi[tmpj],1]=(np.random.uniform(size=(sum(tmpj),))*
                    info['maxY'])
            else:
                prey['pos']=np.delete(prey['pos'],tmpi[tmpj],axis=0)
                prey['num']=prey['num']-sum(tmpj)
        
        # Plot the predator and prey
        
        # Re-visualize
        if info['viz_dyn']:
            # clean up the old visualization and reset
            plt.pause(0.001)
            ax.clear()
            ax.set_xlim([0,info['maxX']])
            ax.set_ylim([0,info['maxY']])
            ax.set_aspect('equal')
            
            tmppred_scan = plt.Circle(predator['pos'][0,:],predator['r'],color='g')
            ax.add_patch(tmppred_scan)
            plt.plot(predator['pos'][:,0],predator['pos'][:,1],'ro')
            plt.setp(tmppred_body, markerfacecolor='r')
            plt.plot(prey['pos'][:,0],prey['pos'][:,1],'k.')
            tmps = 't = {time}'.format(time=t[i])
            plt.text(info['maxX']*1.05,info['maxY'],tmps,fontsize=24)
            if carcassFlag:
                plt.plot(carcassX,carcassY,
                         'bo',markersize=8,markerfacecolor='b')
    plt.show()
        
    return [t,numeaten]            


def restore_pos(curpos,maxX,maxY):
    '''
    function newpos = resotre_pos(curpos,maxX,maxY)
    
    Enusres that the position in curpos (X,Y) is between maxX, maxY,
    by subtracting/adding the maxX/maxY as appropriate
    '''
    newpos = np.array(curpos)
    tmpi = np.where(newpos[:,0]<0)[0]
    if len(tmpi) is not 0:
        newpos[tmpi,0] = newpos[tmpi,0]+maxX
    tmpi = np.where(newpos[:,0]>maxX)[0]
    if len(tmpi) is not 0:
        newpos[tmpi,0] = newpos[tmpi,0]-maxX
    tmpi = np.where(newpos[:,1]<0)[0]
    if len(tmpi) is not 0:
        newpos[tmpi,1] = newpos[tmpi,1]+maxY
    tmpi = np.where(newpos[:,1]>maxY)[0]
    if len(tmpi) is not 0:
        newpos[tmpi,1] = newpos[tmpi,1]-maxY
    return newpos

            
            