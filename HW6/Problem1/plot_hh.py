# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 09:42:06 2018

@author: Alexander Bo Lee
"""

import matplotlib.pyplot as plt

def plot_hh(dyn,pars,savePDF=False,figName=''):
    '''
    function ploth = plot_hh(dyn, pars)
    
    Takes a HH output structure in dyn and plots dynamics of Voltage, gating variables, conductance
    currents and applied current
    '''

    fig, axes = plt.subplots(5,1,figsize = (9,12),gridspec_kw = {'height_ratios':[4,1,1,1,1]})
    
    
    
    # main data goes here
    
    # Applied current on bottomA
    ax4 = axes[4]
    ax4.plot(dyn['t'],dyn['appliedI'],color='k', linestyle='-',linewidth=2)
    ax4.set_xlabel('Time (ms)',fontsize=12)
    ax4.text(0,23,'Applied current',fontsize=20)
    ax4.set_ylim([0,30])
    plt.setp(ax4.spines.values(),linewidth=2)
    ax4.tick_params(labelsize = 20)
    
    # Next up is Currents
    ax3 = axes[3]
    ax3.plot(dyn['t'],dyn['IK'],color='k')
    ax3.plot(dyn['t'],dyn['INa'],color='r')
    ax3.plot(dyn['t'],dyn['IL'],color='b')
    ax3.text(0,300,'Currents',fontsize=20)
    ax3.set_ylim([-900,900])
    plt.setp(ax3.spines.values(),linewidth=2)
    ax3.tick_params(labelsize = 20)
    ax3.set_xticks([])
    
    # Conductance
    ax2 = axes[2]
    ax2.plot(dyn['t'],dyn['gK'],color='k')
    ax2.plot(dyn['t'],dyn['gNa'],color='r')
    ax2.plot(dyn['t'],dyn['gL'],color='b')
    ax2.text(0,23,'Conductances',fontsize=20)
    ax2.set_ylim([0,30])
    plt.setp(ax2.spines.values(),linewidth=2)
    ax2.tick_params(labelsize = 20)
    ax2.set_xticks([])
    
    #Gating
    ax1 = axes[1]
    ax1.plot(dyn['t'],dyn['n'],color='k')
    ax1.plot(dyn['t'],dyn['m'],color='r')
    ax1.plot(dyn['t'],dyn['h'],color='b')
    ax1.text(0,0.7,'Gating',fontsize=20)
    ax1.set_ylim([0,1])
    plt.setp(ax1.spines.values(),linewidth=2)
    ax1.tick_params(labelsize = 20)
    ax1.set_xticks([])
    
    # Voltage
    ax0 = axes[0]
    ax0.plot(dyn['t'],dyn['V'],color='k',linewidth=3)
    ax0.text(0,130,'Membrane voltage (mV)',fontsize=20)
    ax0.set_ylim([-20,150])
    ax0.text(0,pars['EK'],r'$E_K$',fontsize=20)
    ax0.text(0,pars['ENa'],r'$E_{Na}$',fontsize=20)
    ax0.text(0,pars['EL'],r'$E_L$',fontsize=20)
    plt.setp(ax0.spines.values(),linewidth=2)
    ax0.tick_params(labelsize = 20)
    ax0.set_xticks([])
    
    
    plt.tight_layout()
    plt.show()
    if savePDF:
        fig.savefig(figName+'.pdf',bbox_inches='tight')
    return [fig,axes]