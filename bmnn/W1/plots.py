# -*- coding: utf-8 -*-
"""a function of ploting figures."""
import numpy as np
import matplotlib.pyplot as plt


def correlation_plot(Cx_norm_W,w,t_max,labely):
    t=[i for i in range(0,t_max+1)]
    
    plt.figure(figsize=(8, 6), dpi=80)
    
    plt.plot(t,Cx_norm_W[0,:],linestyle="-", color='b',label="gamma/m=0.5w")
    plt.plot(t,Cx_norm_W[1,:],linestyle="-", color='g',label="gamma/m=1.0w")
    plt.plot(t,Cx_norm_W[2,:],linestyle="-", color='r',label="gamma/m=1.5w")
     
    plt.xlabel("t")
    plt.ylabel(labely)
    plt.legend(loc='upper right')
    
    
    
    plt.savefig(labely, dpi=72)
        

          

def mean_plot(v,t,vlabel):
    
    moy=np.mean(v,axis=0)
    
    plt.plot(t,moy,linestyle="-", color='b')
             
    plt.xlabel("t")
    plt.ylabel(vlabel)

def Qmean_plot(v,t,C,labelC,labelv):
    Cline=C*np.ones((len(t),1))
    vmoy2=(np.mean(v,axis=0))**2
    v_var=np.var(v,axis=0)
    v2=v_var+vmoy2
    
    
    plt.plot(t,Cline,linestyle="-", color='r',label=labelC)
    plt.plot(t,v2,linestyle="-", color='b',label=labelv)
             
    plt.xlabel("t")
    plt.legend(loc='center right')


def position(x,t,ci):
    """visualization position"""
    
    # Create a figure of size 8x6 inches, 80 dots per inch
    #plt.figure(figsize=(8, 6), dpi=80)

    # Create a new subplot from a grid of 1x1
    #plt.subplot(1, 1, 1)
    
    
    
    plt.plot(t, x,linestyle="-", color=((1-ci),0.6,ci)) # linewidth=1.0,label="position"), marker="+"

    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.title("position in function of time")
    
    
    #plt.legend(loc='upper left') affiche les label des lignes
    
    # Set x limits
    #plt.xlim(-4.0, 4.0)
    
    # Set x ticks
    #plt.xticks(np.linspace(-4, 4, 9, endpoint=True))    9tickde -4 a 4
    #plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$+\pi/2$', r'$+\pi$']) affiche avec pi


    # Set y limits
    #plt.ylim(-1.0, 1.0)

    # Set y ticks
    #plt.yticks(np.linspace(-1, 1, 5, endpoint=True))


    
    #plt.grid(True)
    
    # Show result on screen
    #plt.show()
    
    # Save figure using 72 dots per inch
    #plt.savefig("position_h", dpi=72)


