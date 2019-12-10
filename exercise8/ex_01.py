"""
Created on Mon Dec  9 15:13:09 2019
@author: dominik
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.grid'] = True


def num_integration(L, T, NT, initial, technique, plot=False):
    dx = L/(len(initial)-1.)
    dt = T/(NT-1.)
    evolution = np.copy(initial)
    
    if plot:
        result = np.zeros((int(np.floor(NT/10.)+2.), len(initial)))
        result[0,:] = np.copy(initial)
        if (technique == 1):
            for i in range(1, NT):
                for k in range(1, (len(initial)-1)):
                    evolution[k] -= (evolution[k+1]-evolution[k-1])/(2*dx)*dt
                if (i%10 == 0):
                    result[int(i/10),:] = np.copy(evolution)
            result[-1,:] = np.copy(evolution)
            return result
        
        elif (technique == 2):
            for i in range(1, NT):
                for k in range(1, len(initial)):
                    evolution[k] -= (evolution[k]-evolution[k-1])/dx*dt
                evolution[-1] = initial[-1]
                if (i%10 == 0):
                    result[int(i/10),:] = np.copy(evolution)
            result[-1,:] = np.copy(evolution)
            return result
        
        elif (technique == 3):
            for i in range(1, NT):
                for k in range(0, (len(initial)-1)):
                    evolution[k] -= (evolution[k+1]-evolution[k])/dx*dt
                evolution[0] = initial[0]
                if (i%10 == 0):
                    result[int(i/10),:] = np.copy(evolution)
            result[-1,:] = np.copy(evolution)
            return result   
        
    else:
        if (technique == 1):
            for i in range(1, NT):
                for k in range(1, (len(initial)-1)):
                    evolution[k] -= (evolution[k+1]-evolution[k-1])/(2*dx)*dt
            return evolution
        
        elif (technique == 2):
            for i in range(1, NT):
                for k in range(1, len(initial)):
                    evolution[k] -= (evolution[k]-evolution[k-1])/dx*dt
                evolution[-1] = initial[-1]
            return evolution
        
        elif (technique == 3):
            for i in range(1, NT):
                for k in range(0, (len(initial)-1)):
                    evolution[k] -= (evolution[k+1]-evolution[k])/dx*dt
                evolution[0] = initial[0]
            return evolution
    
    
if __name__ == "__main__":
    
    L = 10.
    NL = 100
    T = 3.
    NT = 101
    
    bound_low = 1.
    bound_up = 0.
    
    initial = np.zeros(NL+2)
    initial[:int((NL+2)/2.)] = 1.
    initial[0] = bound_low
    initial[-1] = bound_up
    
    evolution_method_1 = num_integration(L, T, NT, initial, 1, True)
    evolution_method_2 = num_integration(L, T, NT, initial, 2, True)
    evolution_method_3 = num_integration(L, T, NT, initial, 3, True)
    
    x = np.linspace(-L/2., L/2., num=(NL+2))
    
    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(131)
    ax1.set_title('symmetrical difference quotient')
    ax1.set_ylabel('q(x,t) in arbitrary units')
    ax1.set_xlabel('x in arbitrary units')
    for i in range(len(evolution_method_1[:,0])):
        ax1.plot(x, evolution_method_1[i,:], label='t = {}'.format(int(10*i)))
  
    ## Stability due to additional diffusion term which damps oscillations.
    ax2 = fig.add_subplot(132)
    ax2.set_title('up-stream difference quotient')
    ax2.set_ylabel('q(x,t) in arbitrary units')
    ax2.set_xlabel('x in arbitrary units')
    for i in range(len(evolution_method_2[:,0])):
        ax2.plot(x, evolution_method_2[i,:], label='t = {}'.format(int(10*i)))
        
    ax3 = fig.add_subplot(133)
    ax3.set_title('down-stream difference quotient')
    ax3.set_ylabel('q(x,t) in arbitrary units')
    ax3.set_xlabel('x in arbitrary units')
    for i in range(len(evolution_method_3[:,0])):
        ax3.plot(x, evolution_method_3[i,:], label='t = {}'.format(int(10*i)))
    
    plt.tight_layout()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    plt.savefig('stability.pdf')
    plt.show()