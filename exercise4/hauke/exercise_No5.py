"""
Created on Wed Nov 13 15:33:12 2019

@author: dominik
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.rcParams['axes.grid'] = True


def linear(x, b, m):
    return (m * x + b)


if __name__ == "__main__":
    
    file1 = pd.read_csv("result.txt", sep=",")
    file1.columns = ["nParticles", "angle", "timeTree", "timeFull", "error", "interaction"]
    
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax1.set_ylabel('Execution time in s')
    ax1.set_xlabel('Number of objects N')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax1.plot(file1['nParticles'].values[np.argwhere(file1['angle'].values==0.4)], file1['timeTree'].values[np.argwhere(file1['angle'].values==0.4)], marker='+', color='blue', markersize=10, linestyle = 'None', label='raw data - Tree algorithm')
    
    popt_Tree, pcov_Tree = curve_fit(linear, np.log(file1['nParticles'].values[np.argwhere(file1['angle'].values==0.4)])[:,0], np.log(file1['timeTree'].values[np.argwhere(file1['angle'].values==0.4)])[:,0])
    x = np.linspace(np.min(file1['nParticles'].values[np.argwhere(file1['angle'].values==0.4)]), np.max(file1['nParticles'].values[np.argwhere(file1['angle'].values==0.4)]), num=1000)
    y_sim = np.zeros(len(x))
    for t in range(len(x)):
        y_sim[t] = linear(np.log(x[t]), popt_Tree[0], popt_Tree[1])
    
    ax1.plot(x, np.exp(y_sim), 'r-', label='linear regression, m = {} $\pm$ {} [68%]'.format(np.round(popt_Tree[1], 2), np.round(np.sqrt(pcov_Tree[1,1]), 2)))
    
    ax2 = fig.add_subplot(122)
    ax2.set_ylabel('Execution time in s')
    ax2.set_xlabel('Number of objects N')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    ax2.plot(file1['nParticles'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))], file1['timeFull'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))], marker='+', color='blue', markersize=10, linestyle = 'None', label='raw data - Direct summation')
    
    popt_Direct, pcov_Direct = curve_fit(linear, np.log(file1['nParticles'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))])[:,0], np.log(file1['timeFull'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))])[:,0])
    x = np.linspace(np.min(file1['nParticles'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))]), np.max(file1['nParticles'].values[np.argwhere((file1['angle'].values==0.4) & (file1['error'].values!=0.))]), num=1000)
    y_sim = np.zeros(len(x))
    for t in range(len(x)):
        y_sim[t] = linear(np.log(x[t]), popt_Direct[0], popt_Direct[1])
    
    ax2.plot(x, np.exp(y_sim), 'r-', label='linear regression, m = {} $\pm$ {} [68%]'.format(np.round(popt_Direct[1], 2), np.round(np.sqrt(pcov_Direct[1,1]), 2)))
    
    plt.tight_layout()
    ax1.legend()
    ax2.legend()
    plt.savefig('Task_e.pdf')
    plt.show()