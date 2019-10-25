"""
Created on Thu Oct 24 17:51:34 2019

@author: Dominik Werner Wolf
@version: 1.0
"""
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams['axes.grid'] = True


##############################################################################
## Visualisation function imported from moodle:
##############################################################################

def frame(phi1, phi2, phi1_arr, phi2_arr, t, num):
    #length of the pendula
    L1 = 2.
    L2 = 1.
    L = L1 + L2
    
    #create a new figure
    fig = plt.figure(figsize=[8,8], facecolor='white')
    axes = fig.gca()

    #proportional x and y axis
    axes.axis("equal")
    #no coordinate labels
    axes.set_axis_off()
    #set the plot range in x and y direction
    axes.set_ylim([-(L+0.2),(L+0.2)])
    axes.set_xlim([-(L+0.2),(L+0.2)])
    
    #calculate position of pendula
    pos1 = np.zeros(2)
    pos2 = np.zeros(2)
    
    pos1[0] = L1 * np.sin(phi1)
    pos1[1] = - L1 * np.cos(phi1)
    pos2[0] = pos1[0] + L2 * np.sin(phi2)
    pos2[1] = pos1[1] - L2 * np.cos(phi2)
    
    #plot the two pendula
    c0 = plt.Circle([0.,0.],.07, color='k', zorder=100)
    c1 = plt.Circle(pos1,.05, color='r', zorder=100)
    c2 = plt.Circle(pos2,.05, color='b', zorder=100)
    axes.add_artist(c0)
    axes.add_artist(c1)
    axes.add_artist(c2)
    
    #plot the two limbs
    axes.plot([0., pos1[0]], [0., pos1[1]], color='k', linewidth=2, zorder=50)
    axes.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='k', linewidth=2, zorder=50)
    
    #add the trajectory in the background
    N = len(phi1_arr)
    pos1_arr = np.zeros((N,2))
    pos2_arr = np.zeros((N,2))
    
    pos1_arr[:,0] = L1 * np.sin(phi1_arr)
    pos1_arr[:,1] = - L1 * np.cos(phi1_arr)
    pos2_arr[:,0] = pos1_arr[:,0] + L2 * np.sin(phi2_arr)
    pos2_arr[:,1] = pos1_arr[:,1] - L2 * np.cos(phi2_arr)
    
    axes.plot(pos1_arr[:,0], pos1_arr[:,1], color='r', linewidth = 0.5, zorder=0)
    axes.plot(pos2_arr[:,0], pos2_arr[:,1], color='b', linewidth = 0.5, zorder=0)
    
    #add a time label
    axes.text(0.7, 0.85, "T = %03.2f"%t, transform = fig.transFigure, fontsize=16, horizontalalignment='left', verticalalignment='center')
        
    #save the figure
    fig.savefig("frame_%04d.jpg"%num, dpi=100)
    
    plt.close(fig)
    
##############################################################################


def f1(phi1, phi2, q1, q2, g=1., l1=2., l2=1., m1=0.5, m2=1.0):
    return (q1-l1/l2*np.cos(phi1-phi2)*q2)/(l1**2*(m1+m2*(1-np.cos(phi2-phi2)**2)))
    
def f2(phi1, phi2, q1, q2, g=1., l1=2., l2=1., m1=0.5, m2=1.0):
    return ((m1+m2)*q2-m2*l2/l1*np.cos(phi1-phi2)*q1)/(m2*l2**2*(m1+m2*(1-np.cos(phi1-phi2)**2)))
    
def f3(phi1, phi2, f1, f2, g=1., l1=2., l2=1., m1=0.5, m2=1.0):
    return (-m2*l1*l2*f1*f2*np.sin(phi1-phi2)-g*l1*(m1+m2)*np.sin(phi1))
    
def f4(phi1, phi2, f1, f2, g=1., l1=2., l2=1., m1=0.5, m2=1.0):
    return (m2*l1*l2*f1*f2*np.sin(phi1-phi2)-m2*g*l2*np.sin(phi2))

def energy(phi1, phi2, q1, q2, g=1., l1=2., l2=1., m1=0.5, m2=1.0):
    dt_phi1 = f1(phi1, phi2, q1, q2)
    dt_phi2 = f2(phi1, phi2, q1, q2)
    return (m1/2.*(l1*dt_phi1)**2+m2/2.*((l1*dt_phi1)**2+(l2*dt_phi2)**2+2*l1*l2*dt_phi1*dt_phi2*np.cos(phi1-phi2))+m1*g*l1*(1-np.cos(phi1))+m2*g*(l1*(1-np.cos(phi1))+l2*(1-np.cos(phi2))))

def RK2(steps, dt, AWP):
    result = np.zeros((4, steps+1))
    result[:,0] = AWP
    f_vec1 = np.zeros((4, 1))
    f_vec2 = np.zeros((4, 1))
    for i in range(steps):
        ## gradient at initial fulcrum:
        f_vec1[0] = f1(result[0,i], result[1,i], result[2,i], result[3,i])
        f_vec1[1] = f2(result[0,i], result[1,i], result[2,i], result[3,i])
        f_vec1[2] = f3(result[0,i], result[1,i], f_vec1[0], f_vec1[1])
        f_vec1[3] = f4(result[0,i], result[1,i], f_vec1[0], f_vec1[1])
        ## gradient at final fulcrum:
        f_vec2[0] = f1(result[0,i]+dt*f_vec1[0], result[1,i]+dt*f_vec1[1], result[2,i]+dt*f_vec1[2], result[3,i]+dt*f_vec1[3])
        f_vec2[1] = f2(result[0,i]+dt*f_vec1[0], result[1,i]+dt*f_vec1[1], result[2,i]+dt*f_vec1[2], result[3,i]+dt*f_vec1[3])
        f_vec2[2] = f3(result[0,i]+dt*f_vec1[0], result[1,i]+dt*f_vec1[1], f_vec2[0], f_vec2[1])
        f_vec2[3] = f4(result[0,i]+dt*f_vec1[0], result[1,i]+dt*f_vec1[1], f_vec2[0], f_vec2[1])
        ## state vector at (i+1):
        result[:,i+1] = result[:,i] + 0.5*dt*np.transpose(f_vec1+f_vec2)
        #input("Press Enter to continue...")
    return result

def RK4(steps, dt, AWP):
    result = np.zeros((4, steps+1))
    result[:,0] = AWP
    f_vec1 = np.zeros((4, 1))
    f_vec2 = np.zeros((4, 1))
    f_vec3 = np.zeros((4, 1))
    f_vec4 = np.zeros((4, 1))
    for i in range(steps):
        ## gradient at initial fulcrum:
        f_vec1[0] = f1(result[0,i], result[1,i], result[2,i], result[3,i])
        f_vec1[1] = f2(result[0,i], result[1,i], result[2,i], result[3,i])
        f_vec1[2] = f3(result[0,i], result[1,i], f_vec1[0], f_vec1[1])
        f_vec1[3] = f4(result[0,i], result[1,i], f_vec1[0], f_vec1[1])
        ## gradient at second fulcrum:
        f_vec2[0] = f1(result[0,i]+dt*f_vec1[0]/2., result[1,i]+dt*f_vec1[1]/2., result[2,i]+dt*f_vec1[2]/2., result[3,i]+dt*f_vec1[3]/2.)
        f_vec2[1] = f2(result[0,i]+dt*f_vec1[0]/2., result[1,i]+dt*f_vec1[1]/2., result[2,i]+dt*f_vec1[2]/2., result[3,i]+dt*f_vec1[3]/2.)
        f_vec2[2] = f3(result[0,i]+dt*f_vec1[0]/2., result[1,i]+dt*f_vec1[1]/2., f_vec2[0], f_vec2[1])
        f_vec2[3] = f4(result[0,i]+dt*f_vec1[0]/2., result[1,i]+dt*f_vec1[1]/2., f_vec2[0], f_vec2[1])
        ## gradient at third fulcrum:
        f_vec3[0] = f1(result[0,i]+dt*f_vec2[0]/2., result[1,i]+dt*f_vec2[1]/2., result[2,i]+dt*f_vec2[2]/2., result[3,i]+dt*f_vec2[3]/2.)
        f_vec3[1] = f2(result[0,i]+dt*f_vec2[0]/2., result[1,i]+dt*f_vec2[1]/2., result[2,i]+dt*f_vec2[2]/2., result[3,i]+dt*f_vec2[3]/2.)
        f_vec3[2] = f3(result[0,i]+dt*f_vec2[0]/2., result[1,i]+dt*f_vec2[1]/2., f_vec3[0], f_vec3[1])
        f_vec3[3] = f4(result[0,i]+dt*f_vec2[0]/2., result[1,i]+dt*f_vec2[1]/2., f_vec3[0], f_vec3[1])
        ## gradient at fourth fulcrum:
        f_vec4[0] = f1(result[0,i]+dt*f_vec3[0], result[1,i]+dt*f_vec3[1], result[2,i]+dt*f_vec3[2], result[3,i]+dt*f_vec3[3])
        f_vec4[1] = f2(result[0,i]+dt*f_vec3[0], result[1,i]+dt*f_vec3[1], result[2,i]+dt*f_vec3[2], result[3,i]+dt*f_vec3[3])
        f_vec4[2] = f3(result[0,i]+dt*f_vec3[0], result[1,i]+dt*f_vec3[1], f_vec4[0], f_vec4[1])
        f_vec4[3] = f4(result[0,i]+dt*f_vec3[0], result[1,i]+dt*f_vec3[1], f_vec4[0], f_vec4[1])
        ## state vector at (i+1):
        result[:,i+1] = result[:,i] + 1./6.*dt*np.transpose(f_vec1+2*f_vec2+2*f_vec3+f_vec4)
        #input("Press Enter to continue...")
    return result


if __name__ == "__main__":
    
    ## starting point of trajectory:
    AWP = np.array([50.*np.pi/180., -120.*np.pi/180., 0., 0.])
    
    steps = 2000
    dt = 0.05
    
    method = input("Please select the method: RK2 or RK4:") 
    
    if (method == 'RK2'):
        results = RK2(steps, dt, AWP)
        
        energy_error = np.zeros(len(results[0,:]))
        for i in range(len(results[0,:])):
            energy_error[i] = (energy(results[0,i], results[1,i], results[2,i], results[3,i]) - energy(results[0,0], results[1,0], results[2,0], results[3,0]))/energy(results[0,0], results[1,0], results[2,0], results[3,0])
        
        time_axis = np.linspace(0.0, steps*dt, num=(steps+1))
        
        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(231)
        ax1.set_ylabel('$\phi_{1}$ in degree units')
        ax1.set_xlabel('time in arbitrary units')
        
        ax1.plot(time_axis, results[0,:]*180./np.pi, 'k.')
        
        ax2 = fig.add_subplot(232)
        ax2.set_ylabel('$\phi_{2}$ in degree units')
        ax2.set_xlabel('time in arbitrary units')
        
        ax2.plot(time_axis, results[1,:]*180./np.pi, 'b.')
        
        ax3 = fig.add_subplot(233)
        ax3.set_ylabel('percentage of relative energy error in %')
        ax3.set_xlabel('time in arbitrary units')
        
        ax3.plot(time_axis, energy_error*100., 'r-')
        ax3.plot(time_axis, np.zeros(len(time_axis)), 'g-')
        
        ax4 = fig.add_subplot(234)
        ax4.set_ylabel('$\dot{\phi}_{1}$ in degree units per time unit')
        ax4.set_xlabel('time in arbitrary units')
        
        ax4.plot(time_axis, results[2,:]*180./np.pi, 'k.')
        
        ax5 = fig.add_subplot(235)
        ax5.set_ylabel('$\dot{\phi}_{2}$ in degree units per time unit')
        ax5.set_xlabel('time in arbitrary units')
        
        ax5.plot(time_axis, results[3,:]*180./np.pi, 'b.')
        
        plt.tight_layout()
        plt.savefig('RK2_dominik.pdf')
        plt.show()
        
        for i in np.arange(steps+1):
            frame(results[0,i], results[1,i], results[0,:][:i+1], results[1,:][:i+1], i*0.05, i)
            
        os.system("ffmpeg -r 24 -i frame_%04d.jpg RK2_movie.mp4")
        
    elif (method == 'RK4'):
        results = RK4(steps, dt, AWP)
        
        energy_error = np.zeros(len(results[0,:]))
        for i in range(len(results[0,:])):
            energy_error[i] = (energy(results[0,i], results[1,i], results[2,i], results[3,i]) - energy(results[0,0], results[1,0], results[2,0], results[3,0]))/energy(results[0,0], results[1,0], results[2,0], results[3,0])
        
        time_axis = np.linspace(0.0, steps*dt, num=(steps+1))
        
        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(231)
        ax1.set_ylabel('$\phi_{1}$ in degree units')
        ax1.set_xlabel('time in arbitrary units')
        
        ax1.plot(time_axis, results[0,:]*180./np.pi, 'k.')
        
        ax2 = fig.add_subplot(232)
        ax2.set_ylabel('$\phi_{2}$ in degree units')
        ax2.set_xlabel('time in arbitrary units')
        
        ax2.plot(time_axis, results[1,:]*180./np.pi, 'b.')
        
        ax3 = fig.add_subplot(233)
        ax3.set_ylabel('percentage of relative energy error in %')
        ax3.set_xlabel('time in arbitrary units')
        
        ax3.plot(time_axis, energy_error*100., 'r-')
        ax3.plot(time_axis, np.zeros(len(time_axis)), 'g-')
        
        ax4 = fig.add_subplot(234)
        ax4.set_ylabel('$\dot{\phi}_{1}$ in degree units per time unit')
        ax4.set_xlabel('time in arbitrary units')
        
        ax4.plot(time_axis, results[2,:]*180./np.pi, 'k.')
        
        ax5 = fig.add_subplot(235)
        ax5.set_ylabel('$\dot{\phi}_{2}$ in degree units per time unit')
        ax5.set_xlabel('time in arbitrary units')
        
        ax5.plot(time_axis, results[3,:]*180./np.pi, 'b.')
        
        plt.tight_layout()
        plt.savefig('RK4_dominik.pdf')
        plt.show()
        
        for i in np.arange(steps+1):
            frame(results[0,i], results[1,i], results[0,:][:i+1], results[1,:][:i+1], i*0.05, i)
            
        os.system("ffmpeg -r 24 -i frame_%04d.jpg RK4_movie.mp4")
        
    else:
        print('Selected method not available!')