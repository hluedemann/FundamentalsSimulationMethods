from tree import *
from random import *
import time
from copy import *
from math import sqrt

import matplotlib.pyplot as plt

#
# Create a set of randomly positioned particles
# For convenience we assume all masses to be 1.
# If we have in reality another mass, we can simply
# rescale our answers.
#
nparticles = 1000

def runForceCalculation(nparticles, angle):


    particles = []
    for i in range(nparticles):
        x = random()
        y = random()
        z = random()
        particles.append([x,y,z])

#
# Now create the tree
#
    q=TreeClass(particles)
    q.insertallparticles()
    q.computemultipoles(0)


    print ("starting tree gravity with {} particles and angle {}".format(nparticles, angle))
    t0 = time.time()
    q.allgforces(angle)
    t1 = time.time()
    treegrav_dt = t1-t0
    print ("done in "+str(treegrav_dt)+" seconds\n")

    

    fapprox = deepcopy(q.forces)
    interaction = np.sum(q.interactionCount)/len(q.interactionCount) 
    print("Interaction: ", interaction)
    print("Totale Number Nodes: ", len(q.nodelist))
    print("\n")

    print ("starting N^2 gravity")

    t0 = time.time()

    # Using an angle of 0 results in the exact summation of all forces
    q.allgforces(0)

    t1 = time.time()
    treegrav_dt = t1-t0


    t1 = time.time()
    fullgrav_dt = t1-t0
    print ("done in "+str(fullgrav_dt)+" seconds\n")

    fexact = deepcopy(q.forces)

    errorPerParticle =  np.sum((np.array(fapprox) - np.array(fexact))**2, axis=1)**(1/2) / np.sum((np.array(fexact))**2, axis=1)**(1/2)
    error = np.sum(errorPerParticle)/len(fexact[:,0])

    print("Mean error: ", error)

    # Save data
    f=open('result.txt','a')
    np.savetxt(f,[[nparticles, angle, treegrav_dt, fullgrav_dt, error, interaction]], delimiter=',')
    f.close()

    # Create histogram of errors
    plt.xlabel("Relative absolute error")
    plt.ylabel("Frequency")
    plt.hist(errorPerParticle, bins=int(nparticles/10))
    plt.savefig("errorPlot_{}_{}.pdf".format(nparticles, angle))
    plt.close()


if __name__ == "__main__":
    

    f=open('result.txt','a')
    f.write("nParticles, angle, timeTree, timeFull, error, interaction\n")
    f.close()

    for i in [10, 20, 30]:
        for j in [0.2, 0.5, 0.6]:
            runForceCalculation(i, j)