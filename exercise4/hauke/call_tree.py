from tree import *
from random import *
import time
from copy import *
from math import sqrt

import matplotlib.pyplot as plt



#
# Function to run the force calculation for the specified number of particles
#

def runTree(q, nparticles, angle, printProgress):

    
    print(50*"#")
    print ("starting tree gravity with {} particles and angle {}".format(nparticles, angle))
    t0 = time.time()
    q.allgforces(angle, printProgress)
    t1 = time.time()
    treegrav_dt = t1-t0
    print ("done in "+str(treegrav_dt)+" seconds\n")

    

    fapprox = deepcopy(q.forces)
    interaction = np.sum(q.interactionCount)/len(q.interactionCount) 
    print("Mean node interaction of particle: ", interaction)
    print("Totale number nodes: ", len(q.nodelist))
    print("\n")

    return treegrav_dt, fapprox, interaction


def runExact(q, nparticles, printProgress):

    print(50*"#")
    print ("starting N^2 gravity with {} particles".format(nparticles))

    t0 = time.time()

    # Using an angle of 0 results in the exact summation of all forces
    q.allgforces(0, printProgress)

    t1 = time.time()

    fullgrav_dt = t1-t0
    print ("done in "+str(fullgrav_dt)+" seconds\n")

    fexact = deepcopy(q.forces)

    print(50*"#")

    return fullgrav_dt, fexact






if __name__ == "__main__":
    

    #######
    ## Run the simulation with a low number of particels as an example
    #######

    nparticles = 500
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

    # Run the exact simulation
    timeExact, forceExact = runExact(q, nparticles, False)

    # Run the tree simulation
    timeTree, forceTree, interaction = runTree(q, nparticles, 0.8, False)

    # Calculate the mean error
    errorPerParticle =  np.sum((np.array(forceTree) - np.array(forceExact))**2, axis=1)**(1/2) / np.sum((np.array(forceExact))**2, axis=1)**(1/2)
    error = np.sum(errorPerParticle)/len(forceExact[:,0])

    print("Mean error: ", error)
    print(50*"#")
    print("\n")



    ##########
    # Run the simulation for multiple parameters
    # Note: This takes very long
    ##########

    if False: # Set this to run multiple simulations

        # Create the header to save the data
        outputPath = "result.txt"
        f=open(outputPath,'a')
        f.write("nParticles,angle,timeTree,timeFull,error,interaction\n")
        f.close()

        # Tarameters to simulate
        angles = np.array([0.2, 0.4, 0.8])
        numberParticles = np.array([5000, 10000, 20000, 40000])


        for nparticles in numberParticles:

            #
            # Create a set of randomly positioned particles
            # For convenience we assume all masses to be 1.
            # If we have in reality another mass, we can simply
            # rescale our answers.
            #

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

            timeExact, forceExact = runExact(q, nparticles, True)

            for angle in angles:

                timeTree, forceTree, interaction = runTree(q, nparticles, angle, True)
                
                errorPerParticle =  np.sum((np.array(forceTree) - np.array(forceExact))**2, axis=1)**(1/2) / np.sum((np.array(forceExact))**2, axis=1)**(1/2)
                error = np.sum(errorPerParticle)/len(forceExact[:,0])

                print("Mean error: ", error)
                print(50*"#")
                print("\n")

                # Save data
                f=open(outputPath,'a')
                np.savetxt(f,[[nparticles, angle, timeTree, timeExact, error, interaction]], delimiter=',')
                f.close()

                # Create histogram of errors
                plt.xlabel("Relative absolute error")
                plt.ylabel("Frequency")
                plt.hist(errorPerParticle, bins=int(nparticles/200))
                plt.savefig("errorPlot_{}_{}.pdf".format(nparticles, angle))
                plt.close()
                
                