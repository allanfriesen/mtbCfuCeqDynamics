import KineticModel as km
import StochasticSimulation as ssa
import numpy as np
import matplotlib.pyplot as plt
from math import *

nTraj = 100

model = km.KineticModel()
model.add_species( 'viable' )
model.add_species( 'dead' )
model.add_reaction( 1, [ 1, 0] )
model.add_reaction( 1, [-1,+1] )
model.add_reaction( 2, [ 0,-1] )

#Uncomment the following three lines for low early rates (zero early death)
#model.add_stage( [0.4223237, 0.,        0.036], 21. )
#model.add_stage( [1.0657719, 0.858831,  0.036], 28. )
#model.add_stage( [0.809613,  0.8654658, 0.036], 500.)

#Uncomment the following three lines for high early rates (high death, maximum reasonable replication rate)
model.add_stage( [0.9976605, 0.5753368, 0.036], 21. )
model.add_stage( [0.9976605, 0.7907196, 0.036], 28. )
model.add_stage( [0.809613,  0.8654658, 0.036], 500.)
model.printSpecies()

#create a simulation object that defines the parameters of the simulation
simulation = ssa.StochasticSimulation( model, nTraj=nTraj, dtRecip = 10, tSim = 77, initialPopulations=[[1,0]]*nTraj, mode='normal' )

#run the simulation and store the results in a dictionary called 'results'
results=simulation.run_simulation()

#numpy array holding the trajectories. They are returned as a python list, but 
#converting to numpy array makes it easier to compute CEQs and CFU/CEQ ratio
trajectories = np.array(results['trajectories'])

#numpy array holding the times
times = np.array(results['times'])

#Get trajetories of CFUs (b) and dead chromosomes (d)
b = trajectories[0]
d = trajectories[1]

#calculate CEQs
q = b + d

#Replace zeros of CEQs, for visualization
for sim in range( np.size(b,axis=0) ):
    for t in range( np.size( b[sim] ) ):
        if( q[sim][t] < 1 ):
            q[sim][t] = 0.1

#compute CFU/CEQ ratio
z = b / q

#plot trajectories of CFUs
fig,ax = plt.subplots()
for sim in range( len( b ) ):
    ax.plot(times,b[sim])
ax.set_yscale( 'log' )
ax.set_ylabel( 'CFUs' )
ax.set_xlabel( 'days since infection' )
plt.show()

#plot trajectories of CEQs
fig,ax = plt.subplots()
for sim in range( len( q ) ):
    ax.plot(times,q[sim])
ax.set_yscale( 'log' )
ax.set_ylabel( 'CEQs' )
ax.set_xlabel( 'days since infection' )
plt.show()

#plot trajectories of CFUs/CEQs
fig,ax = plt.subplots()
for sim in range( len( z ) ):
    ax.plot(times,z[sim])
ax.set_yscale( 'log' )
ax.set_ylabel( 'CFUs/CEQs' )
ax.set_xlabel( 'days since infection' )
plt.show()

