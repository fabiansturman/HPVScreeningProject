import hpvsim as hpv # Import HPVsim 

import numpy as np
import matplotlib.pyplot as plt


#TODO: maybe worth adding a contact layer for 'o'ne off relationships as well as (m) and (c) - but for now I can stick with (m) and (c)


sophie_married_matrix = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [10, 0, 0, 0.08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [15, 0, 0, 0.08, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [20, 0, 0, 0, 0, 0.6, 2, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0],
    [25, 0, 0, 0, 0, 0.6, 1, 2, 0.4, 0.1, 0, 0, 0, 0, 0, 0, 0],
    [30, 0, 0, 0, 0, 0.5, 0.5, 2, 1, 0.5, 0.1, 0, 0, 0, 0, 0, 0],
    [35, 0, 0, 0, 0, 1, 0.5, 1, 2, 1, 0.5, 0.2, 0, 0, 0, 0, 0],
    [40, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0],
    [45, 0, 0, 0, 0, 0.1, 1, 2, 2, 2, 1, 0.5, 0.2, 0.08, 0, 0, 0],
    [50, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 2, 2, 0.5, 0.2, 0.05, 0, 0],
    [55, 0, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 3, 2, 1, 0.3, 0.1, 0.1],
    [60, 0, 0, 0, 0, 0, 0, 0.1, 0.5, 1, 2, 3, 3, 2, 0.5, 0.3, 0.1],
    [65, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 2, 3, 3, 2, 1, 0.2],
    [70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 3, 2, 1],
    [75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3],
]
sophie_married_matrix = np.array(sophie_married_matrix)

sophie_casual_matrix = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
       [5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
       [10,0,0,0,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0],
       [15,0,0,1,2,3,2,1,0.5,0,0,0,0,0,0,0,0],
       [20,0,0,0.15,2,3,2,2,1,0.15,0,0,0,0,0,0,0],
       [25,0,0,0.15,0.25,1,2,2,1,1,0,0,0,0,0,0,0],
       [30,0,0,0,0,0.5,0.5,2,1,0.15,0,0,0,0,0,0,0],
       [35,0,0,0,0,1,0.5,1,2,1,0.5,0,0,0,0,0,0],
       [40,0,0,0,0,1,1,1,1,1,0.5,0.25,0,0,0,0,0],
       [45,0,0,0,0,0.15,1,2,2,2,1,0.5,0.2,0.1,0,0,0],
       [50,0,0,0,0,0,0.15,1,2,3,2,2,0.5,0.2,0.05,0,0],
       [55,0,0,0,0,0,0,0.15,1,2,3,3,2,1,0.25,0.1,0.1],
       [60,0,0,0,0,0,0,0.15,0.15,1,2,3,3,2,0.5,0.25,0.1],
       [65,0,0,0,0,0,0,0,0,0,1,1,2,2,1,0.5,0],
       [70,0,0,0,0,0,0,0,0,0,0,0,0,0.8,1,0.7,0.5],
       [75,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.25]
]
sophie_casual_matrix = np.array(sophie_casual_matrix)


pars = dict(n_agents= 10e3,
                start=1980, end=2017, dt=0.25,
                location='united kingdom', 
                rand_seed=0,
                #verbose=-1,
                debut=dict(f=dict(dist='normal', par1=16.0, par2=3.1), m=dict(dist='normal', par1=16.0, par2=4.1)),
                mixing = {'m':sophie_married_matrix,
                          'c':sophie_casual_matrix},
                condoms = dict(m=0.01, c=0.2),  #condom usage in (m)arried and (c)asual relationships
                #NOTE: I think having too many young cancers is probably caused by the debut age being too low <-or it could just be the whole data of that being in shambles, now we have nowhere near enough younger cancers
                )


sim = hpv.Sim(pars) # Make sim

sim.initialize()

#print(type(sim['mixing']['m']))



print(sim.people.layer_keys())

sim.run()		 # Run the simulation
sim.plot()
