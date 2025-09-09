import hpvsim as hpv
import matplotlib.pyplot as plt
import numpy as np

# Define a range of probabilities for sensitivity analysis
probabilities = np.linspace(0.2, 0.8, 4)

# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    n_years       = 80,#50,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
)

if __name__=="__main__":
    # Create the sim with and without interventions
    orig_sim = hpv.Sim(pars, label='Baseline')

    # Run simulations for different probabilities
    simulation_results = []

    for prob in probabilities:
        vx = hpv.routine_vx(prob=prob, start_year=2015, age_range=[9, 10], sex='f', product='nonavalent')
        label = 'With {:2f} % vaccination Uptake'.format(prob*100) 
        sim = hpv.Sim(pars, interventions=vx, label=label)
        simulation_results.append(sim)

    # Plot the results
    msim = hpv.parallel( orig_sim, *simulation_results)
    msim.plot(['hpv_incidence', 'cancer_incidence'])

