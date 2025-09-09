import hpvsim as hpv
import matplotlib.pyplot as plt

# Define the parameters
pars = dict(
        n_agents      = 20e3,       # Population size
        location      = 'united kingdom',
        n_years       = 50,         # Number of years to simulate
        verbose       = 0,          # Don't print details of the run
        rand_seed     = 2,          # Set a non-default seed
        genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    )

if __name__=="__main__":
    # Create an empty list to store the simulation results
    results = []
    prob_girls= 0.7
    vxf = hpv.routine_vx(prob=prob_girls, start_year=2015, age_range=[9, 10], sex='f', product='nonavalent')
    orig_sim = hpv.Sim(pars, label=f'Baseline')
    sim_f = hpv.Sim(pars, interventions=vxf, label=f'Girls {prob_girls*100}% Vaccination Coverage')

    # Define a range of values for the prob parameter
    prob_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    # Iterate over different prob values
    '''for prob_value in prob_values:
        '''# Create the intervention with the current prob value
    vxm5 = hpv.routine_vx(prob=0.5, start_year=2015, age_range=[9, 14], sex=['f', 'm'], product='nonavalent')
        # Create the sim with the intervention
    '''sim_m = hpv.Sim(pars, interventions=vxm, label=f'With boys - prob={prob_value}')
    '''
    vxm6 = hpv.routine_vx(prob=0.6, start_year=2015, age_range=[9, 14], sex=['f', 'm'], product='nonavalent')
        # Create the sim with the intervention
    ''' sim_m = hpv.Sim(pars, interventions=vxm, label=f'With boys - prob={prob_value}')'''

    vxm7 = hpv.routine_vx(prob=0.7, start_year=2015, age_range=[9, 14], sex=['f', 'm'], product='nonavalent')
    # Create the sim with the intervention
    sim_m5 = hpv.Sim(pars, interventions=vxm5, label=f'Girls and Boys 50% Vaccination Coverage')
    sim_m6 = hpv.Sim(pars, interventions=vxm6, label=f'Girls and Boys 60% Vaccination Coverage')
    sim_m7 = hpv.Sim(pars, interventions=vxm7, label=f'Girls and Boys 70% Vaccination Coverage')

        # Run the simulation
    msim = hpv.parallel(orig_sim, sim_m5, sim_m6, sim_m7, sim_f)
        # Append the simulation results to the list
    results.append(msim)

    msim.plot(['infections','hpv_incidence', 'cancer_incidence'])


'''
import hpvsim as hpv
import matplotlib.pyplot as plt

# Define a range of values for the prob parameter
prob_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# Create an empty list to store the simulation results
results = []
orig_sim = hpv.Sim(pars, label='Baseline')
# Iterate over different prob values
for prob_value in prob_values:
    # Create the intervention with the current prob value
    vx_intervention = hpv.routine_vx(prob=prob_value, start_year=2015, age_range=[9, 14], sex=['f', 'm'], product='nonavalent')
    
    # Create the sim with the intervention
    sim = hpv.Sim(pars, interventions=vx_intervention, label=f'prob={prob_value}')
    
    # Run the simulation
    sim.run()
    
    # Append the simulation results to the list
    results.append(sim)

# Plot the outcomes for each prob value
    
#hpv.plot.plot_scens(results, ['infections', 'hpv_incidence', 'cancer_incidence'])
#plt.show()

# Run and plot
msim = hpv.parallel(orig_sim, sim)
#msim.plot();
msim.plot(['infections','hpv_incidence', 'cancer_incidence'])
'''