import hpvsim as hpv

vx1 = hpv.routine_vx(prob=0.6, start_year=2000, age_range=[9,10], product='bivalent')
# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    location      = 'United Kingdom',
    n_years       = 80,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
)

# Create the sim with and without interventions
orig_sim = hpv.Sim(pars, label='Baseline')
sim1 = hpv.Sim(pars, interventions = vx1, label='With single dose')

# Run and plot
msim = hpv.parallel(orig_sim, sim1)
msim.plot();