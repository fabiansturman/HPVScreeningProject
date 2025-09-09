import hpvsim as hpv

# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    n_years       = 65,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    vaccine_pars  = dict(),     # Initialize the vaccine parameters dictionary
)

if __name__=="__main__":
    # Define the bivalent vaccine intervention
    vx1 = hpv.routine_vx(prob=0.6, start_year=2015, age_range=[9, 10], product='bivalent')

    # Customize immunity parameters for the bivalent vaccine
    pars['vaccine_pars']['bivalent'] = dict(
        rel_immunity = 0.1,          # Relative immunity conferred by the vaccine
        max_dur      = 20,            # Maximum duration of immunity in years
        wane_immunity= 'exp',         # Type of waning immunity ('exp' for exponential)
        eff_duration = 5,             # Effective duration of immunity in years
    )

    # Create the simulation with and without interventions
    orig_sim = hpv.Sim(pars, label='Baseline')
    sim1 = hpv.Sim(pars, interventions=vx1, label='With Single Dose')

    # Run and plot
    msim = hpv.parallel(orig_sim, sim1)
    msim.plot();
