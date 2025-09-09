import hpvsim as hpv

import pandas as pd
my_vaccination_data = pd.DataFrame({'name':'new_vx', 'genotype':['hpv16','hpv18','hi5','ohr','hr','lr'],'rel_imm':[.2,.3,.1,.3,.4,.4]})
my_vaccination = hpv.vx(my_vaccination_data)

def custom_vx(sim):
    if sim.yearvec[sim.t] == 2000:
        target_group = (sim.people.age>9) * (sim.people.age<14)
        sim.people.peak_imm[0, target_group] = 1

pars = dict(
    location = 'United Kingdom', # Use population characteristics for Japan
    n_agents = 10e3, # Have 50,000 people total in the population
    start = 1980, # Start the simulation in 1980
    n_years = 50, # Run the simulation for 50 years
    burnin = 10, # Discard the first 20 years as burnin period
    verbose = 0, # Do not print any output
)

# Running with multisims -- see Tutorial 3
s1 = hpv.Sim(pars, label='Default')
s2 = hpv.Sim(pars, interventions=custom_vx, label='Custom vaccination')
hpv.parallel(s1, s2).plot(['hpv_incidence', 'cancer_incidence']);