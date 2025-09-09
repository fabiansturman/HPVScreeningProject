#The 7 algorithms recommended in the WHO's guidelines for screening and treatment of cervical pre-cancer lesions.

import hpvsim as hpv
import numpy as np
import random
import pandas as pd

# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    n_years       = 90,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18,'hr'],   # Include the two genotypes of greatest general interest
)

# Define a series of interventions to screen, triage, assign treatment, and administer treatment
start_year = 2000
primary_screen_prob = 0.6
triage_screen_prob = 0.9
ablate_prob = 0.9
prob = 0.6
screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt']))

#Algorithm 7. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
hpv_primary7  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary7') #HPV Testing
# Send HPV+ women* for cytology
to_cytology = lambda sim: sim.get_intervention('hpv_primary7').outcomes['positive']
cytology7   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology7')
#Send ASCUS and abnormal cytology results for colpo.
to_colpo7 = lambda sim:list(set(sim.get_intervention('cytology7').outcomes['abnormal'].tolist() + sim.get_intervention('cytology7').outcomes['ascus'].tolist())) #Define who's eligible for colpo
colpo7 = hpv.routine_triage(eligibility=to_colpo7, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
hsils7 = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
ablation7  = hpv.treat_num(eligibility=hsils7, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

algo = hpv.Sim(pars, interventions = [hpv_primary7, cytology7, colpo7, ablation7], label='NHS Cervical Screening Pathway 2023 (3 Year Call Back)')
##### 10 year call back #######
################## 4 Life time screens ####################
screen_eligible4 = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 10 / sim['dt']))

#Algorithm 7. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
hpv_primary74  = hpv.routine_screening(eligibility=screen_eligible4, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary74') #HPV Testing
# Send HPV+ women* for cytology
to_cytology4 = lambda sim: sim.get_intervention('hpv_primary74').outcomes['positive']
cytology74   = hpv.routine_triage(eligibility=to_cytology4, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology74')
#Send ASCUS and abnormal cytology results for colpo.
to_colpo74 = lambda sim:list(set(sim.get_intervention('cytology74').outcomes['abnormal'].tolist() + sim.get_intervention('cytology74').outcomes['ascus'].tolist())) #Define who's eligible for colpo
colpo74 = hpv.routine_triage(eligibility=to_colpo74, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy4') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
hsils74 = lambda sim: sim.get_intervention('colposcopy4').outcomes['hsil'] # Define who's eligible for ablation
ablation74  = hpv.treat_num(eligibility=hsils74, prob=ablate_prob, product='ablation', label='ablation4') # Administer ablation

algo10 = hpv.Sim(pars, interventions = [hpv_primary74, cytology74, colpo74, ablation74], label='10 Year Call Back')

##### 6 year call back #######
screen_eligible6 = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 6 / sim['dt']))
hpv_primary76  = hpv.routine_screening(eligibility=screen_eligible6, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary76') #HPV Testing
# Send HPV+ women* for cytology
to_cytology6 = lambda sim: sim.get_intervention('hpv_primary76').outcomes['positive']
cytology76   = hpv.routine_triage(eligibility=to_cytology6, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology76')
#Send ASCUS and abnormal cytology results for colpo.
to_colpo76 = lambda sim:list(set(sim.get_intervention('cytology76').outcomes['abnormal'].tolist() + sim.get_intervention('cytology76').outcomes['ascus'].tolist())) #Define who's eligible for colpo
colpo76 = hpv.routine_triage(eligibility=to_colpo76, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy6') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
hsils76 = lambda sim: sim.get_intervention('colposcopy6').outcomes['hsil'] # Define who's eligible for ablation
ablation76  = hpv.treat_num(eligibility=hsils76, prob=ablate_prob, product='ablation', label='ablation6') # Administer ablation

algo6 = hpv.Sim(pars, interventions = [hpv_primary76, cytology76, colpo76, ablation76], label='6 Year Call Back')

##### 15 year call back #######
screen_eligible15 = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 15 / sim['dt']))
hpv_primary715  = hpv.routine_screening(eligibility=screen_eligible15, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary715') #HPV Testing
# Send HPV+ women* for cytology
to_cytology15 = lambda sim: sim.get_intervention('hpv_primary715').outcomes['positive']
cytology715   = hpv.routine_triage(eligibility=to_cytology15, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology715')
#Send ASCUS and abnormal cytology results for colpo.
to_colpo715 = lambda sim:list(set(sim.get_intervention('cytology715').outcomes['abnormal'].tolist() + sim.get_intervention('cytology715').outcomes['ascus'].tolist())) #Define who's eligible for colpo
colpo715 = hpv.routine_triage(eligibility=to_colpo715, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy15') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
hsils715 = lambda sim: sim.get_intervention('colposcopy15').outcomes['hsil'] # Define who's eligible for ablation
ablation715  = hpv.treat_num(eligibility=hsils715, prob=ablate_prob, product='ablation', label='ablation15') # Administer ablation

algo15 = hpv.Sim(pars, interventions = [hpv_primary715, cytology715, colpo715, ablation715], label='15 Year Call Back')

# Create the sim with and without interventions
orig_sim = hpv.Sim(pars, label='Baseline')

# Run and plot
#msim = hpv.parallel(algo, algo6, algo10, algo15)
#msim.plot(['cancer_incidence']);

trials_data = []  # Initialize a list to store trial DataFrames

   # Loop to run sim1 30 times with different random seeds
for i in range(20):
    # Generate a random seed for each iteration. You can specify any range for the random seed.
    random_seed = random.randint(1, 10)
    
    # Update the 'rand_seed' parameter with the new random seed
    pars['rand_seed'] = random_seed
    
    # Create and run the simulation with the updated parameters
    algo15 = hpv.Sim(pars, interventions = [hpv_primary715, cytology715, colpo715, ablation715], label='15 Year Call Back')
    #algo = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo7, ablation], label='NHS Cervical Screening Pathway 2023')
    #algo.run()
    algo15.run()
    trial_data = algo15.to_df()  # Assuming this is the correct method to extract data; adjust if necessary
    trial_data['trial_number'] = i + 1
    trials_data.append(trial_data)

# Concatenate and export trial DataFrames as before
all_trials_df = pd.concat(trials_data, ignore_index=True)
all_trials_df.to_excel('algo15_random_seed(20).xlsx', index=False)

print("Exported all trials data with random seeds to sim3_all_trials_random_seed(20).xlsx")
