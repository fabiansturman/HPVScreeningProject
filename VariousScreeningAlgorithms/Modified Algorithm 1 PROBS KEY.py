import hpvsim as hpv
import numpy as np
import random
import pandas as pd

# Define the parameters
pars = dict(
    n_agents      = 20e3,            # Population size
    start         = 2000,            # Start Year
    n_years       = 90,              # Number of years to simulate
    verbose       = 0,               # Don't print details of the run
    rand_seed     = 2,               # Set a non-default seed
    burnin        = 5,               # Burn-in period to stabilise the results
    genotypes     = [16, 18, 'hr'],  # Include the two genotypes of greatest general interest
)

# Define a series of interventions to screen, triage, assign treatment, and administer treatment
primary_screen_prob = 0.6
start_year = 2021
triage_screen_prob = 0.8
ablate_prob = 0.9

#NHS probabilities
nhsprimary_screen_prob = 0.6
nhsstart_year = 2021
nhstriage_screen_prob = 0.9
nhsablate_prob = 0.9

# Adjusted eligibility criteria to include recall periods based on screening results and vaccination status
def screen_eligible_with_recall(sim, unvaccinated_recall_years=3, vaccinated_recall_years=5.5):
    # Determine the recall period for each individual based on vaccination status
    # This requires element-wise operation to handle arrays correctly
    recall_periods = np.where(sim.people.vaccinated, unvaccinated_recall_years, vaccinated_recall_years)
    
    # Calculate the time until next screening for each individual
    time_until_rescreen = recall_periods / sim['dt']
    
    # Determine if it is time for re-screen based on individual recall periods
    is_time_for_rescreen = np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + time_until_rescreen))
    
    # Determine eligibility based on time for re-screen, age, and possibly immunity level
    is_eligible = is_time_for_rescreen & (sim.people.age >= 24) & (sim.people.age <= 65)
    
    return is_eligible

if __name__=="__main__":
    #Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
    cytology   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo = hpv.routine_triage(eligibility=to_colpo, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation  = hpv.treat_num(eligibility=hsils, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Modified Algorithm 2')

    ###### NHS Algorithm For Comparison #######
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt'])) & (sim.people.age >= 20) & (sim.people.age <= 65)

    #Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=nhsprimary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
    cytology   = hpv.routine_triage(eligibility=to_cytology, prob=nhstriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo7 = hpv.routine_triage(eligibility=to_colpo, prob=nhstriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation  = hpv.treat_num(eligibility=hsils, prob=nhsablate_prob, product='ablation', label='ablation') # Administer ablation

    algo1 = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo7, ablation], label='NHS Cervical Screening Pathway 2023')

    # Create and run the simulation
    orig_sim = hpv.Sim(pars, label='No Screening') 

    # Run the simulation
    msim = hpv.parallel(orig_sim, algo1, algo)

    

    # Loop to run sim1 30 times with different random seeds
    for i in range(4):#20):
        # Generate a random seed for each iteration. You can specify any range for the random seed.
        random_seed = random.randint(1, 10)
        
        # Update the 'rand_seed' parameter with the new random seed
        pars['rand_seed'] = random_seed
        
        # Create and run the simulation with the updated parameters
        offset = 0 #adjustable so that we can add to an existing colection of trial data easily 
        run_name = f'Modified Algorithm 2 (MAYBE 1 CHECK) trial {i+offset}' 
        algo = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label=run_name)
        algo.run()

        #Save the data from this simulation - saving into seperate excel files,
        trial_data = algo.to_df()  # Assuming this is the correct method to extract data; adjust if necessary
        trial_data.to_excel(f"{run_name}.xlsx", index=False)

   

    print("Exported all trials data with random seeds to sim3_all_trials_random_seed(20).xlsx")

    # Plot the results: 
    msim.plot('cancer_incidence')