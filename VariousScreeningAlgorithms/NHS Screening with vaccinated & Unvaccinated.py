#Algorithm for the 2023 NHS Cervical Cancer Screening Programme for Anyone With a Cervix Aged 24.5-64

import hpvsim as hpv
import numpy as np


if __name__=="__main__":
    # Define the parameters
    pars = dict(
        n_agents      = 20e3,            # Population size
        start = 2015,                    # Start Year
        n_years       = 110,#70,              # Number of years to simulate
        verbose       = 0,               # Don't print details of the run
        rand_seed     = 2,               # Set a non-default seed
        burnin        = 5,               #Burn in period to stabilise the results
        genotypes     = [16, 18,'hr'],   # Include the two genotypes of greatest general interest
    )

    # Define a series of interventions to screen, triage, assign treatment, and administer treatment
    primary_screen_prob = 0.6
    start_year = 2021
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    prob = 0.6
    # Define eligibility based on screening history, age, and immunity level
    screen_eligible1 = lambda sim: (
        (np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt']))) &
        (sim.people.age >= 24) & (sim.people.age <= 65) &
        (sim.people.sev_imm < 0.5)  # Assuming 'immunity' attribute exists and <0.5 indicates low immunity
    )

    screen_eligible2 = lambda sim: (
        (np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt']))) &
        (sim.people.age >= 24) & (sim.people.age <= 65) &
        (sim.people.sev_imm > 0.5)  # Assuming 'immunity' attribute exists and <0.5 indicates low immunity
    )

    #Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary  = hpv.routine_screening(eligibility=screen_eligible1, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
    cytology   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo = hpv.routine_triage(eligibility=to_colpo, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation  = hpv.treat_num(eligibility=hsils, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Not vaccinated')

    #Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary  = hpv.routine_screening(eligibility=screen_eligible2, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
    cytology   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo = hpv.routine_triage(eligibility=to_colpo, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation  = hpv.treat_num(eligibility=hsils, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo1 = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Vaccinated')

    # Create the sim with and without interventions
    orig_sim = hpv.Sim(pars, label='Baseline') 

    # Run and plot
    msim = hpv.parallel(orig_sim, algo, algo1)
    msim.plot();