#The 7 algorithms recommended in the WHO's guidelines for screening and treatment of cervical pre-cancer lesions.

#Currently the code simulates the current NHS cervical screening algortithm 
#It should have the eligibility criteria as: individuals who have either:
#not been screened or have not been screened in the last 3 years, and they have been vaccinated.
#This should allow me to see the effect of the screening pathway on vaccinated vs unvaccinated individuals.
#Once I am able to run the simulation on vaccinated and unvaccinated populations and plot them together, I will be able to compare.
#I will devise improved cervical screening algorithms which takes into account the vaccinated cohorts of girls who are
#growing up and reaching vaccination age.


import hpvsim as hpv
import numpy as np
import pandas as pd
import sciris as sc

if __name__=="__main__":
    # Define the parameters
    pars = dict(
        n_agents      = 20e3,       # Population size
        n_years       = 35,         # Number of years to simulate
        verbose       = 0,          # Don't print details of the run
        rand_seed     = 2,          # Set a non-default seed
        genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    )

    # Define a series of interventions to screen, triage, assign treatment, and administer treatment
    start_year = 2015
    primary_screen_prob = 0.2
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    prob = 0.6

    #this eligibility only screens people who have not been vaccinated (we can use this for people who will continue to use the NHS pathway9)
    sim = hpv.Sim(pars, label='MySimulation')

    # Define your eligibility criteria
    vaccinated_people_eligibility = lambda sim: sim.people.vaccinated
    non_vaccinated_people_eligibility = lambda sim: ~sim.people.vaccinated

    ######Screening for vaccinated people########
    #Algorithm 7. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary7  = hpv.routine_screening(eligibility=vaccinated_people_eligibility, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary7') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary7').outcomes['positive']
    cytology7   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology7')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo7 = lambda sim:list(set(sim.get_intervention('cytology7').outcomes['abnormal'].tolist() + sim.get_intervention('cytology7').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo7 = hpv.routine_triage(eligibility=to_colpo7, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils7 = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation7  = hpv.treat_num(eligibility=hsils7, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo1 = hpv.Sim(pars, interventions = [hpv_primary7, cytology7, colpo7, ablation7], label='Vaccinated only NHS Cervical Screening Pathway 2023')

    ######Screening for non vaccinated people########

    #Algorithm 7. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary  = hpv.routine_screening(eligibility=non_vaccinated_people_eligibility, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
    cytology   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo = hpv.routine_triage(eligibility=to_colpo, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation= hpv.treat_num(eligibility=hsils, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo2 = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Non Vaccinated only NHS Cervical Screening Pathway 2023')


    # Create the sim with and without interventions
    orig_sim = hpv.Sim(pars, label='Baseline')

    # Run and plot
    msim = hpv.parallel(orig_sim, algo1, algo2)
    msim.plot();