    #The 7 algorithms recommended in the WHO's guidelines for screening and treatment of cervical pre-cancer lesions.

import hpvsim as hpv
import numpy as np

if __name__=="__main__":

    # Define the parameters
    pars = dict(
        n_agents      = 20e3,       # Population size
        location      = 'united kingdom',
        n_years       = 70,#35,         # Number of years to simulate
        verbose       = 0,          # Don't print details of the run
        #rand_seed     = 2,          # Set a non-default seed
        genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
    )

    # Define a series of interventions to screen, triage, assign treatment, and administer treatment
    start_year = 2015
    primary_screen_prob = 0.2
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    prob = 0.6
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 5 / sim['dt']))

    screen  = hpv.routine_screening(eligibility=screen_eligible, start_year=2015, prob=primary_screen_prob, product='via', label='screen') # Routine screening
    to_triage   = lambda sim: sim.get_intervention('screen').outcomes['positive'] # Define who's eligible for triage
    triage      = hpv.routine_triage(eligibility=to_triage, prob=triage_screen_prob, product='hpv', label='triage') # Triage people
    to_treat    = lambda sim: sim.get_intervention('triage').outcomes['positive'] # Define who's eligible to be assigned treatment
    assign_tx   = hpv.routine_triage(eligibility=to_treat, prob=prob, product='tx_assigner', label='assign_tx') # Assign treatment
    to_ablate   = lambda sim: sim.get_intervention('assign_tx').outcomes['ablation'] # Define who's eligible for ablation treatment
    ablation    = hpv.treat_num(eligibility=to_ablate, prob=ablate_prob, product='ablation') # Administer ablation
    to_excise   = lambda sim: sim.get_intervention('assign_tx').outcomes['excision'] # Define who's eligible for excision
    excision    = hpv.treat_delay(eligibility=to_excise, prob=prob, product='excision') # Administer excision

    ##Algorithm 1: Visual Inspection with acetic acid (VIA) as the primary screening test, followed by treatment.

    via_primary  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='via', label='via_primary') #Visual Inspection with acetic acid (VIA)
    via_positive = lambda sim: sim.get_intervention('via_primary').outcomes['positive'] # Define who's eligible for ablation
    ablation1   = hpv.treat_num(eligibility=via_positive, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo1 = hpv.Sim(pars, interventions = [via_primary, ablation1], label='Algorithm 1')

    ##Algorithm 2: HPV testing, then immediate ablation for anyone eligible.
    hpv_primary  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
    hpv_positive = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive'] # Define who's eligible for ablation
    ablation2   = hpv.treat_num(eligibility=hpv_positive, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo2 = hpv.Sim(pars, interventions = [hpv_primary, ablation2], label='Algorithm 2')

    ##Algorithm 3: Cytology testing, triage ASCUS results with HPV, triage all HPV+ and abnormal cytology results with colposcopy/biopsy, then ablation for all HSILs.
    cytology  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='lbc', label='cytology') #Cytology testing
    ascus = lambda sim: sim.get_intervention('cytology').outcomes['ascus'] # Define who's eligible for triage, ASCUS = Atypical Squamous Cells of Undetermined Significance.
    hpv_triage = hpv.routine_triage(eligibility=ascus, prob=triage_screen_prob, product='hpv', annual_prob=False, label='hpv triage') # Triage people
    #Send abnormal cytology results, plus ASCUS that were HPV+, for colpo.
    to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('hpv triage').outcomes['positive'].tolist())) #Define who's eligible for colpo
    colpo = hpv.routine_triage(eligibility=to_colpo, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with ablation
    hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation3   = hpv.treat_num(eligibility=hsils, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo3 = hpv.Sim(pars, interventions = [cytology, hpv_triage, colpo, ablation3], label='Algorithm 3')

    ##Algorithm 4: HPV DNA as the primary screening test, followed by HPV16/18 triage
    #(when already part of the HPV test), followed by treatment,
    #and using VIA triage for those who screen negative for HPV16/18
    hpv_primary4  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='hpv_type', label='hpv_primary4') #HPV DNA testing
    #Those who test + for OHR types are triaged with VIA
    pos_ohr = lambda sim: sim.get_intervention('hpv_primary4').outcomes['positive_ohr'] # Define who's eligible for triage
    via_triage = hpv.routine_triage(eligibility=pos_ohr, prob=triage_screen_prob, product='via', annual_prob=False, label='via triage') # Triage people

    #Determine ablation eligibility for people with 16/18, plus those who test positive from VIA
    to_assign = lambda sim:list(set(sim.get_intervention('hpv_primary4').outcomes['positive_1618'].tolist() + sim.get_intervention('via triage').outcomes['positive'].tolist())) #Define who's eligible for triage
    tx_assigner = hpv.routine_triage(eligibility=to_assign, prob=triage_screen_prob, product='tx_assigner', annual_prob=False, label='tx_assigner') #Send people to triage
    #Ablate anyone eligible for ablation
    to_ablate4 = lambda sim: sim.get_intervention('tx_assigner').outcomes['ablation'] # Define who's eligible for ablation
    ablation4   = hpv.treat_num(eligibility=to_ablate4, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation

    algo4 = hpv.Sim(pars, interventions = [hpv_primary4, via_triage, tx_assigner, ablation4], label='Algorithm 4')

    # Create the sim with and without interventions
    orig_sim = hpv.Sim(pars, label='Baseline')
    sim = hpv.Sim(pars, interventions = [screen, triage, assign_tx, ablation], label='With screen & treat')
    sim1 = hpv.Sim(pars, interventions = [via_primary])

    # Run and plot
    msim = hpv.MultiSim([orig_sim, algo1, algo2, algo3, algo4])
    msim.run(n_runs=5) #this only seems to work if i do a multisim with only 1 sim when i set it up, but i suppose thats fine 
    #msim.mean()
    msim.plot()

    for sim in msim.sims:
        sim.brief()