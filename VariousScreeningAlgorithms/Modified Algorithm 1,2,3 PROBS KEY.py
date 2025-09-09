import hpvsim as hpv
import numpy as np
import pylab as plt
import pandas as pd

# Define the parameters
pars = dict(
    n_agents      = 20e3,            # Population size
    start         = 2015,            # Start Year
    n_years       = 70,              # Number of years to simulate
    verbose       = 0,               # Don't print details of the run
    rand_seed     = 2,               # Set a non-default seed
    burnin        = 5,               # Burn-in period to stabilise the results
    genotypes     = [16, 18, 'hr'],  # Include the two genotypes of greatest general interest
)

# Define a series of interventions to screen, triage, assign treatment, and administer treatment
#algorithm1 pars
aprimary_screen_prob = 0.6
start_year = 2021
atriage_screen_prob = 0.74
aablate_prob = 0.9

#algorithm2 pars
bprimary_screen_prob = 0.6
start_year = 2021
btriage_screen_prob = 0.74#0.8
bablate_prob = 0.9

#algorithm3 pars
cprimary_screen_prob = 0.6#0.55
start_year = 2021
ctriage_screen_prob = 0.74#0.85
cablate_prob = 0.9

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

#Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
hpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=start_year, prob=aprimary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
# Send HPV+ women* for cytology
to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
cytology   = hpv.routine_triage(eligibility=to_cytology, prob=atriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
#Send ASCUS and abnormal cytology results for colpo.
to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
colpo = hpv.routine_triage(eligibility=to_colpo, prob=atriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
ablation  = hpv.treat_num(eligibility=hsils, prob=aablate_prob, product='ablation', label='ablation') # Administer ablation

algo1 = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Modified Algorithm 1')

#Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
bhpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=start_year, prob=bprimary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
# Send HPV+ women* for cytology
bto_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
bcytology   = hpv.routine_triage(eligibility=bto_cytology, prob=btriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
#Send ASCUS and abnormal cytology results for colpo.
bto_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
bcolpo = hpv.routine_triage(eligibility=bto_colpo, prob=btriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
bhsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
bablation  = hpv.treat_num(eligibility=bhsils, prob=bablate_prob, product='ablation', label='ablation') # Administer ablation

algo2 = hpv.Sim(pars, interventions = [bhpv_primary, bcytology, bcolpo, bablation], label='Modified Algorithm 2')

#Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
chpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=start_year, prob=cprimary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
# Send HPV+ women* for cytology
cto_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
ccytology   = hpv.routine_triage(eligibility=cto_cytology, prob=ctriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
#Send ASCUS and abnormal cytology results for colpo.
cto_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
ccolpo = hpv.routine_triage(eligibility=cto_colpo, prob=ctriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
chsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
cablation  = hpv.treat_num(eligibility=chsils, prob=cablate_prob, product='ablation', label='ablation') # Administer ablation

algo3 = hpv.Sim(pars, interventions = [chpv_primary, ccytology, ccolpo, cablation], label='Modified Algorithm 3')

###### NHS Algorithm For Comparison #######
screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt'])) & (sim.people.age >= 20) & (sim.people.age <= 65)

#Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
nhpv_primary  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=nhsprimary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
# Send HPV+ women* for cytology
nto_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive']
ncytology   = hpv.routine_triage(eligibility=nto_cytology, prob=nhstriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
#Send ASCUS and abnormal cytology results for colpo.
nto_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
ncolpo7 = hpv.routine_triage(eligibility=nto_colpo, prob=nhstriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
#After colpo, treat HSILs with Ablation
nhsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
nablation  = hpv.treat_num(eligibility=nhsils, prob=nhsablate_prob, product='ablation', label='ablation') # Administer ablation

algo = hpv.Sim(pars, interventions = [nhpv_primary, ncytology, ncolpo7, nablation], label='NHS Cervical Screening Pathway 2023')

# Create and run the simulation
orig_sim = hpv.Sim(pars, label='No Screening') 
'''
# Run the simulation
msim = hpv.parallel(algo1, algo, algo2, algo3)
msim.compare() #Gives quantitative results

#msim.to_excel('my-sim.xlsx')
#df = pd.read_excel('my-sim.xlsx')
#print(df)

# Plot the results
msim.plot('cancer_incidence')
'''

if __name__=="__main__":
    msim = hpv.parallel(orig_sim, algo, algo1, algo2, algo3) #makes and runs a multisim

    msim.compare()

    
    msim.plot()

    input("Stopping here")

    i=0
    for sim in msim.sims:
        print(sim)
        sim.to_excel(f'MultiSimSummary{i}_1.xlsx')
        i+=1

    

    #orig_sim.run()
    orig_sim.to_excel('algobasic.xlsx')
    df = pd.read_excel('algobasic.xlsx')
    print("Algobasic")
    print(df)



    #algo.run()
    algo.to_excel('algoNHS.xlsx')
    df = pd.read_excel('algoNHS.xlsx')
    print("AlgoNHS")
    print(df)

    #algo1.run()
    algo1.to_excel('algo1.xlsx')
    df = pd.read_excel('algo1.xlsx')
    print("Algo1")
    print(df)

    #algo2.run()
    algo2.to_excel('algo2.xlsx')
    df = pd.read_excel('algo2.xlsx')
    print("Algo2")
    print(df)

    #algo3.run()
    algo3.to_excel('algo3.xlsx')
    df = pd.read_excel('algo3.xlsx')
    print("Algo3")
    print(df)
