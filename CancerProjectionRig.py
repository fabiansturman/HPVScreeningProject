import hpvsim as hpv
import numpy as np
import pylab as plt
import pandas as pd

#TODO: confirm this function is only used for those who are having ROUTINE recall (that is, nothing abnormal found last time)
#Function to determine who is eligible for screening according to their last time being screened and their recall period - which here depends on vaccination status
def screen_eligible_with_recall(sim, unvaccinated_recall_years=3, vaccinated_recall_years=6):#5.5):
    # Determine the recall period for each individual based on vaccination status
    recall_periods = np.where(sim.people.vaccinated, unvaccinated_recall_years, vaccinated_recall_years) #element-wise operation to handle arrays correctly
    
    # Calculate the time until next screening for each individual
    time_until_rescreen = recall_periods / sim['dt']
    
    # Determine if it is time for re-screen based on individual recall periods
    is_time_for_rescreen = np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + time_until_rescreen))
    
    # Determine eligibility based on time for re-screen, age, and possibly immunity level
    is_eligible = is_time_for_rescreen & (sim.people.age >= 24) & (sim.people.age <= 65)
    
    return is_eligible

############################

## Define parameters for our underlying HPVsim simulation
pars = dict(
    n_agents      = 20e3,
    start         = 1980,            
    end           = 2080,
    verbose       = 0,               
    rand_seed     = 2,               # TODO: Set seeds smarter; iterating thorugh them or smth
    burnin        = 30,              
    genotypes     = [16, 18, 'hr'],  # Include all high risk HPV genotypes
)
#TODO: update parameters with results from calibration(s)!



## Define parameters for the interventions
#We share these parameters between all interventions for fair comparison between them
intervention_start_year = 2021
primary_screen_prob = 0.6 #TODO: justify the choice of this prob

#NHS' Current Algorithm: unique parameters
nhstriage_screen_prob = 0.9
nhsablate_prob = 0.9

#Modified Algorithm 1: unique parameters
atriage_screen_prob = 0.74
aablate_prob = 0.9

#Modified Algorithm 2: unique parameters
btriage_screen_prob = 0.8
bablate_prob = 0.9

#Modified Algorithm 3: unique parameters
    #NOTE: for some reason, this was given a primary screen prob of 0.55 by Sophie, the only of the algorithms to have this
ctriage_screen_prob = 0.85
cablate_prob = 0.9




#TODO: does it make sense for the above parameters not to be shared? Surely for a fair comparison they should be shared??

# Define a series of interventions to screen, triage, assign treatment, and administer treatment



"""
Modified Algorithm 1:
- Primary Screening Test: HPV DNA Testing  (self sampled, or collected by clinicians)
- IF hrHPV POSITIVE:
  - Refer to Cervical Cyotology Triage:
  - IF CYOTOLOGY ABNORMAL:
    - Refer to colposcopy; further management based on coloscopy/histopathology diagnosis
  - IF CYOTOLOGY NORMAL:
    - Repeat screening (HPV DNA test) in 12 months
- IF hrHPV NEGATIVE:
  - Routine rescreening in {6 if vaccinated, else 3} years
"""
hpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=intervention_start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing

to_cytology = lambda sim: sim.get_intervention('hpv_primary').outcomes['positive'] #defines eligibility function, where we take a sim and return all people who have had a hpv_primary intervention come back as positive for hrhpv
cytology   = hpv.routine_triage(eligibility=to_cytology, prob=atriage_screen_prob, product='lbc', annual_prob=False, label='cytology')
    #{annual_prob} refers to whether probabilities are annual or per-timestep. By setting False, probabilities are per-timestep

to_colpo = lambda sim:list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist())) #Send ASCUS (Atypical Squamous Cells of Undetermined Significance) and abnormal cytology results for colposcopy
colpo = hpv.routine_triage(eligibility=to_colpo, prob=atriage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') 

hsils = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] #After colpo, treat HSILs (high-grade squamous intraepithelial lesions) with Ablation
ablation  = hpv.treat_num(eligibility=hsils, prob=aablate_prob, product='ablation', label='ablation') # Administer ablation

algo1 = hpv.Sim(pars, interventions = [hpv_primary, cytology, colpo, ablation], label='Modified Algorithm 1')

#TODO: this has not accounted for when we have a NORMAL coloscopy, has it? Check the HPVsim tutorial's examples of implementing WHO algoriuthms to see if they also do this and omit this part, which implies maybe we don't need to put that in, but otherwise i think we do!
#TODO: this doesnt account for the complexities of cancer treatments after a positive colposcopy - it just treats them all with ablation. See what else we can do, if anything. THis can of course just be a modelling simplification though - especially if we only model #cancers, rather than #deaths from cancer, this sort of is all we need ; I suppose careful modelling of stuff after a cancer diagosis isnt much
# ^ but IS THIS(?) a cancer diagnosis. Perhaps HPVsim doesnt log this as a detected cancer, so just really double check the examples of implemented algorithms 




"""
Modified Algorithm 2:
- Primary Screening Test: HPV DNA Testing  (self sampled, or collected by clinicians)
- IF hrHPV POSITIVE:
  - Refer to Cervical Cyotology Triage:
  - IF CYOTOLOGY ABNORMAL:
    - Refer to colposcopy; further management based on coloscopy/histopathology diagnosis
  - IF CYOTOLOGY NORMAL:
    - Repeat screening (HPV DNA test) in {12 if unvaccinated, else 24} months
- IF hrHPV NEGATIVE:
  - Routine rescreening in 3 years
"""
bhpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=intervention_start_year, prob=primary_screen_prob, product='hpv', label='bhpv_primary') #HPV Testing

bto_cytology = lambda sim: sim.get_intervention('bhpv_primary').outcomes['positive']
bcytology   = hpv.routine_triage(eligibility=bto_cytology, prob=btriage_screen_prob, product='lbc', annual_prob=False, label='bcytology')

bto_colpo = lambda sim:list(set(sim.get_intervention('bcytology').outcomes['abnormal'].tolist() + sim.get_intervention('bcytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
bcolpo = hpv.routine_triage(eligibility=bto_colpo, prob=btriage_screen_prob, product='colposcopy', annual_prob=False, label='bcolposcopy') #Send people to colposcopy

bhsils = lambda sim: sim.get_intervention('bcolposcopy').outcomes['hsil'] # Define who's eligible for ablation
bablation  = hpv.treat_num(eligibility=bhsils, prob=bablate_prob, product='ablation', label='bablation') # Administer ablation

algo2 = hpv.Sim(pars, interventions = [bhpv_primary, bcytology, bcolpo, bablation], label='Modified Algorithm 2')

#TODO: I dont think this modified algorithm 2 is any different to modified algorithm 1? was it just the different probabilities that made the differences????


"""
Modified Algorithm 3:
- Primary Screening Test: HPV DNA Testing  (self sampled, or collected by clinicians)
- IF hrHPV POSITIVE:
  - Refer to Cervical Cyotology Triage:
  - IF CYOTOLOGY ABNORMAL:
    - Refer to colposcopy; further management based on coloscopy/histopathology diagnosis
  - IF CYOTOLOGY NORMAL:
    - Repeat screening (HPV DNA test) in {12 if unvaccinated, else 24} months
- IF hrHPV NEGATIVE:
  - Routine rescreening in {6 if unvaccinated, else 3} years
"""
chpv_primary  = hpv.routine_screening(eligibility=screen_eligible_with_recall, start_year=intervention_start_year, prob=primary_screen_prob, product='hpv', label='chpv_primary') #HPV Testing

cto_cytology = lambda sim: sim.get_intervention('chpv_primary').outcomes['positive']
ccytology   = hpv.routine_triage(eligibility=cto_cytology, prob=ctriage_screen_prob, product='lbc', annual_prob=False, label='ccytology')

cto_colpo = lambda sim:list(set(sim.get_intervention('ccytology').outcomes['abnormal'].tolist() + sim.get_intervention('ccytology').outcomes['ascus'].tolist())) #Define who's eligible for colpo
ccolpo = hpv.routine_triage(eligibility=cto_colpo, prob=ctriage_screen_prob, product='colposcopy', annual_prob=False, label='ccolposcopy') #Send people to colposcopy

chsils = lambda sim: sim.get_intervention('ccolposcopy').outcomes['hsil'] # Define who's eligible for ablation
cablation  = hpv.treat_num(eligibility=chsils, prob=cablate_prob, product='ablation', label='cablation') # Administer ablation

algo3 = hpv.Sim(pars, interventions = [chpv_primary, ccytology, ccolpo, cablation], label='Modified Algorithm 3')

###### NHS Algorithm For Comparison #######
screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt'])) & (sim.people.age >= 20) & (sim.people.age <= 65)

#Algorithm. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
nhpv_primary  = hpv.routine_screening(eligibility=screen_eligible, start_year=intervention_start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary') #HPV Testing
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
