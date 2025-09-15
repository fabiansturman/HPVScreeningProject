"""
This file contains an model of the current NHS England Cervical Screening Algorithm.
    ^https://www.gov.uk/government/publications/cervical-screening-care-pathway/cervical-screening-care-pathway (version as of 1st July 2025)

File usage:
    -> Import this file to easily integrate screening algorithm model into a HPVsim simulation
    -> Run this file (as __main__) to run a HPVsim simulation with this screening and debug changes to the algorithm

Notable modelling assumptions:
    -> HIV is not being modelled in this study, justified by low HIV prevalence in the UK, particularly among those with cervixes. Therefore, we do not include branches of the pathway which relate to HIV, assuming all individuals are HIV -ve.

Below are sources to justify modelling decisions.

Treatment effectiveness:
    -> HPV DNA testing has been modelled with 96% sensitivity to HPV (https://doi.org/10.1002/ijc.31094)
    -> Cytology has been modelled with 70% sensitivity to low grade CIN, 85% for high grade CIN, and 92% for asymptomatic cancer (https://doi.org/10.1002/ijc.31094)
    -> POST CYTOLOGY DIAGNOSTICS/TREATMENT has been modelled by https://doi.org/10.1002/ijc.31094 as:
            - if positive cytology with state='high-grade CIN', 90% assumed treated successfully (allows for some treatment failure and some women not to return for treatment)
            - if positive cytology with state='low-grade CIN', 5% assmed treated successfully (i.e. preventing possible future progression to high grade CIN without further infection)
            - of the women with three consecutive HPV +ve, cytology -ve results, they have a colposcopy (like we model), and 
                    -90% of high-grade CIN treated successfully ONLY


    -> Excition:
            - the HPVsim paper itself says they have modelled them as  95% effective at removing lesions (which, from some googling, should incldue early-stage cervical cancer), and 80% effective at clearing active HPV infection. I have made it 80% effective at traeting preCIN and 95% at CIN, and 0% for cancerous for now because HPVsim doesnt differentiate between cancer stages very well from what I can tell. FOr some reason, this isnt what HPVsim added into teh code even though they say it in the paper, so it doesnt really need a reference 

Compliance:
    -> Cervical screening coverage 2013-2023 ranges from 68-74%: (73.93386232,74.16190589,73.45176973,72.6770978,72.04060386,71.41738301,71.89922063,72.17636968,70.21941735,69.9326474,68.70463328)% (https://digital.nhs.uk/data-and-information/publications/statistical/cervical-screening-annual/england-2022-2023, published Nov 2023. See the screenshot from Excel that I have saved for details on how I have calculated these figures from the NHS data)
            - from the same data source, we see that between 2017 and 2023, 15-20% of screens each year are not prompted by the programme (i.e. booked in for a screening otherwise). 
                ^should i be modelling this?
            

"""


#--- Imports ---#
import hpvsim as hpv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

from copy import copy

## Modelling the current NHS Cervical Cancer Screening Strategy


#Trackers - i think they should initialise like the sim.people properties but directly as properties of the sim
colposcopy_eligible_with_cancer = set([]) #this will store all the people id's of those with cancer who (are eligible to) undergo cytology.
#TODO: justify (to myself, then written) that this is a fair proxy for #cancers caught (at least for comparative purposes)
with_cancer = set([]) #this tracks all individuals who have cancer at some point 

# Below, keeping track of {person id: (timepoint, scale factor for the agent at that timepoint) } for A) when agent with id pid first gets cancer, B) the timepoint at which they first get cancer
time_with_first_cancer = dict()             #^by tracking scale factor, we can see how many people an agent represents at that point in time
time_with_first_cancer_colposcopy = dict()


start_year = 1980#2015
primary_screen_prob = 1#0.8#1     #TODO: this probability should encompass A)all women who have opted out of screening (reminders); B)the probability of a woman who has not opted out of screening (reminders) actually taking up a given appointment oppertunity; C)appointment availability
secondary_screen_prob = 1#0.8#1
third_screen_prob = 1#0.8#1
triage_screen_prob = 1#0.8#1      #TODO: as above
colpo_prob = 1#0.8#1
ablate_prob = 1#0.8#1             #TODO: as above

# Adding attributes to our sim's people so that we can track the outcomes of their interventions
def add_last_hpv_result(sim):
    #Add 'last_hpv_result' to sim.people as a dictionary attribute 1=positive, -1=negative.
    sim.people.last_hpv_result = {} #dictionary {agent_id: most_recent_HPV_result} to handle dynamic populations

def add_needs_consec_screening_2(sim):
    #Adds a property to sim.people that contains a timepoint from which we request the patient undergoes a 2nd consecutive screening. If None, no 2nd consec screening needed.
    sim.people.needs_consec_screening_2 = {}

def add_needs_consec_screening_3(sim):
    #Adds a property to sim.people that contains a timepoint from which we request the patient undergoes a 3rd consecutive screening. If None, no 3rd consec screening needed.
    sim.people.needs_consec_screening_3 = {}

#Custom interventions to update attributes that keep track of the other interventions' results, persistently across years
def update_last_hpv_result(sim, test_name="routine_screening"):
    #For everyone who has undergone a screening of some sort, update this result
    hpv_test_names = ["routine_screening", "second_consecutive_screening", "third_consecutive_screening"]
        #the above list is in order of administration of the tests within a given timepoint; noone should undergo several of the same tests at any point, but just in case we want to use the most recent result!
    for name in hpv_test_names:
        hpv_test = sim.get_intervention(name)
        for positive in hpv_test.outcomes.get('positive'):
            sim.people.last_hpv_result[positive] = +1
        for negative in hpv_test.outcomes.get('negative'):
            sim.people.last_hpv_result[negative] = -1
update_last_hpv_result.label="update_last_hpv_result" #by adding an intervention as a function, when searching 'get_intervention', interventions without a label attribute cause an error

def update_needs_consec_screening_2(sim, name="first_cytology"):
    #If a woman has a -ve result from her first cytology, she is eligible for a second consecutive screening in 12 months
    for i in sim.get_intervention(name).outcomes['normal']: #iterate over all negative outcomes from first cytology
        sim.people.needs_consec_screening_2[i] = sim.t+1/sim["dt"]
    
    
    #All women who have undergone their second conseucutive screening have this updated in their record (no more need for 2nd consecutive screening)
    for i in sim.get_intervention("second_consecutive_screening").outcomes['positive']:
        sim.people.needs_consec_screening_2[i] = None
    for i in sim.get_intervention("second_consecutive_screening").outcomes['negative']:
        sim.people.needs_consec_screening_2[i] = None
    for i in sim.get_intervention("second_consecutive_screening").outcomes['inadequate']:
        sim.people.needs_consec_screening_2[i] = None
        
update_needs_consec_screening_2.label="update_needs_consec_screening_2"

def update_needs_consec_screening_3(sim, name="second_cytology"):
    #If a woman has a -ve result from her second cytology, she is eligible for a third consecutive screening in 12 months
    for i in sim.get_intervention(name).outcomes['normal']: #iterate over all negative outcomes from first cytology
        sim.people.needs_consec_screening_3[i] = sim.t+1/sim["dt"]
        


    #All women who have undergone their second conseucutive screening have this updated in their record (no more need for 2nd consecutive screening)
    for i in sim.get_intervention("third_consecutive_screening").outcomes['positive']:
        sim.people.needs_consec_screening_3[i] = None
    for i in sim.get_intervention("third_consecutive_screening").outcomes['negative']:
        sim.people.needs_consec_screening_3[i] = None
    for i in sim.get_intervention("third_consecutive_screening").outcomes['inadequate']:
        sim.people.needs_consec_screening_3[i] = None
update_needs_consec_screening_3.label="update_needs_consec_screening_3"

def update_with_cancer(sim):
    #Keep track of my set which tracks who is with cancer
    for pid in range(sim.n):
        if sim.people.alive[pid] and sim.people.cancerous.any(axis=0)[pid]:
            if pid not in with_cancer:
                with_cancer.add(pid)
                time_with_first_cancer[pid] = (sim.t, sim.people.scale[pid]) #we note there a cancerous agent represents 10x as many people as a normal agent
update_with_cancer.label="update_with_cancer"



# Setting up routine screening, according to last HPV screening outcome (if any), time since last screening, and age
def routine_screen_eligible(sim):
    #A woman is eligible for a routine screening {3 IF last hpv test +ve and under 50 years old ELSE 5} years after their most recent screening
    for_screening = np.isnan(sim.people.date_screened) & (sim.people.age>=24) & (sim.people.age<=60) #eligible for screening if they have not yet been screened

        #^I think I can safely conclude that indices are retained for PEOPLE, and when new people are born they are given fresh indices. In which case, I can absoltuely use a dictionary to store treatment outcomes for people, indexed by thier sim.people index!
        #TODO: confirm the above
    
    for i in range(sim.n): #sim.n==len(sim.people)
        if i in sim.people.last_hpv_result.keys():
            last_hpv_result = sim.people.last_hpv_result[i]
            age = sim.people.age[i]
            date_screened = sim.people.date_screened[i]

            if last_hpv_result==1 and 24<=age<50 and sim.t>date_screened+3/sim['dt']:
                for_screening[i] = True
            if last_hpv_result==-1 and 24<=age and sim.t > date_screened+5/sim['dt']:
                for_screening[i] = True
            
    return for_screening
    #TODO: verify that the above function definitely works as I want


routine_screening  = hpv.routine_screening(eligibility=routine_screen_eligible, 
                                    start_year=start_year, 
                                    prob=primary_screen_prob, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='routine_screening') 



# First (Consecutive) Cyotology for women with positive routine screening
    #[We will refer to actions as 'consecutive' if they can all be traced back to the same initial routine screening]
to_first_cytology = lambda sim: sim.get_intervention('routine_screening').outcomes['positive']
first_cytology = hpv.routine_triage(eligibility=to_first_cytology,
                                    prob = triage_screen_prob,
                                    product='lbc',
                                    annual_prob=False,
                                    label="first_cytology")

#If first cytology negative, send to a screening in a year (i.e. a second consecutive screening)
def second_screening_eligible(sim):
    for_screening = np.array([False,]*sim.n)
    for i in sim.people.needs_consec_screening_2.keys():
        if sim.people.needs_consec_screening_2[i] is not None and sim.people.needs_consec_screening_2[i] >= sim.t:
            for_screening[i]=True
    return for_screening

second_consecutive_screening  = hpv.routine_screening(eligibility=second_screening_eligible, 
                                    start_year=start_year, 
                                    prob=secondary_screen_prob, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='second_consecutive_screening') 

#-ve from second screening are returned to standard screening - their date_screened is updated from the second screening, so nothing needed for them
#+ve from second screening do second cytology
to_second_cytology = lambda sim: sim.get_intervention('second_consecutive_screening').outcomes['positive']
second_cytology = hpv.routine_triage(eligibility=to_second_cytology,
                                    prob = triage_screen_prob,
                                    product='lbc',
                                    annual_prob=False,
                                    label="second_cytology")

#-ve from second cytology do a third consecutive HPV screening in 12 months
def third_screening_eligible(sim):
    for_screening = np.array([False,]*sim.n)
    for i in sim.people.needs_consec_screening_3.keys():
        if sim.people.needs_consec_screening_3[i] is not None and sim.people.needs_consec_screening_3[i] >= sim.t:
            for_screening[i]=True
    return for_screening

third_consecutive_screening  = hpv.routine_screening(eligibility=third_screening_eligible, 
                                    start_year=start_year, 
                                    prob=third_screen_prob, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='third_consecutive_screening') 

#-ve from third screening are returned to standard screening. As their latest screening time is updated in {date_screened}, no need to do anything here
#+ve from third screening perform a third cyotology followed by a mandatory colposcopy
to_third_cytology = lambda sim: sim.get_intervention('third_consecutive_screening').outcomes['positive']
third_cytology = hpv.routine_triage(eligibility=to_third_cytology,
                                    prob = triage_screen_prob,
                                    product='lbc',
                                    annual_prob=False,
                                    label="third_cytology")

#Eligible for a colposcopy if any of the following hold: A)tested not normal for any cytology; B) went all the way to a third screening/cytology
colposcopy_eligible = lambda sim:list(set(sim.get_intervention('first_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('first_cytology').outcomes['ascus'].tolist() +
                                        sim.get_intervention('second_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('second_cytology').outcomes['ascus'].tolist() +
                                        sim.get_intervention('third_cytology').outcomes['normal'].tolist() + sim.get_intervention('third_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('third_cytology').outcomes['ascus'].tolist()))

def colposcopy_eligible_and_track_cancers(sim):
    #Does the same thing as colposcopy_eligible, but also tracks the ids of all cancerous patients that go through colposcopy process
    eligible=colposcopy_eligible(sim)
    for pid in range(sim.n):
        if sim.people.alive[pid] and (pid in eligible)  and sim.people.cancerous.any(axis=0)[pid]:
            if pid not in colposcopy_eligible_with_cancer:
                colposcopy_eligible_with_cancer.add(pid)
                time_with_first_cancer_colposcopy[pid] = (sim.t, sim.people.scale[pid])
    return eligible

colposcopy =  hpv.routine_triage(eligibility=colposcopy_eligible_and_track_cancers,#colposcopy_eligible, 
                                prob=colpo_prob, 
                                product='colposcopy', 
                                annual_prob=False, 
                                label='colposcopy')

#TODO: improve the below - have a more nuanced response than just ablating every positive colposcopy, especially w.r.t colposcopies able to diagnose cancer! Perhaps design a treatment for that, just based on average succeses of cancer treatments!
hsil_by_colpo = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] #i.e. covers all women diagnosed with HSIL by a colposcopy
ablation = hpv.treat_num(eligibility=hsil_by_colpo, 
                        prob=ablate_prob, 
                        product='ablation', 
                        label='ablation') #TODO: can model a max capacity with treat_num too, but as it stands we pretend infinite capacity





if __name__=="__main__":
    #HPVsim simulation parameters
    pars = dict(
            n_agents      = 20e3,       
            start=1970, end=2200,#2080,
            #verbose       = 0,   
            rand_seed     = 1, 
            genotypes     = [16,18, 'hi5'],
            burnin=30,
            location='united kingdom',
            )


    #Define simulation without interventions
    sim_basic = hpv.Sim(pars,
                        label='Sim, no interventions')

    #Define simulation with interventions
    sim_interventions =  hpv.Sim(pars, interventions = [ 
                                        routine_screening, first_cytology, second_consecutive_screening, 
                                            second_cytology, third_consecutive_screening, third_cytology, 
                                            colposcopy, ablation, 
                                        update_last_hpv_result, update_needs_consec_screening_2,  update_needs_consec_screening_3,
                                        update_with_cancer,
                                        ]
                                 , #interventions can just be functions of a sim
                                            label='NHS Cervical Screening Pathway 2025')
    sim_interventions.initialize() #need to initialize so that people are created and then we can access them! below, we define new properties of people that are persistent over timesteps 
    add_last_hpv_result(sim_interventions)
    add_needs_consec_screening_2(sim_interventions)
    add_needs_consec_screening_3(sim_interventions)

    #Run and plot sims
    reduced_sims = []
    for sim in [sim_basic, sim_interventions]:
        msim = hpv.MultiSim(sim, n_runs=10)  
        msim.run()         
        reduced_sims.append(msim.reduce(use_mean=True, output=True))  


    msim = hpv.MultiSim(reduced_sims)
    msim.plot()