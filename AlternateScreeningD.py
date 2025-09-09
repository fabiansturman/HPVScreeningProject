#Implementing Alternate Screening Pathway (ASP) D: Normal NHS until 2025. Starting 2025, vaccinated women have a 10year/10year callback (rather than 5 year/3year)

#MODELLIGN ASSUMPTION: I THINK I AM NOT MODELLING HIV ALONGSIDE DUE TO NOT HAVING A HIV DATAFILE TO HAND! THAT SIMPLIFIES THE ALG SOMEWHAT

#TODO: as noted in study journal, I want to get some more intricate measures when I redo results properly 
#   -> track intervention usage + outcomes over time, to compare between algs (esp in terms of screening useage/efficiency, etc) :a function used as a custom intervention can monitor intervention.outcomes for all desired interventions, and the sum over all of these is total intervention usage in that timestep for example
        #^ especially needed to properly see what is going on with ASP C, to see how many fewer later screenings we do when we restrict it only to the unvaccinated
#   -> split up all results by vaccination status, so that we can check for fairness between the two (esp comparatively)

import hpvsim as hpv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

from copy import copy

#the following is the correct timepoint at which to switch screening methods, as long as we start simulating at 1970 and have dt=0.25
tp_switch = int((2025-1970)/0.25)

#the seeds for which we want to generate results
seeds = list(range(50))

#the different simulation parameters over which (for each seed) we want to generate our results
base_parss = [dict(
            n_agents      = 20e3,       
            start=1970, end=2080,
            #verbose       = 0,   
            rand_seed     = 1,   
            genotypes     = [16, 18, 'hr'],
            burnin=30,
            location='united kingdom',
            ),
            ]

#cartesian product of seeds and base parameters for our sim
parss = []
for base_pars in base_parss:
    for seed in seeds:
        pars = copy(base_pars)
        pars['rand_seed']=seed
        parss.append(pars)


if __name__=="__main__":
    run_number=-1
    for pars in parss:
        #Unique filename under which to save these results
        run_number+=1
        base_parss_number = run_number // len(seeds)
        seed_number = run_number % len(seeds)
        run_name = f"results/AlternateScreeningD/DResults_seed{seeds[seed_number]}_basepars{base_parss_number}"


        #Trackers
        colposcopy_eligible_with_cancer = set([]) #this will store all the people id's of those with cancer who (are eligible to) undergo cytology.
        #TODO: justify (to myself, then written) that this is a fair proxy for #cancers caught (at least for comparative purposes)
        with_cancer = set([]) #this tracks all individuals who have cancer at some point 

        # Below, keeping track of {person id: (timepoint, scale factor for the agent at that timepoint) } for A) when agent with id pid first gets cancer, B) the timepoint at which they first get cancer
        time_with_first_cancer = dict()             #^by tracking scale factor, we can see how many people an agent represents at that point in time
        time_with_first_cancer_colposcopy = dict()

        

        ## Modelling the ongoing NHS HPV Vaccination Strategy
        vacc_prob = 0.8

        vx1 = hpv.routine_vx(prob=vacc_prob, start_year=2008,end_year=2020, age_range=[12,13], product='quadrivalent') #vaccinate {vacc_prob} of girls age 12-13 every year, with Gardasil (quadrivalent) vaccine
        vx2 = hpv.routine_vx(prob=vacc_prob, start_year=2020, age_range=[12,13], sex=[0,1], product='nonavalent') #start to vaccinate boys as well
            #TODO: refine this with info on Gardasil and Gardasil 9, in terms of efficacy against strains and waning (if notable)


        ## Modelling the current NHS Cervical Cancer Screening Strategy
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

            
            for i in range(sim.n): #sim.n==len(sim.people)
                if i in sim.people.last_hpv_result.keys():
                    #vaccinated after 2025 get longer callback intervals 
                    callback_intervals = (3,5) if (sim.t < tp_switch or not sim.people.vaccinated[i]) else (10,10)

                    last_hpv_result = sim.people.last_hpv_result[i]
                    age = sim.people.age[i]
                    date_screened = sim.people.date_screened[i]

                    if last_hpv_result==1 and 24<=age<50 and sim.t>date_screened+ callback_intervals[0]/sim['dt']:
                        for_screening[i] = True
                    if last_hpv_result==-1 and 24<=age and sim.t > date_screened+  callback_intervals[1]/sim['dt']:
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
        

        
        algo = hpv.Sim(pars, interventions = [vx1,vx2, 
                                            routine_screening, first_cytology, second_consecutive_screening, 
                                                second_cytology, third_consecutive_screening, third_cytology, 
                                                colposcopy, ablation, 
                                            update_last_hpv_result, update_needs_consec_screening_2,  update_needs_consec_screening_3,
                                            update_with_cancer,
                                            ], #interventions can just be functions of a sim
                                                label='NHS Cervical Screening Pathway 2025')
        
        algo.initialize()#need to initialize so that people are created and then we can access them!
        add_last_hpv_result(algo) #so we can keep track of outcomes from hpv screening tests persistently across years
        add_needs_consec_screening_2(algo)
        add_needs_consec_screening_3(algo)



       

        algo.run()

        """
        print(algo.t)
        cutoff_timepoint = float(input("Enter cutoff timepoint:"))


        print(f"Agents with cancer by end: {len(with_cancer)}, of which {len(colposcopy_eligible_with_cancer)} had a colposcopy while having cancer.")
        cum_cancer_people = 0
        cum_colposcopy_cancer_people = 0

        for pid in with_cancer:
            tp,scale = time_with_first_cancer[pid]
            if tp<=cutoff_timepoint:
                cum_cancer_people += scale 
        for pid in colposcopy_eligible_with_cancer:
            tp, scale = time_with_first_cancer_colposcopy[pid]
            if tp<=cutoff_timepoint:
                cum_colposcopy_cancer_people += scale

        print(f"People with cancer by end: {cum_cancer_people}, of which {cum_colposcopy_cancer_people} had a colposcopy while having cancer.")
        """
        with open(f'{run_name}.pickle', "wb") as file:
            pickle.dump({'final_timepoint':algo.t,
                        'with_cancer':with_cancer,
                        'colposcopy_eligible_with_cancer':colposcopy_eligible_with_cancer,
                        'time_with_first_cancer': time_with_first_cancer,
                        'time_with_first_cancer_colposcopy':time_with_first_cancer_colposcopy}, file=file)


        algo.to_excel(f'{run_name}.xlsx')
        
        