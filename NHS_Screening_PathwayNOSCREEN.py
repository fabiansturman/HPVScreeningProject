#No screening, only vaccination

#MODELLIGN ASSUMPTION: I THINK I AM NOT MODELLING HIV ALONGSIDE DUE TO NOT HAVING A HIV DATAFILE TO HAND! THAT SIMPLIFIES THE ALG SOMEWHAT

import hpvsim as hpv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

from copy import copy

#the seeds for which we want to generate results
seeds = list(range(10,50)) 

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
        run_name = f"results/NHSPathwayResultsNOSCREEN/NHSPathwayResults_seed{seeds[seed_number]}_basepars{base_parss_number}"

        #Trackers
        with_cancer = set([]) #this tracks all individuals who have cancer at some point 

        # Below, keeping track of {person id: (timepoint, scale factor for the agent at that timepoint) } for A) when agent with id pid first gets cancer, B) the timepoint at which they first get cancer
        time_with_first_cancer = dict()             #^by tracking scale factor, we can see how many people an agent represents at that point in time
        
        

        ## Modelling the ongoing NHS HPV Vaccination Strategy
        vacc_prob = 0.8

        vx1 = hpv.routine_vx(prob=vacc_prob, start_year=2008,end_year=2020, age_range=[12,13], product='quadrivalent') #vaccinate {vacc_prob} of girls age 12-13 every year, with Gardasil (quadrivalent) vaccine
        vx2 = hpv.routine_vx(prob=vacc_prob, start_year=2020, age_range=[12,13], sex=[0,1], product='nonavalent') #start to vaccinate boys as well
        vx3 = TODO: from 2023 onwards, only 1 dose!!
            #TODO: refine this with info on Gardasil and Gardasil 9, in terms of efficacy against strains and waning (if notable)

        
        def update_with_cancer(sim):
            #Keep track of my set which tracks who is with cancer
            for pid in range(sim.n):
                if sim.people.alive[pid] and sim.people.cancerous.any(axis=0)[pid]:
                    if pid not in with_cancer:
                        with_cancer.add(pid)
                        time_with_first_cancer[pid] = (sim.t, sim.people.scale[pid]) #we note there a cancerous agent represents 10x as many people as a normal agent
        update_with_cancer.label="update_with_cancer"


        
        algo = hpv.Sim(pars, interventions = [vx1,vx2, 
                                            update_with_cancer,
                                            ], #interventions can just be functions of a sim
                                                label='NHS Cervical Screening Pathway 2025')
        
        algo.initialize()#need to initialize so that people are created and then we can access them!



       

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
                        'time_with_first_cancer': time_with_first_cancer,
                        }, file)


        algo.to_excel(f'{run_name}.xlsx')
        
        