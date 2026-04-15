"""
This file offers a function which, if added to the end of an intervention list, allows for logging of each agent over time
It also offers some functionality to process/interpret logs - but the use of logs to generate the core results of this project are instead in resultGeneration.py
"""
#Imports
from tqdm import tqdm
import numpy as np
import pickle
from copy import deepcopy

from math import ceil

import hpvsim as hpv

import matplotlib.pyplot as plt

from InterventionAlgorithms.simulationLogging import *


if __name__=="__main__":
    from basePars import base_pars
    from InterventionAlgorithms import NHS_2025_lambdamu, NHS_Vacc , GlobalScreeningParameters, NHS_2025_lambdamuForVaccCohortsOnly

    start = 1980
    end= 2100#2055 

    logger, log = makeLoggingIntervention()

    base_pars = deepcopy(base_pars)

    base_pars['rand_seed'] = 1
    base_pars['start'] = start
    base_pars['end'] = end
    base_pars['n_agents']= 10_000#200_000
    base_pars['interventions'] = NHS_2025_lambdamu.get_interventions(l=1, m=1) + [logger]
    #base_pars['interventions'] = NHS_2025_lambdamu.get_interventions(l=1, m=1) + NHS_Vacc.vaccinations + [logger]
 #   base_pars['interventions'] = NHS_Vacc.vaccinations + [logger]
   # base_pars['interventions'] = NHS_2025_lambdamuForVaccCohortsOnly.get_interventions(l=1, m=1) + NHS_Vacc.vaccinations + [logger]

    sim = hpv.Sim(base_pars) 
    sim.run()
    #sim.plot()

    """
    ###-RUN THIS CODE TO INVESTIGATE VACCINATION ROLLOUT-###
    doses = log[len(log.keys())-1]['vx_doses']
    counts = {}
    for dose in doses:
        if dose not in counts.keys():
            counts[dose] = 0
        counts[dose] +=1
    print(f"(x:num people with x doses): {counts}")
    

    v1 = get_vacinnated_proportion_timeseries(log, 1, progress_bar=True)
    v2 = get_vacinnated_proportion_timeseries(log, 2, progress_bar=True)
    v3 = get_vacinnated_proportion_timeseries(log, 3, progress_bar=True)

    vacc_proportion_timeseries_1dose = [v1[i] for i in range(len(log.keys()))]
    vacc_proportion_timeseries_2doses = [v2[i] for i in range(len(log.keys()))]
    vacc_proportion_timeseries_3doses = [v3[i] for i in range(len(log.keys()))]

    xs = np.arange(start, end+1, base_pars['dt'])

    plt.plot(xs, vacc_proportion_timeseries_1dose, label="1 dose") 
    plt.plot(xs, vacc_proportion_timeseries_2doses, label = "2 doses")
    plt.plot(xs, vacc_proportion_timeseries_3doses, label = "3 doses")
    plt.hlines(GlobalScreeningParameters.projected_teen_vaccination_uptake, start, end, color='green', linestyle='dotted')
    plt.ylim((0,1))
    plt.legend()
    plt.show()

    cumulative_doses = get_cumulative_doses_timeseries(log)
    cumulative_doses_for_plotting = [cumulative_doses[i] for i in range(len(log.keys()))]

    plt.plot(xs, cumulative_doses_for_plotting, label='cumulative doses over time')
    plt.show()

    #NOTE: we see the behaviour we expect - it may take some additional years (e.g. running to year 2200 just for good measure) to get to equilibruim in the proportion of the vaccinated adult population, but once we are there it is all correct.
    """



    #"""
    ###-RUN THIS CODE FOR GENERATING SANKEY DIAGRAMS just take out the quotations at top and bottom to get a basic Sankey diagram, and add filters to focus on alternative subsets of the population)-###
    
    #ALT FOCUS 0: PRESENTING PATH FOR NON-CANCEROUS AGENTS#
   # log = filter_intervention_by_cancerous(log, False, 'routine_screening_under50', relative_search_zone=[-1], progress_bar=True)
   # log = filter_intervention_by_cancerous(log, False, 'routine_screening_50andover', relative_search_zone=[-1], progress_bar=True)
    #ALT FOCUS 0: PRESENTING PATH FOR NON-CANCEROUS AGENTS#
        #NOTE: running ALT FOCUS 0, particularly comparing a run with our modelled colposcopy specificity compared to setting specificity to 1, we see that the manu failed ablations are just from there being so many more healthy agents getting screened that most results are a healthy agents getting an ablation and are victims of a low test specificity (not the ablation itself giving a fail result)

    #ALT FOCUS 1: PRESENTING PATH FOR CANCEROUS AGENTS#
  #  log = filter_intervention_by_cancerous(log, True, 'routine_screening_under50', relative_search_zone=[-1])
  #  log = filter_intervention_by_cancerous(log, True, 'routine_screening_50andover', relative_search_zone=[-1])
    #ALT FOCUS 1: PRESENTING PATH FOR CANCEROUS AGENTS#

    #ALT FOCUS 1: PRESENTING PATH FOR CANCEROUS AGENTS#
  #  log = filter_intervention_by_cancerous(log,True, 'routine_screening_under50', relative_search_zone=[-1])
  #  log = filter_intervention_by_cancerous(log,True, 'routine_screening_50andover', relative_search_zone=[-1])
    #ALT FOCUS 1: PRESENTING PATH FOR CANCEROUS AGENTS#

    #ALT FOCUS 2: PRESENTING PATH FOR CANCEROUS, OR SOON TO BE CANCEROUS, AGENTS#
 #   log = filter_intervention_by_cancerous(log,True, 'routine_screening_under50', relative_search_zone=[-1, 0, 1, 2,3,4], progress_bar=True)
 #   log = filter_intervention_by_cancerous(log,True, 'routine_screening_50andover', relative_search_zone=[-1, 0, 1, 2,3,4], progress_bar=True)
    #ALT FOCUS 2: PRESENTING PATH FOR CANCEROUS, OR SOON TO BE CANCEROUS, AGENTS#

    #ALT FOCUS 3: PRESENTING PATH FOR CIN AGENTS#
 #   log = filter_intervention_by_cin(log,True,  'routine_screening_under50', relative_search_zone=[-1])
 #   log = filter_intervention_by_cin(log, True, 'routine_screening_50andover', relative_search_zone=[-1])
    #ALT FOCUS 3: PRESENTING PATH FOR CIN AGENTS# - its reassuring to see that ablation does succeed consistently with this subgroup! it just must not succeed that great overall becasue both cancers and non-cins end up getting passed to it quite a bit?

    #ALT FOCUS 4: PRESENTING PATH FOR CIN, OR SOON TO BE CIN, AGENTS#
 #   log = filter_intervention_by_cin(log, True, 'routine_screening_under50', relative_search_zone=[-1, 0, 1, 2,3,4], progress_bar=True)
 #   log = filter_intervention_by_cin(log, True, 'routine_screening_50andover', relative_search_zone=[-1, 0, 1, 2,3,4], progress_bar=True)
    #ALT FOCUS 4: PRESENTING PATH FOR CIN, OR SOON TO BE CIN, AGENTS#

    #ALT FOCUS 5: PRESENTING PATH JUST FOR UNVACCINATED AGENTS#
    #log = filter_intervention_by_vacc(log, keep_vacc=False)
    #ALT FOCUS 5: PRESENTING PATH JUST FOR UNVACCINATED AGENTS#

    #ALT FOCUS 6: PRESENTING PATH JUST FOR VACCINATED AGENTS#
    #log = filter_intervention_by_vacc(log, keep_vacc=True)
    #ALT FOCUS 6: PRESENTING PATH JUST FOR VACCINATED AGENTS#

    #plot_sankey_soft_filter(log, assert_interventions_cleared = False, front_years=20,  tail_years = 6, referral_time_cutoff=8) #20) #8)
    plot_sankey_soft_filter(log, assert_interventions_cleared = False, front_years=20,  tail_years = 10, referral_time_cutoff=8) #20) #8)

        #We wait {referral_time_cutoff/4} years after an agent is administered an intervention until it times-out of a pathway. We further truncate pathways at the intervention before an agent undergoes a primary screening, meaning referral_time_cutoff is not bounded by regular screening callback intervals
        #Meanwhile, the longest acceptable time for an agnet to stay within a single pathway, when all results and followups other than 2nd and 3rd consecutive DNA tests are instantaneous (i.e. occur within a single timestep) is {2*(1 + abandon_followup_invites_threshold)} years, or {8*(1 + abandon_followup_invites_threshold)} timesteps, as followup DNA test invites are sent out with a delay of 1 year, and on top of that we allow for some time until the invite is taken up.
        #Therefore, a primary screening starting at timestep {T}, which takes as long as possible, will continue looking for events at timestep {T + referral_time_cutoff + 8*(1 + abandon_followup_invites_threshold)}. 
        #If we want to not have any 'Simulation Truncated' nodes, we must therefore have additional padding of {referral_time_cutoff + 8*(1 + abandon_followup_invites_threshold)} timesteps, or {referral_time_cutoff/4 + 2*(1 + abandon_followup_invites_threshold)} years, at the end of our log where we introduce no new routine screens in those timesteps.
        #This introudces the inequality {tail_buffer_years >= referral_time_cutoff/4 + 2*(1 + abandon_followup_invites_threshold) } which is necessary  for the existance of 'Simulation Truncated' nodes in a Sankey diagram to be indicative of a bug in the algorithm.
            #-> e.g. with abandon_followup_invites_threshold=0, this reduces to {tail_buffer_years>=referral_time_cutoff/4+2}
        #Furthermore, referral_time_cutoff needs to be able to capture all modelled treatment delays: {referral_time_cutoff >= 4*(1+abandon_followup_invites_threshold)} (each wait is at most this number of timesteps)
        #As there is no benefit to making referral_time_cutoff too large (because why look at longer intervals than are modelled between any two treatments), we can make the optimal equality {referral_time_cutoff=4*(1+abandon_followup_invites_threshold)}. This has the additional benefit that, provided all followups in our log should have an ancestor primary screening in our log (i.e. none are orpahned by primary screen filtering processes), any remaining followups that have not contributed to flows indicate a bug/misunderstanding of the modelled followup times
        #This yields inequality {tail_buffer_years >= 3*(1 + abandon_followup_invites_threshold)}, with {referral_time_cutoff=4*(1+abandon_followup_invites_threshold)} 
            #-> e.g. with abandon_followup_invites_threshold=1, we need referral_time_cutoff:=8 and tail_buffer_years>=6 (with tail_buffer_years=6 for optimality, not generating additional simulated data that is never used!)
        
                #TODO ^ is the reasoning above correct? I think the bound maybe isnt as tight as absolutely needed but still!
    #"""



