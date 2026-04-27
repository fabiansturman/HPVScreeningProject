"""
This file defines the functionality with which we log information from a HPVsim instance.
It can be used stand-alone for tracking the state of a HPVsim simulation over time, or can be used in combination with project_modelling/resultGeneration.py to run a HPVsim simulation and get finalised results from it.



Functions offered (helper functions that are only defined to support the functions listed below are ommitted):
-> makeLoggingIntervention: returns a log dictionary, and 'logger' which can be added to a simulation's intervention list to log simulation state at each timestep
        - each log is a dictionary of states at each timestep, indexed by timestep, i.e. {t: 'simulation state at time t'},
            where 'simulation state at time t' = {name of quantity: value of quantity}. 
            See docstring of log_current_state for info.

-> soft_filter_by_time: filters a log by removing all (outcomes from) interventions other than the primary screens that happened between specified timepoints, and the interventions which occured as a followup from these primary screens (with follow-up defined recursively)
-> filter_intervention_by_partition: for each logged intervention, filters out all outcomes from the intervention for an agent not in the specifed partition (partitions can be agents which are 'alive', 'susceptible', 'cancerous' etc)
    -> filter_intervention_by_cancerous: the partition of included agents are those with/without (as specified in parameter) cancer (by default: checking the cancer state in previous timestep)
    -> filter_intervention_by_cin: the partition of included agents are those in/not in (as specified in parameter) CIN state (by default: checking the CIN state in previous timestep)
    -> filter_intervention_by_vacc: the partition of included agents are those who have/have not (as specified in parameter) been vaccinated 

-> get_timeseries: Extracts a timeseries from a log, given a 'functional_snapshot_extractor' (i.e. a function (log, timepoint)->value we want to track at that timepoint)
    -> get_boolean_partition_size_timeseries: Extracts timeseries of the number of living agents labelled as 'True' for the specified boolean property (where that boolean property needs to be recorded in our log, as in 'log_current_state', e.g. 'alive', 'inactive', etc) 
    -> get_cumulative_doses_timeseries: Extracts timeseries of the total number of HPV vaccine doses (of any type) administered so far
    -> get_vaccinated_proportion_timeseries: Extracts timeseries of the proportion of the adult population that has been administered at least a specified number of doses of some HPV vaccine
    -> get_intervention_timeseries: For the intervention with specified label, returns two timeseries to characterise the intervention over time:
                                            i)  Timeseries of the totla number of intervention usages at each timepointr
                                            ii) Timseries where at each timepoint we have the number of outcomes of each type from this intervention (which, according to a parameter, may be normalised to form a categorical distribution at each timepoint)

-> plot_sankey: plots a Sankey diagram of the flow of agents through the specified interventions, as logged in the provided log
-> plot_sankey_soft_filter: first does a soft-filter of the log to clip away some years from the start and end (as specified in parmeters), with the idea that this will stop us biasing our Sankey diagram with 'incomplete flows' at the end in particular
"""

#Imports
from tqdm import tqdm
import numpy as np
from copy import deepcopy
from math import ceil
import hpvsim as hpv
import matplotlib.pyplot as plt


###-FUNCTIONALITY FOR WRITING A LOG-###
#Logging funcationality - {log} is the dictionary we are logging to (python passes dictionaries by reference)
def log_current_state(sim, log, trackable_interventions_by_label = None): 
    """
    trackable_interventions_by_label: UNION(None, [str]) - either the names (label attributes) of the interventions which we want to track the outcomes of over the simulation, else (if left None) we track all interventions

    This logging schema keeps track of the state of the simulation at each timestep by tracking the events in each timestep (as opposed to saving the full agent state at each timestep, which would take too long and too much space)
    log is a dictonary indexed by timestep as {t: 'data at time t'}, 
        where 'data at time t' is a dictionary with the following keys and values:
            -> n (int): number of agents so far modelled by the sim at this timepoint
            -> alive (np.array[boolean], shape(n,) ): list of boolean values where alive[pid] is the liveliness state of agent pid at this timestep
            -> sex (np.array[int], shape(n,) ): list of int values \in \({0,1}\) where sex[pid] is 0 iff agent pid is a woman
            -> age (np.array[float], shape(n,)): list of floats where age[pid] is the age of agent pid at this timestep
            -> deaths_by_cause: {'other': boolean list (np.array[boolean], shape(n,) ) for whether each agent has died of cause 'other' by completion of this timestep,
                                 'cancer': as above, for cause 'cancer',
                                 'emigrated': as above, for cause 'emigrated' (for HPVsim's purposes, modelled the same as death),
                                 'hiv': as above, for cause 'hiv'}
            -> scale (np.array[float]): list of floats where scale[pid] is the relative scale of agent pid (i.e. if scale=[1,10], then agent 1 represents 10x as many simulated people as agent 0.)

            -> susceptible (np.array[boolean], shape(num_genotypes,n) ): boolean array where susceptible[g,pid] stores, for agent pid, whether the agent is susceptible to (i.e. not infected by) genotype g
            -> infectious (np.array[boolean], shape(num_genotypes,n) ): boolean array where infectious[g,pid] stores, for agent pid, whether the agent is infectious with (i.e. has active infection) genotype g
            -> inactive (np.array[boolean], shape(num_genotypes,n) ): boolean array where inactive[g,pid] stores, for agent pid, whether the agent has an inactive infection with genotype g
           
            -> normal (np.array[boolean], shape(num_genotypes,n) ): boolean array where normal[g,pid] stores, for agent pid, whether the agent has no, or fewer than detectable, abnormal cells due genotype g
            -> cin (np.array[boolean], shape(num_genotypes,n) ): boolean array where cin[g,pid] stores, for agent pid, whether the agent has a detectable number of abnormal cells due to infection from genotype g progressing to the CIN stage, but this has not progressed to cervical cancer
            -> cancerous (np.array[boolean], shape(num_genotypes,n) ): boolean array where cancerous[g,pid] stores, for agent pid, whether the agent has cervical cancer due to an infection with genotype g

            -> vx_doses (np.array[int], shape(n,) ): 'vx_doses[pid]=k' <=> counting any HPV vaccine type, agent pid has recieved k HPV vaccine doses

            -> interventions: {intervention_label: intervention.outcomes} (with a key for each intervention which has recorded outcomes)
                ^note that intervention.outcomes is itself a dictionary {outcome:[list of pid values for people who underwent this outcome]}
                
    """
    track_all_interventions = trackable_interventions_by_label is None

    current_info = {} #'data at time t's
    
    #Record demographic info
    current_info['n'] = sim.n
    current_info['alive'] = np.copy(sim.people.alive)
    current_info['sex'] = np.copy(sim.people.sex) #0 means woman
    current_info['age'] = np.copy(sim.people.age)
    current_info['deaths_by_cause'] = {'other'      : np.copy(sim.people.dead_other),
                                       'cancer'     : np.copy(sim.people.dead_cancer),
                                       'emigrated'  : np.copy(sim.people.emigrated),
                                       'hiv'        : np.copy(sim.people.dead_hiv)
                                       }
    current_info['scale'] = np.copy(sim.people.scale)
    
    #Record viral (i.e. infection) info
    current_info['susceptible'] = np.copy(sim.people.susceptible) #NOTE: I can fully infer susceptible from infectious and inactive; so should I not record this (which will have the most Trues) and then store the other two only, and store them as sparse arrays
    current_info['infectious'] = np.copy(sim.people.infectious)
    current_info['inactive'] = np.copy(sim.people.inactive)

    #Record cell state info
    current_info['normal'] = np.copy(sim.people.normal) #NOTE: I can fully infer normal from cin and cancerous; so should I not record this (which will have the most Trues) and then store the other two only, and store them as sparse arrays
    current_info['cin'] = np.copy(sim.people.cin)
    current_info['cancerous'] = np.copy(sim.people.cancerous)

    #Record vaccination info
    current_info['vx_doses'] = np.copy(sim.people.doses)

    #Record screening pathway info
    intervention_outcomes = {}
    for intervention in sim.interventions:
        if 'outcomes' in dir(intervention) and 'label' in dir(intervention):
            if track_all_interventions or intervention.label in trackable_interventions_by_label: #guards against unwanted tracking of interventions
                intervention_outcomes[intervention.label] = {outcome_key:intervention.outcomes[outcome_key] #possible outcome from this intervention: list of all the pid's who had this outcome from this intervention in this timestep
                                                                for outcome_key in intervention.outcomes.keys()}
    current_info['interventions'] = intervention_outcomes

    #Save all this current information to the log for this timestep
    log[sim.t] = current_info

#Need to define this seperately, so a user can create their logging intervention with everything saving to the correct file each time it is called
def makeLoggingIntervention(trackable_interventions_by_label = None):
    """
    returns (function to pass as an intervention with which to do logging, dictionary info will be logged to)
    """ 
    log = {}
    logger = lambda sim: log_current_state(sim, log=log, trackable_interventions_by_label=trackable_interventions_by_label)
    logger.label="logger" #logger needs a 'label' attribute, so that it can be safely treated as an intervention within the Sim object
    return logger, log

###-FUNCATIONALITY FOR EDITING A LOG-###
def soft_filter_by_time(log, start_timepoint, end_timepoint):
    """
    Returns a copy of {log}, where outcomes from primary screenings ('routine_screening_under50' , 'routine_screening_50andover') AND THE OUTCOMES FROM FOLLOWUP INTERVENTIONS are not recorded before start_timepoint, and outcomes at or beyond end_timepoint at recorded for a given agent pid until pid encounters a primary screen
    This log can be interpreted as one where no screens happen before start_timepoint and none happen at or beyond end_timepoint - and followup interventions are filtered accordingly too
        -> NOTE: interventions can still happen beyond end_timepoint, they just need to be incited by a screening that happened before end_timepoint, hence this is a 'soft' filter
    """
    start_intervention_labels = [ 'routine_screening_under50' , 'routine_screening_50andover' ]
    followup_intervention_labels = list( set(list(log[0]['interventions'].keys())) - set(start_intervention_labels) )

    log = deepcopy(log)

    #Filter the start - iterate from the very first to very last timepoint, ammassing a set of 'allowed pids' (who have had a screning at, or later than, start_timepoint) who get to stay in outcomes - and all others are scrubbed
    allowed_pids = set()
    for t in range(len(log.keys())):
        if t<start_timepoint:
            #Scrub all outcomes from all interventions at this timepoint
            for intervention_label in log[t]['interventions'].keys():
                for outcome in log[t]['interventions'][intervention_label].keys():
                    log[t]['interventions'][intervention_label][outcome] = np.array([])
        else:
            #Update the set of agents who are eligible for tracking (they have had a primary screen at or later than start_timepoint)
            for label in start_intervention_labels:
                for outcome in log[t]['interventions'][label].keys():
                    allowed_pids = allowed_pids.union( set(list(log[t]['interventions'][label][outcome])) )
            #Filter all followup interventions at this step to only include allowed agents (so that followups to primary interventions that happened too early are scrubbed)
            for label in followup_intervention_labels:
                for outcome in log[t]['interventions'][label].keys():
                    outcome_pids = set(list(log[t]['interventions'][label][outcome]))
                    filtered_outcome_pids = outcome_pids.intersection(allowed_pids)
                    log[t]['interventions'][label][outcome] = np.array(list(filtered_outcome_pids))  #must convert the set to a list first, so that if it is empty we get a 1d array with 0 elements (works for iteration) rather than a 0d array (does not work for iteration)

    #Filter the end - iterate from first 'primary screen end' timepoint to the final timepoint in log, keeping track of all the agents whp have had a primary screen at or beyond the 'primary screen end' so far
    banned_pids = set()
    for t in range(ceil(end_timepoint), len(log.keys())):
        #Find which agents to add to the banned-list (and also remove the agents from the primary screenings too)
        for label in start_intervention_labels:
            for outcome in log[t]['interventions'][label].keys():
                banned_pids = banned_pids.union( set(list(log[t]['interventions'][label][outcome])) )
                log[t]['interventions'][label][outcome] = np.array([])
        #Filter all other (i.e. followup) interventions according to the banned-list
        for label in followup_intervention_labels:
            for outcome in log[t]['interventions'][label].keys():
                outcome_pids = set(list(log[t]['interventions'][label][outcome]))
                filtered_outcome_pids = outcome_pids - banned_pids
                log[t]['interventions'][label][outcome] = np.array(list(filtered_outcome_pids)) #must convert the set to a list first, so that if it is empty we get a 1d array with 0 elements (works for iteration) rather than a 0d array (does not work for iteration)

    return log

def filter_intervention_by_partition(log, intervention_label, partition_name, collapser = None, relative_search_zone = [0], progress_bar = False):
    """
    partition_name <- {'alive', 'susceptible', 'infectious', 'inactive', 'normal', 'cin', 'cancerous'}: the quantity we are filtering agents for
    collapser: function to collapse down a row of booleans for a parition defined over genotypes to a single boolean (e.g. needs to be identity for 'alive', needs to be ALL for suscpetible usually, and usually ANY for cancerous - unless we want to filter for a specific genotype)
    relative_search_zone: a list of numbers that say where to look for existance in a partition. e.g. if [-1,0,1], we include an agent in the outcomes of intervention_label at time t iff at time t-1 OR t OR t+1 it is in the parition
    progress_bar (bool): whether to show a tqdm process bar on the timesteps done so far
    
    Returns a copy of {log} where outcomes from {intervention_label} are only recorded if the given agent, at that same timepoint, is in the 'True' parition of {partition_name} (where, if the parition is defined by genotype, we collapse bool^n_genotypes->bool with the provided function )
    """
    log = deepcopy(log)

    if collapser is None:
        collapser = lambda x:x

    for t in log.keys() if not progress_bar else tqdm(log.keys()):
        for outcome in log[t]['interventions'][intervention_label].keys():
            pids = []
            for pid in log[t]['interventions'][intervention_label][outcome]:
                for i in relative_search_zone:
                    if t+i>=0 and t+i<len(log.keys()) and pid<len(collapser(log[t+i][partition_name])) and collapser(log[t+i][partition_name])[pid]: #First three conditions check for validity of timepoint. First and second conditions check the timepoint exists. Second condition to check the agent exists at timepoint t+i (e.g. if i<0, we could query an agent before it is born, which would yield an error): collapser(log[t+i][partition_name]) has one boolean for every agent that exists (or had existed by) timepoint {t+i}.
                        pids.append(pid)
            pids = np.array(list(set(pids))) #remove duplicates
            log[t]['interventions'][intervention_label][outcome] = pids
    
    return log

def filter_intervention_by_cancerous(log, keep_cancerous, intervention_label, relative_search_zone = [-1], progress_bar = False):
    """
    keep_cancerous (bool): if True, filters out all but cancerous agents; if False, filters out only cancerous agents
    relative_search_zone: a list of numbers that say where to look for existance in a partition. e.g. if [-1,0,1], we include an agent in the outcomes of intervention_label at time t iff at time t-1 OR t OR t+1 it is cancerous. Default is 'if the agent was cancerous just before, or at, this timepoint', noting that they can be modelled to test postivie for HPV and follow through to have cancer cleared all in one timestep.
    progress_bar (bool): whether to show a tqdm process bar on the timesteps done so far

    Returns a copy of {log} where outcomes from {intervention_label} are only recorded if the given agent, at that same timepoint, is in the 'True'/'False' (according to 'keep_cancerous') parition of 'cancerous' for any genotype
    """
    if keep_cancerous:
        return filter_intervention_by_partition(log, intervention_label, 'cancerous', lambda bools:np.any(bools, axis=0), relative_search_zone=relative_search_zone, progress_bar=progress_bar)
    else:
        return filter_intervention_by_partition(log, intervention_label, 'cancerous', lambda bools:~np.any(bools, axis=0), relative_search_zone=relative_search_zone, progress_bar=progress_bar)

def filter_intervention_by_cin(log, keep_cin, intervention_label, relative_search_zone = [-1], progress_bar = False):
    """
    keep_cin (bool): if True, filters out all but cin agents; if False, filters out only cin agents
    relative_search_zone: a list of numbers that say where to look for existance in a partition. e.g. if [-1,0,1], we include an agent in the outcomes of intervention_label at time t iff at time t-1 OR t OR t+1 it is cancerous. Default is 'if the agent was cancerous just before, or at, this timepoint', noting that they can be modelled to test postivie for HPV and follow through to have cancer cleared all in one timestep.
    progress_bar (bool): whether to show a tqdm process bar on the timesteps done so far
    
    Returns a copy of {log} where outcomes from {intervention_label} are only recorded if the given agent, at that same timepoint, is in the 'True'/'False' (as specified by 'keep_cin') parition of 'cin' for any genotype
    """
    if keep_cin:
        return filter_intervention_by_partition(log, intervention_label, 'cin', lambda bools:np.any(bools, axis=0), relative_search_zone=relative_search_zone, progress_bar=progress_bar)
    else:
        return filter_intervention_by_partition(log, intervention_label, 'cin', lambda bools:~np.any(bools, axis=0), relative_search_zone=relative_search_zone, progress_bar=progress_bar)

def filter_intervention_by_vacc(log, keep_vacc=True, progress_bar = False):
    """
    progress_bar (bool): whether to show a tqdm process bar on the timesteps done so far

    Returns a copy of {log} where either only agents administered some HPV vaccine in thier lifetime are included in intervention outcomes iff keep_vacc, else the unvaccinated only are included 
    """
    log = deepcopy(log)
    timepoints = sorted(list(log.keys()))
    final_timepoint = timepoints[-1]

    agent_vaccinated = log[final_timepoint]['vx_doses']>0
    keep_agent = agent_vaccinated if keep_vacc else ~agent_vaccinated

    for timepoint in timepoints if not progress_bar else tqdm(timepoints):
        for intervention_label in log[timepoint]['interventions'].keys():
            for outcome in log[timepoint]['interventions'][intervention_label].keys():
                pids = []
                for pid in log[timepoint]['interventions'][intervention_label][outcome]:
                    if keep_agent[pid]:
                        pids.append(pid)
                pids = np.array(pids)
                log[timepoint]['interventions'][intervention_label][outcome] = pids
    
    return log


###-FUNCTIONALITY FOR READING A LOG-###
def get_timeseries(log, functional_snapshot_extractor, use_tqdm = False):
    """
    log: Log from which to extract a timeseries
    functional_snapshot_extractor: as function (log, timepoint)->the value of the quantity we want to track at that timepoint
    use_tqdm: whether to show progress of the iteration through timepoints (may be useful if functional_snapshot_extractor is slow)
    """
    ts = {}
    for t in tqdm(log.keys()) if use_tqdm else log.keys():
        ts[t] = functional_snapshot_extractor(log, t)
    return ts

def get_boolean_partition_size_timeseries(log, parition_name, collapser = None):
    """
    parition_name <- {'alive', 'susceptible', 'infectious', 'inactive', 'normal', 'cin', 'cancerous'}: the quantity we are tracking
    collapser: function to collapse down a row of booleans for a parition defined over genotypes to a single boolean (e.g. needs to be identity for 'alive', needs to be ALL for suscpetible usually, and usually ANY for cancerous - unless we want to filter for a specific genotype)
    
    Returns timeseries of the number of (living) people labelled as part of the positive class of a named boolean partition (e.g. alive/dead, susceptible/not_susceptible)  - collapsed with {collapser} if necessary
    """
    if collapser is None:
        collapser = lambda x:x

    return get_timeseries(log, 
                        lambda l, t:np.dot(np.logical_and(collapser(l[t][parition_name]), l[t]['alive']) , #We only count living agents at each timepoint
                                           l[t]['scale']) #Dot product between [bool] and [float]: True is interpreted as 1.0, False as 0.0. 
    )

def get_cumulative_doses_timeseries(log):
    """
    Returns timeseries of the cumulative number of HPV vaccine doses (of any type) administered 
    """
    return get_timeseries(log, 
                        lambda l, t:np.dot(l[t]['vx_doses'] , l[t]['scale']) #We count all agents that have lived up to each timepoint - not filtering out for living only
    )

def get_vacinnated_proportion_timeseries(log, doseage_requirement:int = 1, progress_bar = False):
    """
    dosage_requirement (int): An agent is considered to have completed their course of vaccination iff they have had at least dosage_requirement doses of (any) HPV vaccine
    progress_bar (bool): If True, outputs a tqdm progress bar when iterating through timesteps

    Returns timeseries of the proportion of the adult population that has been administered at least {doseage_requirement} doses of some HPV vaccine 
    """
    vaccinated_population_sizes = get_timeseries(log,
                                                 lambda l,t:np.dot(np.logical_and(l[t]['vx_doses']>=doseage_requirement, 
                                                                                  np.logical_and(l[t]['age']>=16,l[t]['alive'])), 
                                                                   l[t]['scale'])
                                                 )
    
    
    total_population_sizes = get_timeseries(log, 
                                            lambda l, t: np.dot(np.logical_and(l[t]['age']>=16,l[t]['alive']),
                                                                l[t]['scale'])
                                            )
    
    #vaccinated_population_sizes = get_timeseries(log, #changing to only include those susceptive to cervical cancer
     #                                            lambda l,t:np.dot(np.logical_and(l[t]['sex']==0,
      #                                                                            np.logical_and(l[t]['vx_doses']>=doseage_requirement, 
       #                                                                           np.logical_and(l[t]['age']>=16,l[t]['alive']))), 
        #                                                           l[t]['scale'])
         #                                        )

    #total_population_sizes = get_timeseries(log, 
     #                                       lambda l, t: np.dot(np.logical_and(l[t]['sex']==0,np.logical_and(l[t]['age']>=16,l[t]['alive'])),
      #                                                          l[t]['scale'])
       #                                     )
    
    ts = {}
    for t in vaccinated_population_sizes.keys() if not progress_bar else tqdm(vaccinated_population_sizes.keys()):
        ts[t] = vaccinated_population_sizes[t]/total_population_sizes[t]

    return ts

def get_intervention_timeseries(log, intervention_label, normalize:bool = False):
    """
    intervention_label (str): Name of the intervention whose distribution of outcomes over time we will track
    normalize (bool): If False, each outcome's timeseries is the total number of people with each outcome at each timepoint, if True the timeseries(es) are normalised to form an outcome distribution at each point 
    
    Returns timeseries(es) relating to this intervention (a 2-tuple, (uptake_timeseries, outcome_timeseries)):
        ->info about intervention uptake: a timeseries of the total number of intervention usages at each point
        ->info about intervention outcomes: a dictionary of outcomes containing {n_outcomes} timeseries(es) where intervention {intervention_label} has {n_outcomes} outcomes, tracking the number of people recieving this outcome from this intervention at each timepoint. This dictionary is indexed by outcome label.
    """
    outcomes = list(log[0]['interventions'][intervention_label].keys())
    uptake_timeseries = {t:None for t in log.keys()} #to be populated
    outcome_timeseries = {outcome:None for outcome in outcomes} #to be populated

    #Populate outcome_timeseries 
    for outcome in outcomes:
        #Define a function to get the total number of agents recieving outcome {outcome} at each timepoint
        def get_total_outcome_at_t(l, t):
            outcomed_agents = l[t]['interventions'][intervention_label][outcome] #list of the pids of the agents who recieved the specified outcome
            outcomed_agents_scales = l[t]['scale'][outcomed_agents] #using Numpy 'fancy indexing' to get the specific scales of the 'outcomed' agents
            return sum(outcomed_agents_scales)
        
        outcome_timeseries[outcome] = get_timeseries(log,get_total_outcome_at_t) 

    #Populate uptake_timeseries
    for timepoint in uptake_timeseries.keys():
        uptake_timeseries[timepoint] = sum([outcome_timeseries[outcome][timepoint] for outcome in outcomes])
    
    #If required, normalise at each timepoint by total outcomes of any type
    if normalize:
        for timepoint in log.keys():
            for outcome in outcomes:
                outcome_timeseries[outcome][timepoint] = outcome_timeseries[outcome][timepoint] / uptake_timeseries[timepoint]
    
    return uptake_timeseries, outcome_timeseries

def plot_sankey(log, included_interventions_by_label = None, referral_time_cutoff = 8, assert_interventions_cleared = True, show_progress = True):
    """
    included_interventions_by_label: the set of interventions we are plotting in our Sankey diagram. If None, includes all interventions in the algorithm that I am tracking (hardcoded in)
    referral_time_cutoff (float): Between subsequent interventions A and B for a fixed agent, if the time difference between administration of A and B is at most {referral_time_cutoff}, we consider this agent to have 'flowed' from A to B. Else, the agent does not contribute to any flow between the two interventions (in same unit of time as 'timepoint' index in log; HPVsim default measures in years)
    assert_interventions_cleared (bool): If True, runs a check to ensure all intervention outcomes have been accounted for after populating Sankey diagram. As diagram is populated using forward passes (i.e. populating flows by tracking agents as they flow through the screening algorithm from primary screen), if filtering logs by time in a way that isn't 'soft' (i.e. leaving orphaned followup interventions where the log has no record of the ancestor primary screen which led to it), this will likely fail. It will also fail if the primary screens in the log have been filtered to restrict our population of analysis to, say, agents with cancer in the timestep before a primary screening.
    show_progress (bool): if True, shows a progress bar over timepoints when we iterate over timepoints to compute the Sankey diagram - showing a progress bar over the most costly part of making the diagram

    Plots a Sankey diagram of the flow between interventions in a simulated screening algorithm.
    Notes: 
        -> referral_time_cutoff = 8 as in our modelling, all followups are first available within 4 timesteps (1 year) and {GlobalScreeningParameters.abandon_followup_invites_threshold} is set to 1 year, meaning that any followups which could be modelled will have occured within 8 timesteps of the inciting previous intervention
        -> Agents may develop cancer while flowing from intervention A to B, thereby changing their scale. For consistency across our diagram, we each agent to have a fixed scale throughout the period of their flow on the diagram, determined by their scale at the final timepoint captured in the diagram. For example, an agent that flows (start, scale=1->)A->B->C(->end, scale=10) will be considered to have had scale=10 throughout, that is, cancer throughout. This is justified as (i)the length of a continuous flow can be at most a few years, so few relatively few cancers will develop during flow, especially because (ii) broadly speaking, those agents which develop cancer and made it through to even 2 interventions will have done so as a result of cancer and would have had cancer the whole time, again so the amount of distortion is small.
            #NOTE: this biases the scale of flows across the diagram towards agents with cancer but is not a biggie as the diagram is only to check the right flows are happening and not meant for any deeper conclusions than debugging/demonstrating the implementation works
    """
    import plotly.graph_objects as go

    #Populate list of interventions we are tracking (if needed)
    if included_interventions_by_label is None:
        included_interventions_by_label = ['routine_screening_under50', 'routine_screening_50andover', 'first_cytology', 'second_consecutive_screening', 'second_cytology', 'third_consecutive_screening', 'third_cytology', 'colposcopy', 'ablation', 'general_cancer_treatment']

    #Copy out the timeseries of interventions. intervention_log[intervention_label][t] := "dictionary of outcomes for this intervention at this timepoint"
    timepoints = list(log.keys())
    intervention_log = {intervention_label: {t: {outcome:list(log[t]['interventions'][intervention_label][outcome])
                                                 for outcome in log[t]['interventions'][intervention_label].keys()}
                                             for t in timepoints} 
                        for intervention_label in included_interventions_by_label}

    #Define node properties for plotting - nodes alternate in layers between interventions and outcomes
        #node_labels consists of hard-coded names for each potential node in the Sankey diagram (to correspond with interventions and outcomes defined in the logged simulation)
    node_labels = ['routine_screening_under50', 'routine_screening_50andover', 'routine_screening', #level 1
                    'first_cytology',#level 3
                    'second_consecutive_screening', #level 5
                    'second_cytology', #level 7
                    'third_consecutive_screening', #level 9
                    'third_cytology', #level 11
                    'colposcopy', #level 13
                    'ablation', 'general_cancer_treatment'#level 15
                    ] + [
                    "screening:positive", "screening:inadequate", "screening:negative", #1st screening outcomes, level 2
                    "first_cytology:abnormal", "first_cytology:ascus", "first_cytology:inadequate", "first_cytology:normal", #1st cytology outcomes, level 4
                    "second_consecutive_screening:positive", "second_consecutive_screening:inadequate", "second_consecutive_screening:negative", #2nd screening outcomes, level 6
                    "second_cytology:abnormal", "second_cytology:ascus", "second_cytology:inadequate", "second_cytology:normal",#2nd cytology outcomes, level 8
                    "third_consecutive_screening:positive", "third_consecutive_screening:inadequate", "third_consecutive_screening:negative",#3rd screening outcomes, level 10
                    "third_cytology:abnormal", "third_cytology:ascus", "third_cytology:inadequate", "third_cytology:normal",#3rd cytology outcomes, level 12
                    "colposcopy:cancer", "colposcopy:hsil", "colposcopy:lsil", "colposcopy:ascus", "colposcopy:normal",#colposcopy outcomes, level 14
                    "ablation:successful", "ablation:unsuccessful",#ablation outcomes, level 16
                    "general_cancer_treatment:cleared", "general_cancer_treatment:extended", "general_cancer_treatment:failed" #general cancer treatment outcomes,level 16
                    ] + [
                        'Simulation Truncated', 'Left Pathway',#level 3
                        'Simulation Truncated', 'Left Pathway',#level 5
                        'Simulation Truncated', 'Left Pathway',#level 7
                        'Simulation Truncated', 'Left Pathway',#level 9
                        'Simulation Truncated', 'Left Pathway',#level 11
                        'Simulation Truncated', 'Left Pathway',#level 13
                        'Simulation Truncated', 'Left Pathway',#level 15
                        'Simulation Truncated', 'Left Pathway',#level 17
                    ]
        #get_node_from_label allows us to map internal representations of nodes (eitehr an intervention label, an intervention label alongside the itnervention outcome, or a way to leave the pathway at a particular level) to the corresponding node's index in node_labels 
    get_node_from_label =  {#'routine_screening_under50':0, 'routine_screening_50andover':1,
                            'routine_screening_under50':2, 'routine_screening_50andover':2, #UNCOMMENT THIS TO TREAT ALL INITIAL ROUTINE SCREENS THE SAME
                            'first_cytology':3,
                            'second_consecutive_screening':4, 
                            'second_cytology':5, 
                            'third_consecutive_screening':6,
                            'third_cytology':7, 
                            'colposcopy':8, 
                            'ablation':9, 'general_cancer_treatment':10,
                            "routine_screening_under50:positive":11, "routine_screening_under50:inadequate":12, "routine_screening_under50:negative":13, 
                            "routine_screening_50andover:positive":11, "routine_screening_50andover:inadequate":12, "routine_screening_50andover:negative":13, 
                            "first_cytology:abnormal":14, "first_cytology:ascus":15, "first_cytology:inadequate":16, "first_cytology:normal":17,
                            "second_consecutive_screening:positive":18, "second_consecutive_screening:inadequate":19, "second_consecutive_screening:negative":20,
                            "second_cytology:abnormal":21, "second_cytology:ascus":22, "second_cytology:inadequate":23, "second_cytology:normal":24,
                            "third_consecutive_screening:positive":25, "third_consecutive_screening:inadequate":26, "third_consecutive_screening:negative":27,
                            "third_cytology:abnormal":28, "third_cytology:ascus":29, "third_cytology:inadequate":30, "third_cytology:normal":31,
                            "colposcopy:cancer":32, "colposcopy:hsil":33, "colposcopy:lsil":34, "colposcopy:ascus":35, "colposcopy:normal":36,
                            "ablation:successful":37, "ablation:unsuccessful":38,
                            "general_cancer_treatment:cleared":39, "general_cancer_treatment:extended":40, "general_cancer_treatment:failed":41, 
                            'Simulation Truncated 3':42, 'Left Pathway 3':43,
                            'Simulation Truncated 5':44, 'Left Pathway 5':45,
                            'Simulation Truncated 7':46, 'Left Pathway 7':47,
                            'Simulation Truncated 9':48, 'Left Pathway 9':49,
                            'Simulation Truncated 11':50, 'Left Pathway 11':51,
                            'Simulation Truncated 13':52, 'Left Pathway 13':53,
                            'Simulation Truncated 15':54, 'Left Pathway 15':55,
                            'Simulation Truncated 17':56, 'Left Pathway 17':57,
    }
        #set x-coordinate for each node, according to their 'trophic level'
    xs = [1/17, 1/17, 1/17, #level 1
                    3/17,#level 3
                    5/17, #level 5
                    7/17, #level 7
                    9/17, #level 9
                    11/17, #level 11
                    13/17, #level 13
                    15/17, 15/17#level 15
                    ] + [
                    2/17, 2/17, 2/17, #1st screening outcomes, level 2
                    4/17, 4/17, 4/17, 4/17, #1st cytology outcomes, level 4
                    6/17, 6/17, 6/17, #2nd screening outcomes, level 6
                    8/17, 8/17, 8/17, 8/17,#2nd cytology outcomes, level 8
                    10/17, 10/17, 10/17,#3rd screening outcomes, level 10
                    12/17, 12/17, 12/17, 12/17,#3rd cytology outcomes, level 12
                    14/17, 14/17, 14/17, 14/17,14/17,#colposcopy outcomes, level 14
                    16/17, 16/17,#ablation outcomes, level 16
                    16/17, 16/17, 16/17 #general cancer treatment outcomes,level 16
                    ] + [
                        3/17, 3/17-1/34,#level 3
                        5/17, 5/17-1/34,#level 5
                        7/17, 7/17-1/34,#level 7
                        9/17, 9/17-1/34,#level 9
                        11/17, 11/17-1/34,#level 11
                        13/17, 13/17-1/34,#level 13
                        15/17, 15/17-1/34,#level 15
                        17/17, 17/17-1/34,#level 17
                    ]
        #set y-coordinate for each node, according initially to number of nodes in each 'trophic level', then adjusted as necessary
    ys = [1/3, 2/3, 1/2, #level 1
                    2/3,#level 3
                    1/2, #level 5 - second screening
                    2/3, #level 7
                    1/2, #level 9
                    2/3, #level 11
                    0.85, #level 13
                    8/10, 9/10#level 15
                    ] + [
                    3/4, 3/3, 1/4, #1st screening outcomes, level 2
                    9/10,2/4,3/4, 1/2, #1st cytology outcomes, level 4
                    2/3, 3/3, 1/3, #2nd screening outcomes, level 6
                    8/10,2/4,3/4, 1/2,#2nd cytology outcomes, level 8
                    2/3, 3/3, 1/3,#3rd screening outcomes, level 10
                    2/3+0.05,2/4,3/4, 2/3-0.05,#3rd cytology outcomes, level 12
                    9/10,8/10,3/5, 4/5, 7/10,#colposcopy outcomes, level 14
                    7/10,8/10,#ablation outcomes, level 16
                    9/10, 4/5, 0.95 #general cancer treatment outcomes,level 16
                    ] + [
                        0.01+0/3, 0.18,#level 3
                        0.01+0/3, 0.18,#level 5
                        0.01+0/3, 0.18,#level 7
                        0.01+0/3, 0.18,#level 9
                        0.01+0/3, 0.18,#level 11
                        0.01+0/3, 0.18,#level 13
                        0.01+0/4, 0.18,#level 15
                        0.01+0/2, 0.18,#level 17
                    ]
     #node colors
    colors_by_node_type = {
        'screening':'cornflowerblue',
        'cytology':'yellow',
        'colposcopy':'mediumvioletred',
        'ablation':'pink',
        'generalcancertreatment':'lawngreen',

        'goodoutcome':'green',
        'mediumoutcome':'orangered',
        'pooroutcome':'red',

        'leftscreeningpathway':'cyan',
    
        'unexpected':'black' #for nodes we don't expect to see, colour them black to easily show if something is wrong
    }
    node_colors = [colors_by_node_type['screening'], colors_by_node_type['screening'], colors_by_node_type['screening'], #level 1
                    colors_by_node_type['cytology'],#level 3
                    colors_by_node_type['screening'], #level 5
                    colors_by_node_type['cytology'], #level 7
                    colors_by_node_type['screening'], #level 9
                    colors_by_node_type['cytology'], #level 11
                    colors_by_node_type['colposcopy'], #level 13
                    colors_by_node_type['ablation'], colors_by_node_type['generalcancertreatment']#level 15
                    ] + [
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'], colors_by_node_type['goodoutcome'], #1st screening outcomes, level 2
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'], #1st cytology outcomes, level 4
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'], 
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'], 
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'], 
                    colors_by_node_type['pooroutcome'], colors_by_node_type['unexpected'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'], 
                    colors_by_node_type['pooroutcome'], colors_by_node_type['mediumoutcome'], colors_by_node_type['unexpected'], colors_by_node_type['unexpected'],  colors_by_node_type['goodoutcome'],#colposcopy outcomes, level 14
                     colors_by_node_type['goodoutcome'], colors_by_node_type['pooroutcome'],#ablation outcomes, level 16
                     colors_by_node_type['goodoutcome'], colors_by_node_type['unexpected'], colors_by_node_type['pooroutcome']#general cancer treatment outcomes,level 16
                    ] + [
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 3
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 5
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 7
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 9
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 11
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 13
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 15
                        colors_by_node_type['unexpected'], colors_by_node_type['leftscreeningpathway'],#level 17
                    ]
    
    #Define unpopulated edge properties for populating (i.e. start with no edges)
    edge_value_by_source_and_target = {} #dictionary of form (source index, target index): value, so that we can accumulate over values and easily work out whether we need to add to an existing dictionary entry, or make a new one, while following the paths of agents

    #Capture all left->right movement that starts with 'routine_screening_under50' or 'routine_screening_50andover'  with a forward pass
    """
    We will iterate timestep-by-timstep (counter: i (in the code, 'timepoint')). At each timestep:
        -> the agents of interest are those who have an outcome in one of these two intiial routine screenings
        -> make a list of the agents of interest
        -> iterate over each agent of interest 
            -> make a list node indices to trace out the path this agent takes; instantiate it with the index for the kind of routine screening this agent started with followed by the index for the outcome
                ->then delete this agent from the relevant outcome list (to avoid double counting in the future)
                -> initialise time_since_last_intervention to 0
            ->  We use a while loop to follow the agent over time (while time_since_last_intervention<=referral_time_cutoff)
                ->in each iteration, move to next timestep (a counter j over interval timesteps, initialised at whatever i our overall timestep counter is, would work) and update time_sice_last_intervention with teh change in timestep
                ->iterate over all interventions at this timepoint, within each intervention iterating through all its outcome lists to try and find the current agent.
                    ->if we do not find the current agent, we can move on to the next iteration in the loop
                    ->if we do find the current agent, we add the indecies of the nodes corresponding to the intervention and the outcome to our little list of nodes for this agent, and delete this agent from the list of outcomes that we just found it in, and we update time_since_last_intervention to 0
            -> once we leave the while loop:
                -> do a sanity check that our outcome list isnt too long (look at the algoroithm to find a reasonable upper bound), to check we havent been counting non-followups as followups, or just implemented something wrong somehwere
                -> add index   len(node_labels)-1    to the end of our list of nodes, to show the agent has left the system for this round
                -> search between timepoints i and j (i.e. the timepoints bounding this agent's flow that we are currently considering) for its largest scale in this time, and record that to a variable 'agent flow size'. Our Sankey diagram will be internally consistent as long as we have consistnet flow for each agent over its flow through the diagram, otherwise things won't add up
                ->for each consecutive pair in our list of nodes (pairs overlap, of course), see if the pair already is a key in edge_value_by_source_and_target. If so, add our saved agent flow size to it, else, intiialise a new entry in the dictionary with this agent flow size

    By deleting entries from our log at each timestep (noting we are using a copy of the log so we aren't losing information), we know if we have missed some flow - if any agents remain in any outcomes
    
    WE NOTE THAT this approach imposes no assumptions on how the algorithm should move agents between interventions, it simply keeps track of all agents for {referral_time_cutoff} timesteps after the most recent intervention it was involved in, and says the agent has 'flowed' from one intervention to the next if it has reached that next intervention within that number of timesteps (e.g. if referral_time_cutoff==3, B can happen just after A in the same timestep, or up to 3 timesteps after A (inclusive) and we will define in all cases some flow from A to B).
        -> this has the advantage that we can know for sure that agents are flowing in the right ways between interventions, by confirming that there are no flows that violate what we intend the screening algorithm to be doing. ]
        -> to avoid long referral-time-cutoffs causing us to bleed into discovering outcomes from the next screening algorithm, and therefore getting misleading loops in our Sankey diagram, when we discover that there is an outcome from a primary intervention (i.e. value in start_intervention_labels) in our forward sweep, we prematurely leave the while loop marking it as 'agent left pathway'
    """
    start_intervention_labels = [ 'routine_screening_under50' , 'routine_screening_50andover' ]
    #PRECONDITION: interventions are stored in the sim s.t. members of 'start_intervention_labels' are all ordered before any other 'included_interventions_by_label'.
        #^ this precondition is required to ensure when we forward-pass from an agent's routine screening through the other interventions, we know to stop when we hit another routine screening for the same agent (rather than just append all the routine screenings to one misleading and unsightly looping flow)
    #PRECONDITION (supercedes precondition above): included_interventions_by_label is defined in the (not strict) order we expect them to occur 
        #^ this precondition is required to ensure that when multiple interventions happen to an agent sequentially in the same timestep, we introduce the flows in the correct order
    
    iterator = range(len(timepoints))
    for timepoint_index in iterator if not show_progress else tqdm(iterator):
        timepoint = timepoints[timepoint_index]

        #Iterate over the interventions that are the global sources of flow in the Sankey diagram - (at a given timepoint) we complete the flow for all agents starting in one intervention before beginning to deal with the next intervention
        for intervention_label in start_intervention_labels:
            intervention_node = get_node_from_label[intervention_label]#  node_labels.index(intervention_label) #stores the name of the node correspoding to this intervention

            #Collect the agents of interest into a list
            pids = []
            for outcome in intervention_log[intervention_label][timepoint].keys():
                pids += intervention_log[intervention_label][timepoint][outcome]
            pids = list(set(pids)) #get unique pids

            #Iterate over each agent of interest, to add their contribution to the flow. Note this may not be an agent's sole contribution to the flow over their lifetime, as if they are primary-screened at a different timepoint then they will also travel through the screening pathway, and contribute to the flow, then - and that will be accounted for when {timepoint} is a different value)
            for pid in pids:
                node_path = [] #node_path will store the agent's path through Sankey diagram nodes (alternating between interventions and outcomes)
                corresponding_timepoints_for_path = [timepoint, timepoint] #this is recorded to output debug information if we generate a node_path that is longer than expected to be possible

                #Add the initial two nodes corresponding to the routine screening and its outcome to node_path
                node_path.append(intervention_node)
                for outcome in intervention_log[intervention_label][timepoint].keys():
                    if pid in intervention_log[intervention_label][timepoint][outcome]:
                        node_path.append(get_node_from_label[f"{intervention_label}:{outcome}"]) 
                        intervention_log[intervention_label][timepoint][outcome].remove(pid)
                            #NOTE: ASSUMPTION + SANITY CHECK (when we see the diagram we can see if the assumption holds, which it should) we are assuming that any given pid cannot, in the same intervention at the same timestep, appear in two different outcome lists at the same time. If it can, then this structure will introduce flow in the Sankey diagram directly between outcomes (in the arbitary direction the outcomes are iterated over here).  

                #Sweep over this and later timepoints to find any later interventions that follow on from this intial screening (may be at same timepoint, e.g. a cytology should follow a screening in the same timepoint)
                sweeping_timepoint_index = timepoint_index #if sweeping_timepoint_index == len(timepoints), then our loop will not run because we have reached the end of the simulation
                time_since_last_intervention = 0
                new_primary_screen_found = False


                while sweeping_timepoint_index<len(timepoints) and time_since_last_intervention <= referral_time_cutoff and not new_primary_screen_found:
                    """
                    At the start of each while loop iteration, upon satisfaction of the condition:
                        -> sweeping_timepoint_index == "the next timepoint that has been unexplored", and this timepoint exists (i.e. we arent at the end of the timepoints)
                        -> time_since_last_intervention is the time between the last intervention that pid underwent and the current sweeping_timepoint_index 
                    """
                    sweeping_timepoint = timepoints[sweeping_timepoint_index]
                    intervention_found = False

                    #Iterate over all outcomes of all interventions to find if pid has been administered an intervention at this (sweeping) timepoint
                    for candidate_intervention_label in included_interventions_by_label: #start_intervention_labels needs to be included in included_interventions_by_label
                        for candidate_outcome in intervention_log[candidate_intervention_label][sweeping_timepoint].keys():
                            if pid in intervention_log[candidate_intervention_label][sweeping_timepoint][candidate_outcome]:
                                if not new_primary_screen_found and candidate_intervention_label not in start_intervention_labels:
                                    #In this case, we have found an intervention that leads on from the last intervention, and we have not yet run into a primary screen on this sweeping timestep, not is it a primary screen itself - so we add to this flow
                                    intervention_found = True
                                    node_path.append(get_node_from_label[candidate_intervention_label])
                                    node_path.append(get_node_from_label[f"{candidate_intervention_label}:{candidate_outcome}"])
                                    corresponding_timepoints_for_path.append(sweeping_timepoint); corresponding_timepoints_for_path.append(sweeping_timepoint)

                                    intervention_log[candidate_intervention_label][sweeping_timepoint][candidate_outcome].remove(pid)
                                else: 
                                    if candidate_intervention_label == intervention_label and sweeping_timepoint_index == timepoint_index:
                                        pass #In this case, we have just found the intervention that started off this flow, so can safely ignore (else we get wrongly looping flows around out initial screens!)
                                    else:
                                        #In this case, we have found that the agent has had a new primary screen, and all following interventions should be atributed to that - and will be caught when we increase timepoint_index to it, and so we can safely stop this sweeping while loop
                                        new_primary_screen_found = True 
                                            #NOTE: this statment could run many times after a new primary screen is found - one repeat for each additional intervention in the same timestep, coming off the new screen - but it is OK to set the same thing to True many times
                                    
                                    #NOTE: ASSUMPTION + SANITY CHECK (when we see the diagram we can see if the assumption holds, which it should) we are assuming that any given pid cannot have two outcomes from the same intervention at the same timestep (as above)

                    if not new_primary_screen_found:
                        #If we have not run into a new primary screen undergone by the agent at this timestep, update loop variables normally
                        sweeping_timepoint_index+=1
                        if intervention_found:
                            time_since_last_intervention = 0
                        else: 
                            if sweeping_timepoint_index<len(timepoints):
                                time_since_last_intervention += timepoints[sweeping_timepoint_index] - timepoints[sweeping_timepoint_index-1] 
                            else:
                                time_since_last_intervention = None

                #We have left the while loops: either (i) sweeping_timepoint_index==len(timepoints) meaning we have reached the end of the simulation, or (ii) time_since_last_intervention >referral_time_cutoff meaning the agent has left the followups tpo return to regular screening
                # ... Therefore, we find the final level and exit at the next level's node (found by finding the x coordinate of the final node, and multiplying by 17)
                final_level = round(xs[node_path[-1]]*17)
                if sweeping_timepoint_index==len(timepoints) and not new_primary_screen_found:
                    print("outputting a truncated path:")
                    for index in range(len(node_path)):
                        print(f"Node {node_labels[node_path[index]]} at timepoint {corresponding_timepoints_for_path[index]}")
                    print("------")

                    node_path.append(get_node_from_label[f"Simulation Truncated {final_level+1}"])
                else:  #The agent has left the screening pathway iff (new_primary_screen_found OR time_since_last_intervention > referral_time_cutoff)
                    node_path.append(get_node_from_label[f"Left Pathway {final_level+1}"])

                #Sanity check that our outcome list is not too long
                if len(node_path)>17:
                    print("A path is too long:")
                    for node in node_path:
                        print(node_labels[node])
                    print("----")
                assert len(node_path) <= 17 #3x [screening, cytology] + [colposcopy]+ [ablation OR general cancer treatment] + [left pathway OR simuilation truncated], then x2 as we need nodes for both the interventions and the outcomes

                #Find agent scale
                max_scale = 0
                for tp in timepoints[timepoint_index:sweeping_timepoint_index]:
                    scale = log[tp]['scale'][pid]
                    if scale>max_scale:
                        max_scale = scale
                
                #Add to revelant flows in the Sankey diagram itself
                for i in range(len(node_path)-1):
                    source = node_path[i]; target = node_path[i+1]
                    if (source, target) in edge_value_by_source_and_target.keys():
                        edge_value_by_source_and_target[(source, target)] += max_scale
                    else:
                        edge_value_by_source_and_target[(source, target)] = max_scale

    #Sanity check that we have accounted for all intervention information we wanted to track (if we want to perform this check)
    for intervention_label in intervention_log.keys():
        for timepoint in intervention_log[intervention_label].keys():
            for outcome in intervention_log[intervention_label][timepoint].keys():
                if len(intervention_log[intervention_label][timepoint][outcome]) > 0:
                    print(f"THIS INTERVENTION REMAINS: intervention_log[{intervention_label}][{timepoint}][{outcome}] = {intervention_log[intervention_label][timepoint][outcome]}")
                    if assert_interventions_cleared:
                        assert False #failed assertion that all interventions should be cleared off!
                
    #Format edge info such that it can be read by the Sankey diagram
    sources = []
    targets = []
    values = []
    for source, target in edge_value_by_source_and_target.keys():
        sources.append(source)
        targets.append(target)
        values.append(edge_value_by_source_and_target[(source, target)])

    #Format coordinates - our xs and ys need to be filtered to remove any nodes from the list that are not in the diagram, as coordinates refer to presented nodes only
    included_nodes = list(set(sources).union(set(targets)))
    xs = list(np.array(xs)[included_nodes]) #use Numpy fancy indexing to get just the indices we want
    ys = list(np.array(ys)[included_nodes])

    #Plot Sankey diagram itself
    fig = go.Figure(data=[go.Sankey(
        arrangement="freeform",
        #arrangement="fixed",
        node=dict(
            pad=20,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=["" for _ in node_labels], #dont print labels by default - code below makes them appear when you hover over them
            customdata = node_labels,
            hovertemplate="%{customdata}",
            x = xs,
            y= ys,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            customdata = [f"{node_labels[sources[i]]} -> {node_labels[targets[i]]}" for i in range(len(sources))],
            hovertemplate=(
        "%{customdata}<br>"
        "Value: %{value}<br>"
            "<extra></extra>"
    )
        )
    )])

    fig.update_layout(title_text="My Sankey Diagram", font_size=12)
    fig.show()

def plot_sankey_soft_filter(log, front_years, tail_years, dt=0.25, **kwargs):
    """
    Plots Sankey diagram from provided {log}, after soft-filtering {front_years} from the start and {tail_years} years from the end of the log (to allow for completion of followups in the diagram+analysis)
    dt: timestep of sim (0.25 by default in HPVsim 2.2.4)
    """
    log = soft_filter_by_time(log, front_years/dt, len(log.keys())-tail_years/dt) 
    plot_sankey(log, **kwargs)


###-DEMO (AND TESTING) CODE FOUND IN testingSimulationsUsingLogging.py-###