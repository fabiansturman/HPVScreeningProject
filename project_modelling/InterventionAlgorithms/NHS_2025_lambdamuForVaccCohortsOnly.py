"""
This file contains an model of the current NHS England Cervical Screening Algorithm, with adjustments to allow for changing regular screening callback intervals to differ between individuals according to whether their 'cohort' was offered vaccination at all (inc just girls in their cohort being offered vaccination) or not
The algorithm switches from the current NHS England Cervical Screening Algorithm to the amended one at the switch year (2026).
    ^https://www.gov.uk/government/publications/cervical-screening-care-pathway/cervical-screening-care-pathway (version as of 1st July 2025)

File usage:
    -> Import this file to easily integrate screening algorithm model into a HPVsim simulation
        -> Use the function get_interventions(  l, #the screening scale factor for cohorts not offered any vaccine
                                                m, #the screening scale factor fro cohorts offered some vaccine
                                                end_younger_screening_at_switch_year=False #whether to stop all screening beyond the switch year (2026) for all cohorts offered (or at least partially offered, that is, when it was just girls offered) a vaccine
                                                ):
    -> Run this file (as __main__) to run a HPVsim simulation with this screening and debug changes to the algorithm
"""

#--- Imports ---#
import hpvsim as hpv
import numpy as np
from copy import copy

#NOTE: the imports below are a little unstable, as they change depending on whether we are running things from within the directory InterventionAlgorithms, or outside of it.

from .generalCancerTreatment import  GeneralCancerTreatment, treat_num_cancer, treat_num_canceragnostic  #USE THIS iff running code which imports this from outside of Intervention Algorithms
from InterventionAlgorithms import GlobalScreeningParameters as GlobalScreeningParameters  #USE THIS iff running code which imports this from outside of Intervention Algorithms

#from generalCancerTreatment import  GeneralCancerTreatment, treat_num_cancer, treat_num_canceragnostic #USE THIS iff running this code itself, or other code which is inside and imports this code
#import GlobalScreeningParameters as GlobalScreeningParameters #USE THIS iff running this code itself, or other code which is inside and imports this code


def init_intervention_trackers(sim):
    """
    For t=0, sets up tracking apparatus within the sim for other interventions in this file to run.
    WARNING: Must be added to a sim before all other interventions defined in this file. 
    """
    if sim.t == 0:
        #Tracks properties of agents as dictionaries {agent_id: tracked property} to handle dynamic populations; agent_id is permanent for any given agent
        sim.people.last_hpv_result = {} #track whether an agent's last hpv test result was positive or not
            #TODO: could I refactor this s.t. I just keep track of a list of agent_id's whose last hpv test results were positive, as we dont need to clog things up with keeping track of who has had negative results
        sim.people.needs_consec_screening_2 = {} #keeps track of who needs a second round of consecutive screening, if contains (pid: t) then agent pid needs screening at earliest time t
        sim.people.needs_consec_screening_3 = {} #as above, but for a third round of consecutive screening
init_intervention_trackers.label="init_intervention_trackers"


#Custom interventions to update attributes that keep track of the other interventions' results, persistently across years
def update_last_hpv_result(sim):
    #For everyone who has undergone a screening of some sort, update this result
    hpv_test_names = ["routine_screening_under50","routine_screening_50andover", "second_consecutive_screening", "third_consecutive_screening"]
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
        sim.people.needs_consec_screening_2.pop(i)
    for i in sim.get_intervention("second_consecutive_screening").outcomes['negative']:
        sim.people.needs_consec_screening_2.pop(i)
    for i in sim.get_intervention("second_consecutive_screening").outcomes['inadequate']:
        sim.people.needs_consec_screening_2.pop(i)
update_needs_consec_screening_2.label="update_needs_consec_screening_2"

def update_needs_consec_screening_3(sim, name="second_cytology"):
    #If a woman has a -ve result from her second cytology, she is eligible for a third consecutive screening in 12 months
    for i in sim.get_intervention(name).outcomes['normal']: #iterate over all negative outcomes from first cytology
        sim.people.needs_consec_screening_3[i] = sim.t+1/sim["dt"]
        
    #All women who have undergone their second conseucutive screening have this updated in their record (no more need for 2nd consecutive screening)
    for i in sim.get_intervention("third_consecutive_screening").outcomes['positive']:
        sim.people.needs_consec_screening_3.pop(i)
    for i in sim.get_intervention("third_consecutive_screening").outcomes['negative']:
        sim.people.needs_consec_screening_3.pop(i)
    for i in sim.get_intervention("third_consecutive_screening").outcomes['inadequate']:
        sim.people.needs_consec_screening_3.pop(i)
update_needs_consec_screening_3.label="update_needs_consec_screening_3"


#NOTE: routine screening is split into two eligibility functions and two interventions, to target under 50s and 50-and-overs seperately, due to their differently-modelled probabilities of screening uptake (this models the split easily as two instances of the same screening intervention, which then feed into the same followup process) 

# Setting up routine screening, according to last HPV screening outcome (if any), time since last screening, and age
def routine_screen_eligible_under50(sim, l, m, 
                                    switch_year = GlobalScreeningParameters.switch_year, 
                                    end_younger_screening_at_switch_year = False):
    #l = screning interval for cohorts not offered vaccine, m=screening interval for cohorts offered vaccine
    
    #A woman under 50 is eligible for routine screening {3 IF last hpv test +ve else 5} years after their most recent screening
    for_screening = np.isnan(sim.people.date_screened) & (sim.people.age>=24) & (sim.people.age<50) #eligible for screening if they have not yet been screened

    #add further eligibility for those who have been screened before, but are due a followup
    for i in range(sim.n): #sim.n==len(sim.people)
        if i in sim.people.last_hpv_result.keys(): # check if the agent has been screened before
            last_hpv_result = sim.people.last_hpv_result[i]
            age = sim.people.age[i]
            date_screened = sim.people.date_screened[i]

            skip_screening = False #Track whether we want to skip screening for this person because we aren't screening anyone from their cohort beyond current timepoint

            #Determine the screening interval scale factor for this individual according to the year (i.e. whether we are beyond the switch year) and their vaccination status; OR whether to skip screening altogether for this individual 
            if sim.t>(switch_year-sim['start'])/sim['dt']:
                #If we are at/beyond the year where we switch, the person's screening interval then depends on whether when they were 12 (whether it was (2008 or later), or not)
                current_year = sim.t*sim['dt']+sim['start'] #NOTE: for switch_year==2026, the current year will be >2008
                year_when_12 = current_year - (age-12)

                if year_when_12 >= 2008:
                    #verified this runs, but after 3 runs it bricks. but each run of this should get the endyoungerscsnrenefj to be True and skip screening? so why false
                    if end_younger_screening_at_switch_year:
                        #We arent screening such young people from 2026 onwards, so mark as to be skipped - guarded by year_when_12>=2008 s.t. beyond the switch year we only are not screening those who had the oppertunity to be vaccinated as a child
                        skip_screening = True
                    else:
                        #In this case, the individual was in a cohort offered vaccines, and we are beyond the switch year, so rescale by m (for cohorts offered vaccine)
                        screen_interval_scale_factor = m 
                else:
                    #In this case we are at the switch year, but the individual was not in a cohort offered a vaccine (was aged 12 too early, before 2008), so rescale by l
                    screen_interval_scale_factor = l 
            else:
                #If we have not yet reached the year at which we switch, we do standard NHS 2025 screening algoorithm; multiply all screening intervals by 1 to keep them the same
                screen_interval_scale_factor = 1


            if (24<=age<50) and (not skip_screening): #restrict by age (50-and-overs get dealt with by the other regular screening intervention) AND only set for_screening to True for this agent if we are not to skip screening for them
                screen_interval = 3 if last_hpv_result==1 else 5*screen_interval_scale_factor

                if sim.t>=date_screened+screen_interval/sim['dt']:
                    for_screening[i] = True
            
    return for_screening





def routine_screen_eligible_50andover(sim, l, m, 
                                      switch_year = GlobalScreeningParameters.switch_year,
                                      end_younger_screening_at_switch_year = False):
    #A woman 50 or over is eligible for a routine screening 5 years after their most recent screening
    for_screening = np.isnan(sim.people.date_screened) & (sim.people.age>=50) & (sim.people.age<=65) #eligible for screening if they have not yet been screened


    #add further eligibility for those who have been screened before, but are due a followup
    for i in range(sim.n): #sim.n==len(sim.people)
        if i in sim.people.last_hpv_result.keys(): # check if the agent has been screened before (this implicitly filters out agents too young to be screened, so we only need to check age<65 later)
            age = sim.people.age[i]
            date_screened = sim.people.date_screened[i]

            skip_screening = False #Track whether we want to skip screening for this person because we aren't screening anyone beyond current timepoint

            #Determine the screening interval scale factor for this individual according to the year (i.e. whether we are beyond the switch year) and their vaccination status; OR whether to skip screening altogether for this individual 
            if sim.t>(switch_year-sim['start'])/sim['dt']:
                #If we are at/beyond the year where we switch, the person's screening interval then depends on whether when they were 12 (whether it was 2008 or later)
                current_year = sim.t*sim['dt']+sim['start'] #NOTE: for switch_year==2026, the current year will be >2008
                year_when_12 = current_year - (age-12)

                if year_when_12 >= 2008:
                    if end_younger_screening_at_switch_year:
                        #We arent screening such young people from 2026 onwards, so mark as to be skipped - guarded by year_when_12>=2008 s.t. beyond the switch year we only are not screening those who had the oppertunity to be vaccinated as a child
                        skip_screening = True
                    else:
                        #In this case, the individual was in a cohort offered vaccines, and we are beyond the switch year, so rescale by m (for cohorts offered vaccine)
                        screen_interval_scale_factor = m 
                else:
                    #In this case we are at the switch year, but the individual was not in a cohort offered a vaccine (was aged 12 too early, before 2008), so rescale by l
                    screen_interval_scale_factor = l 
            else:
                #If we have not yet reached the year at which we switch, we do standard NHS 2025 screening algoorithm; multiply all screening intervals by 1 to keep them the same
                screen_interval_scale_factor = 1


            if 50<=age<=65 and not skip_screening: #age filtering (under 50s get dealt with by the other regular screening intervention) AND only set for_screening to True for this agent if we are not to skip screening for them
                screen_interval = 5*screen_interval_scale_factor

                if sim.t>=date_screened+screen_interval/sim['dt']:
                    for_screening[i] = True
            
    return for_screening




def get_routine_screening_interventions(l, m, 
                                        end_younger_screening_at_switch_year = False):
    routine_screening_under50  = hpv.routine_screening(eligibility=lambda sim: routine_screen_eligible_under50(sim, l, m, end_younger_screening_at_switch_year=end_younger_screening_at_switch_year), 
                                    start_year=GlobalScreeningParameters.screening_start_year, 
                                    age_range=[5,150],#overriding the default age range of [30,50] with a very large one - so as not to interfere with the age range I define in my eligibility checking
                                    prob=GlobalScreeningParameters.primary_screen_prob_under50, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='routine_screening_under50') 
    
    routine_screening_50andover  = hpv.routine_screening(eligibility=lambda sim: routine_screen_eligible_50andover(sim, l, m, end_younger_screening_at_switch_year=end_younger_screening_at_switch_year), 
                                    start_year=GlobalScreeningParameters.screening_start_year, 
                                    age_range=[5,150],#overriding the default age range of [30,50] with a very large one - so as not to interfere with the age range I define in my eligibility checking
                                    prob=GlobalScreeningParameters.primary_screen_prob_50andover, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='routine_screening_50andover') 
    return [routine_screening_under50, routine_screening_50andover]



# First (Consecutive) Cyotology for women with positive routine screening
    #[We will refer to actions as 'consecutive' if they can all be traced back to the same initial routine screening]
to_first_cytology = lambda sim: list(set(
    sim.get_intervention('routine_screening_under50').outcomes['positive'].tolist() +
    sim.get_intervention('routine_screening_50andover').outcomes['positive'].tolist()
))
first_cytology = hpv.routine_triage(eligibility=to_first_cytology,
                                    prob = GlobalScreeningParameters.triage_screen_prob,
                                    age_range=[5,150], #I think routine traige doesnt introduce a constricted age range in by default, but in case it does, overriding it here
                                    product='lbc',
                                    annual_prob=False,
                                    label="first_cytology")

#If first cytology negative, send to a screening in a year (i.e. a second consecutive screening)
def second_screening_eligible(sim):
    for_screening = np.array([False,]*sim.n)
    for i in sim.people.needs_consec_screening_2.keys():
        if sim.people.needs_consec_screening_2[i] is not None and sim.people.needs_consec_screening_2[i] <= sim.t <= sim.people.needs_consec_screening_2[i] + GlobalScreeningParameters.abandon_followup_invites_threshold:
            for_screening[i]=True 
    return for_screening

second_consecutive_screening  = hpv.routine_screening(eligibility=second_screening_eligible, 
                                    start_year=GlobalScreeningParameters.screening_start_year,                                     
                                    age_range=[5,150],#overriding the default age range of [30,50] with a very large one - so as not to interfere with the age range I define in my eligibility checking
                                    prob=GlobalScreeningParameters.secondary_screen_prob, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='second_consecutive_screening') 

#-ve from second screening are returned to standard screening - their date_screened is updated from the second screening, so nothing needed for them
#+ve from second screening do second cytology
to_second_cytology = lambda sim: sim.get_intervention('second_consecutive_screening').outcomes['positive']
second_cytology = hpv.routine_triage(eligibility=to_second_cytology,
                                    prob = GlobalScreeningParameters.triage_screen_prob,
                                    age_range=[5,150], #I think routine traige doesnt introduce a constricted age range in by default, but in case it does, overriding it here
                                    product='lbc',
                                    annual_prob=False,
                                    label="second_cytology")

#-ve from second cytology do a third consecutive HPV screening in 12 months
def third_screening_eligible(sim):
    for_screening = np.array([False,]*sim.n)
    for i in sim.people.needs_consec_screening_3.keys():
        if sim.people.needs_consec_screening_3[i] is not None and sim.people.needs_consec_screening_3[i] <= sim.t <= sim.people.needs_consec_screening_3[i] + GlobalScreeningParameters.abandon_followup_invites_threshold:
            for_screening[i]=True
    return for_screening

third_consecutive_screening  = hpv.routine_screening(eligibility=third_screening_eligible, 
                                    start_year=GlobalScreeningParameters.screening_start_year, 
                                    age_range=[5,150],#overriding the default age range of [30,50] with a very large one - so as not to interfere with the age range I define in my eligibility checking
                                    prob=GlobalScreeningParameters.third_screen_prob, 
                                    product='hpv',        #Screening: DNA HPV Testing
                                    label='third_consecutive_screening') 

#-ve from third screening are returned to standard screening. As their latest screening time is updated in {date_screened}, no need to do anything here
#+ve from third screening perform a third cyotology followed by a mandatory colposcopy
to_third_cytology = lambda sim: sim.get_intervention('third_consecutive_screening').outcomes['positive']
third_cytology = hpv.routine_triage(eligibility=to_third_cytology,
                                    prob = GlobalScreeningParameters.triage_screen_prob,
                                    age_range=[5,150], #I think routine traige doesnt introduce a constricted age range in by default, but in case it does, overriding it here
                                    product='lbc',
                                    annual_prob=False,
                                    label="third_cytology")

#Eligible for a colposcopy if any of the following hold: A)tested not normal for any cytology; B) went all the way to a third screening/cytology
colposcopy_eligible = lambda sim:list(set(sim.get_intervention('first_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('first_cytology').outcomes['ascus'].tolist() +
                                        sim.get_intervention('second_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('second_cytology').outcomes['ascus'].tolist() +
                                        sim.get_intervention('third_cytology').outcomes['normal'].tolist() + sim.get_intervention('third_cytology').outcomes['abnormal'].tolist() + sim.get_intervention('third_cytology').outcomes['ascus'].tolist()))


colposcopy =  hpv.routine_triage(eligibility=colposcopy_eligible, 
                                prob=GlobalScreeningParameters.colpo_prob, 
                                age_range=[5,150], #I think routine traige doesnt introduce a constricted age range in by default, but in case it does, overriding it here
                                product='colposcopy', 
                                annual_prob=False, 
                                label='colposcopy')

hsil_by_colpo = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] #i.e. covers all women diagnosed with HSIL by a colposcopy
ablation = treat_num_canceragnostic(eligibility=hsil_by_colpo, 
                        prob=GlobalScreeningParameters.ablate_prob, 
                        #No need to add an age range here, as by default subclasses of BaseTreatment have an age range of [0,99] - that is, they don't restrict the age range beyond any calculations in the provided eligibility function 
                        product='ablation', 
                        label='ablation') #TODO: can model a max capacity with treat_num too, but as it stands we pretend infinite capacity

cancer_by_colpo = lambda sim: sim.get_intervention('colposcopy').outcomes['cancer'] 
general_cancer_treatment = GeneralCancerTreatment(p_clear=GlobalScreeningParameters.cancer_treatment_effectiveness, p_extend=0.00, extend_dur=dict(dist='normal', par1=18.0, par2=2.))
gct = treat_num_cancer(eligibility=cancer_by_colpo, 
                        prob=GlobalScreeningParameters.generalcancertreatment_prob, 
                        #No need to add an age range here, as by default subclasses of BaseTreatment have an age range of [0,99] - that is, they don't restrict the age range beyond any calculations in the provided eligibility function
                        product=general_cancer_treatment, 
                        label='general_cancer_treatment') 






#TODO: neaten the below! and in checking all of this ofc I will need to clear it up adn probs make it all neater with a big restructure anyways

def get_interventions(l, m, 
                      end_younger_screening_at_switch_year=False):
    #l = screening interval for cohorts not offered vaccine, m=screening interval for cohorts offered vaccine
    interventions = [init_intervention_trackers] + get_routine_screening_interventions(l,m,end_younger_screening_at_switch_year) + [
                            first_cytology, second_consecutive_screening, 
                            second_cytology, third_consecutive_screening, 
                            third_cytology, 
                            colposcopy, 
                            ablation, gct, 
                        update_last_hpv_result, update_needs_consec_screening_2,  update_needs_consec_screening_3,
                                       ]
    return interventions



if __name__=="__main__":
    #HPVsim simulation parameters
    pars = dict(
            n_agents      = 5e3,#20e3,       
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
    sim_interventions =  hpv.Sim(pars, interventions = get_interventions(1,1)
                                 , #interventions can just be functions of a sim
                                            label='NHS Cervical Screening Pathway 2025')


    
    #Run and plot sims
    reduced_sims = []
    for sim in [sim_basic, sim_interventions]:
        msim = hpv.MultiSim(sim, n_runs=10)  
        msim.run()         
        reduced_sims.append(msim.reduce(use_mean=True, output=True))  


    msim = hpv.MultiSim(reduced_sims)
    msim.plot()