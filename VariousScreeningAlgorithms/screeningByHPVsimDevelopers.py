'''
Construct the 7 screen and treat algorithms recommended by the WHO
See documentation here: https://www.ncbi.nlm.nih.gov/books/NBK572308/

Looks here in the HPVsim tutorial like these are "for screening and treatment of cervical pre-cancer lesions" 
    ^ does this mean, that this does NOT actually deal with diagnosis and treatment of the cancer itself? If so, I need to add something in for that for later stage developement, surely?
'''

import hpvsim as hpv
import numpy as np

debug = 1

def make_sim(seed=0, interventions = [], label="MySim"):
    ''' Make a single sim '''

    # Parameters
    pars = dict(
        n_agents        = 50e3,
        dt              = 0.25,
        start           = 1960,
        end             = 2080,
        #ms_agent_ratio  = 10,
        #burnin          = 45,
        rand_seed       = seed,
    )
    sim = hpv.Sim(pars=pars, interventions=interventions, label=label)
    print(sim)
    return sim


def make_algorithms(sim=None, seed=0, debug=debug):


    # Shared parameters
    primary_screen_prob = 0.2
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    start_year = 2025

    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 5 / sim['dt']))
        #in a given timepoint, a woman is eligible for screening iff (they have never been screened before OR their last screening was more than 5 years ago)


    ####################################################################
    #### Algorithm 7 (https://www.ncbi.nlm.nih.gov/books/NBK572308/)
    # HPV DNA as the primary screening test, followed by cytology triage,
    # followed by colposcopy and treatment
    ####################################################################


    """
    Algorithm 7:
    HPV DNA TESTING -> negative => rescreen with HPV test in 5-10 years for general women and 3-5 for those with HIV
                    -> positive => CYTOLOGY TRIAGE -> negative => REPEAT HPV TEST in 3 years for general women population and 1 year for those living with HIV -> negative => back to start (rescreen in 5-10 years normal, 3-5 if hiv)
                                                                                                                                                               -> positive => COLPOSCOPY
                                                   -> ascus or worse => COLPOSCOPY

    ...and after a colposcopy, there is further management based on colposcopy diagnosis or histopathology diagnosis
    """
    #This is not quite what is happening with algorithm 7 though, it is a simpliciation of this processsss!



    hpv_primary7 = hpv.routine_screening(
        product='hpv', #hpv product is a very effective HPV DNA test
        prob=primary_screen_prob,
        eligibility=screen_eligible,
        start_year=start_year,
        label='hpv primary',
    )

    # Send HPV+ women for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv primary').outcomes['positive'] 
    cytology7 = hpv.routine_triage(
        product='lbc',
        annual_prob=False,
        prob=triage_screen_prob,
        eligibility=to_cytology,
        label='cytology',
    )

    # Send ASCUS and abnormal cytology results for colpo
    to_colpo = lambda sim: list(set(sim.get_intervention('cytology').outcomes['abnormal'].tolist() + sim.get_intervention('cytology').outcomes['ascus'].tolist()))
    colpo7 = hpv.routine_triage(
        product='colposcopy',
        annual_prob=False,
        prob=triage_screen_prob,
        eligibility=to_colpo,
        label='colpo',
    )

    # After colpo, treat HSILs with ablation
    hsils = lambda sim: sim.get_intervention('colpo').outcomes['hsil']
    ablation7 = hpv.treat_num(
        prob = ablate_prob,
        product = 'ablation',
        eligibility = hsils,
        label = 'ablation'
    )

    algo7 = [hpv_primary7, cytology7, colpo7, ablation7]
    #for intv in algo7: intv.do_plot=False


    ####################################################################
    #### Set up scenarios to compare algoriths
    ####################################################################

    # Create, run, and plot the simulations
    sim0 = make_sim(label='No screening')
    sim7 = make_sim(interventions=algo7, label='Algorithm 7')
    msim = hpv.parallel([sim0, sim7])
    #sim7.run()
    #sim7.plot()


    msim.compare()

    msim.plot()

    return msim




#%% Run as a script
if __name__ == '__main__':

    msim = make_algorithms()
