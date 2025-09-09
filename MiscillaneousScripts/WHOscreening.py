#The 7 algorithms recommended in the WHO's guidelines for screening and treatment of cervical pre-cancer lesions.

import hpvsim as hpv
import numpy as np

debug = 1

def make_sim(seed=0):
    ''' Make a single sim '''
seed = 0
# Parameters
pars = dict(
    n_agents        = [50e3,5e3][debug],
    dt              = [0.5,1.0][debug],
    start           = [1975,2000][debug],
    end             = 2060,
    ms_agent_ratio  = 10,
    burnin          = [45,0][debug],
    rand_seed       = seed,
    )
sim = hpv.Sim(pars=pars)
    #return sim


def make_algorithms(sim=None, seed=0, debug=debug):

    if sim is None: sim = make_sim(seed=seed)

    # Shared parameters
    primary_screen_prob = 0.2
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    start_year = 2025
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 5 / sim['dt']))


    ####################################################################
    #### Algorithm 1 (https://www.ncbi.nlm.nih.gov/books/NBK572308/)
    # Visual inspection with acetic acid (VIA) as the primary screening test, followed by treatment
    ####################################################################

primary_screen_prob = 0.6
screen_eligible = 25
start_year = 2015

via_primary = hpv.routine_screening(
    product='via',
    prob=primary_screen_prob,
    eligibility=screen_eligible,
    start_year=start_year,
    label='via primary',
    )

ablate_prob = 0.8

via_positive = lambda sim: sim.get_intervention('via primary').outcomes['positive']
ablation1 = hpv.treat_num(
    prob = ablate_prob,
    product = 'ablation',
    eligibility = via_positive,
    label = 'ablation'
    )

algo1 = [via_primary, ablation1]
for intv in algo1: intv.do_plot=False

    ####################################################################
    #### Set up scenarios to compare algoriths
    ####################################################################

    # Create, run, and plot the simulations
sim0 = hpv.Sim(label='No screening')
sim1 = hpv.Sim(interventions=algo1, label='Algorithm 1')

msim = hpv.parallel([sim0, sim1])
msim.compare()

#return msim

#%% Run as a script
if __name__ == '__main__':

    msim = make_algorithms()