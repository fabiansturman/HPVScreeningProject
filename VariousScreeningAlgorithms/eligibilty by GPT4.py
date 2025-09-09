import hpvsim as hpv
import numpy as np

# ... [Your existing code for defining functions and parameters]
# Define the parameters with an additional attribute for vaccination status
# Define the parameters with an additional attribute for vaccination status
pars = dict(
    n_agents      = int(20e3),  # Ensure this is an integer
    n_years       = 60,         
    verbose       = 0,          
    rand_seed     = 2,          
    genotypes     = [16, 18, 'hr']
)

# Create the simulation instance
orig_sim = hpv.Sim(pars, label='Baseline')

# Initialize the simulation with vaccinated status
vaccinated_num = int(orig_sim.pars['n_agents'] * 0.75)  # Ensure this is an integer

# Create an array of agent indices and then select vaccinated indices
agent_indices = np.arange(orig_sim.pars['n_agents'])  # Array of agent indices
np.random.seed(orig_sim.pars['rand_seed'])  # Set the random seed
vaccinated_indices = np.random.choice(agent_indices, size=vaccinated_num, replace=False)

orig_sim.people.vaccinated = np.zeros(orig_sim.pars['n_agents'], dtype=bool)
orig_sim.people.vaccinated[vaccinated_indices] = True
# Modify the eligibility function to include vaccination status
def screen_eligible_vaccinated(sim):
    return (sim.people.vaccinated) & (np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt'])))

# ... [Modify other eligibility functions similarly]

# Define interventions with the new eligibility function
start_year = 2024
primary_screen_prob = 0.6
triage_screen_prob = 0.9
ablate_prob = 0.9
import numpy as np


# Use screen_eligible_vaccinated instead of screen_eligible in your interventions
hpv_primary7  = hpv.routine_screening(eligibility=screen_eligible_vaccinated, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary7')
# ... [Modify other interventions similarly]

# Create the simulation instance with the updated eligibility
algo = hpv.Sim(pars, interventions=[hpv_primary7, ...], label='NHS Cervical Screening Pathway 2023 (3 Year Call Back)')
# ... [Define other simulations similarly]

# Run and plot
msim = hpv.parallel(orig_sim, algo, ...)
msim.plot(['cancer_incidence']);