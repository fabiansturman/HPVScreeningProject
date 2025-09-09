import hpvsim as hpv
import matplotlib.pyplot as plt
import numpy as np

# Define a range of probabilities for sensitivity analysis
probabilities = np.linspace(0.1, 0.9, 9)

# Define the parameters
pars = dict(
    n_agents      = 20e3,       # Population size
    n_years       = 50,         # Number of years to simulate
    verbose       = 0,          # Don't print details of the run
    rand_seed     = 2,          # Set a non-default seed
    genotypes     = [16, 18],   # Include the two genotypes of greatest general interest
)

# Create the sim with and without interventions
orig_sim = hpv.Sim(pars, label='Baseline')
orig_sim.run()
print("Baseline Simulation Summary:", orig_sim.summary)

# Run simulations for different probabilities
simulation_results = []

for prob in probabilities:
    vx = hpv.routine_vx(prob=prob, start_year=2015, age_range=[9, 10], sex=['f', 'm'], product='nonavalent')
    sim = hpv.Sim(pars, interventions=vx, label=f'With vaccination, prob={prob:.2f}')
    sim.run()
    simulation_results.append(sim)

print(f"simultation results = {simulation_results}")

# Calculate relative reduction in HPV prevalence
baseline_prevalence = orig_sim.summary.hpv_prevalence
relative_reduction = [[(baseline_prevalence - sim.summary.hpv_prevalence) / baseline_prevalence * 100,]
                       for sim in simulation_results]

print(baseline_prevalence)
print(relative_reduction)
print(probabilities)
print([f'{prob:.2f}' for prob in probabilities])

# Whisker diagram
plt.boxplot(relative_reduction, labels=[f'{prob:.2f}' for prob in probabilities])
plt.xlabel('Vaccination Probability')
plt.ylabel('Relative Reduction in HPV Prevalence (%)')
plt.title('Sensitivity Analysis: Relative Reduction in HPV Prevalence')
plt.show()
