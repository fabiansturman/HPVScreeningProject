#Imports
import pickle
from resultGeneration import run_sim_and_save_raw_result
import project_modelling.scenario_configuration as scenario_configuration


#Load in properties of this round of result-generation from scenario_parameterisation.py
seeds = scenario_configuration.seeds
alg_name = scenario_configuration.alg_name
pars = scenario_configuration.adapted_pars
cal_filename = scenario_configuration.cal_filename

print(f"Loading calibrated parameters from {cal_filename}")

if __name__=="__main__":   
    #Load in calibrated parameter sets
    with open(cal_filename, 'rb') as file:
        loadeddata = pickle.load(file)
    final_cal_data = loadeddata['final_cal_data'] #this is a list of parameter-tuples, so has a fixed ordering
    par_labels = loadeddata['par_labels']

    #Generate results for each seed for the algorithm defined in this script
    for seed in seeds:
        run_sim_and_save_raw_result(pars, seed,
                                    final_cal_data, par_labels,
                                    alg_name)


        
        