#Imports
import pickle
from resultGeneration import initialise_blank_ts_file, run_sim_and_save_raw_result
import scenario_configuration as scenario_configuration

import time

#Load in properties of this round of result-generation from scenario_parameterisation.py
seeds = scenario_configuration.seeds
alg_name = scenario_configuration.alg_name
pars = scenario_configuration.adapted_pars
cal_filename = scenario_configuration.cal_filename


if __name__=="__main__":   
    #Load in calibrated parameter sets
    print(f"Loading calibrated parameters from {cal_filename}")
    with open(cal_filename, 'rb') as file:
        loadeddata = pickle.load(file)
    final_cal_data = loadeddata['final_cal_data'] #this is a list of parameter-tuples, so has a fixed ordering
    par_labels = loadeddata['par_labels']


    #Initialise a blank file to contain the results generated for this algorithm (the function provides the user the option to not override an existing file, if an existing file is detected)
    initialise_blank_ts_file(algs = [alg_name])

    #print("sleeping before getting started, to avoid RAM overload :)")
    #time.sleep( 1.75*60*60 * 5.5 )
    #print("now getting started!")

    #Generate the raw results
    for seed in seeds:
        run_sim_and_save_raw_result(pars, seed,
                                    final_cal_data, par_labels,
                                    alg_name)


        
        