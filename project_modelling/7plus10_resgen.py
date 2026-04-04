#Imports
import copy
import pickle
from basePars import base_pars
from InterventionAlgorithms import NHS_2025_lambdamuForVaccCohortsOnly, NHS_Vacc
from resultGeneration import run_sim_and_save_raw_result

#Define the seeds we are using
seeds = [0,1,2,3,4]

#Define the algorithm we are using in this script
l=1.4
m=2
alg_name = "7plus10_HV" #need alg name to save results from simulation runs to the correct pickle file

base_pars = copy.deepcopy(base_pars)
base_pars['interventions'] =  NHS_2025_lambdamuForVaccCohortsOnly.get_interventions(l=l, m=m) + NHS_Vacc.vaccinations

if __name__=="__main__":   
    #Load in calibrated parameter sets
    with open('finalCalibration.pickle', 'rb') as file:
        loadeddata = pickle.load(file)
    final_cal_data = loadeddata['final_cal_data'] #this is a list of parameter-tuples, so has a fixed ordering
    par_labels = loadeddata['par_labels']

    #Generate results for each seed for the algorithm defined in this script
    for seed in seeds:
        run_sim_and_save_raw_result(base_pars, seed,
                                    final_cal_data, par_labels,
                                    alg_name)

        
        