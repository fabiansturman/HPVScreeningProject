
# Imports
import copy
import pandas as pd
import numpy as np
from basePars import base_pars
from InterventionAlgorithms import screeningAlgorithms, NHS_Vacc, GlobalScreeningParameters


##THE FOLLOWING MAY BE CHANGED TO SET DIFFERENT INTERVENTIONS##
screening_algorithm_stub_name = 'V10U5A5' #set the name we will be using for screening algorithms
screening_algorithm = screeningAlgorithms.get_interventions(v=10,u=5,a=5)
"""                                                         ^
                                                            v = screening interval for vaccinated individuals
                                                            u = screening interval for unvaccinated vaccine-elgiible individuals
                                                            a = screening interval for non-vaccine-eligible individuals
"""


vaccination_algorithm = NHS_Vacc.vaccinations


##THE FOLLOWING MAY BE CHANGED TO SET DIFFERENT SETS OF SEEDS OVER WHICH TO GENERATE FINAL RESULTS##
seeds = [0,1]#,2,3] # Define seeds we are generating the results over (these raw seeds will be offset by parameterisation number, so that we get a much wider set of seeds when iterating over different parameters from the calibration )


##THE FOLLOWING MAY BE CHANGED *WITH CAUTION* TO DEFINE DIFFERENT ENCODING OF SCENARIOS IN FILENAMES##
# Extract parameter values from GlobalScreeningParameters and interventions themselves, to append to naming of result-files to specify their scenarios (for sensitivity analysis)
uptk = GlobalScreeningParameters.projected_teen_vaccination_uptake 
cte = GlobalScreeningParameters.cancer_treatment_effectiveness 
psp1 = GlobalScreeningParameters.primary_screen_prob_under50 
psp2 = GlobalScreeningParameters.primary_screen_prob_50andover  
ssp = GlobalScreeningParameters.secondary_screen_prob 
tsp = GlobalScreeningParameters.third_screen_prob 
ablate_prob = GlobalScreeningParameters.ablate_prob 
gtp = GlobalScreeningParameters.generalcancertreatment_prob 

dx = pd.read_csv('hpvsim/data/products_dx.csv')
tx = pd.read_csv('hpvsim/data/products_tx.csv')
vx = pd.read_csv('hpvsim/data/products_vx.csv')
ohr_rel_imm = np.array(vx[(vx["name"] == "nonavalent") & 
                          (vx["genotype"] == "ohr")]['rel_imm']).item() 
cytology_cin_sensitivity = np.array(dx[(dx["name"]=="lbc") &
                                       (dx["state"]=="cin") &
                                       (dx["result"]=="abnormal")]["probability"]).item()
cytology_specificity = np.array(dx[(dx["name"]=="lbc") &
                                       (dx["state"]=="susceptible") &
                                       (dx["result"]=="normal")]["probability"]).item()

sum_male_init_prev = np.sum(base_pars['init_hpv_prev']['m']) #for the cases we are dealing with, it is sufficient to add up all the numbers in our initial male HPV distribution to uniquely identify which of the two distributions we are dealing with


def twosig(value):
    #We specify all these numerical scenario parameters by the first two significant figures as 2 digits on their own - this is sufficient to tell apart all cases we are considering in our sensitivity analysis
    s = str(value).replace('.', '').replace('-', '') #convert to string, remove decimal point (and any minus sign, though this is irrelevant for this, included for completeness of this function)
    s = s.lstrip('0') #remove all leading zeros
    return s[:2] #print two significant digits (or 1 if there is only 1 digit left)
    

    #{calibration_code} uniquely identifies the calibrated set of parameters to use for this sensitivity analysis
calibration_code = f"_{twosig(cte)}_{twosig(psp1)}_{twosig(psp2)}_{twosig(ssp)}_{twosig(tsp)}_{twosig(ablate_prob)}_{twosig(gtp)}_{twosig(ohr_rel_imm)}_{twosig(cytology_specificity)}_{twosig(cytology_cin_sensitivity)}_{twosig(sum_male_init_prev)}"
    #{scenario_code} uniquely identifies the full scenario being modelled (i.e. the modelling parameters relevant pre-2026 that needs a calibratiob for it, and the assumed later vaccine uptake)
scenario_code = f"_{twosig(uptk)}{calibration_code}"


##DO NOT CHANGE CODE BEYOND THIS POINT##


# Define the algorithm we are using in this script and add to basepars
alg_name = screening_algorithm_stub_name + scenario_code

adapted_pars = copy.deepcopy(base_pars)
adapted_pars['interventions'] =  screening_algorithm + vaccination_algorithm

# Define the name of the calibration relevant for this
cal_filename = f"project_modelling/calibration_results/FinalPooledCalibration{calibration_code}.pickle"



        
        