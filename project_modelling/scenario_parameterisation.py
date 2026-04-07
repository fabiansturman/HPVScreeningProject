"""
This file defines the parameterisation for a round of result-generation. 
It builds on the parameterisation defined in base_pars, the intervention csv files, and calibrated parameters by adding on specific interventions we are modelling for results in this scenario.

USAGE:
    -> Set the parameters in this file and then run 'result_generation.py' to populate a file of results.
    -> Import this file (so that it is run and {cal_filename} is generated) when a calibration is performed to name the calibration according to the sensitivity-analysis scenario we are calibrating to
"""



# Imports
import copy
import pandas as pd
import numpy as np
from basePars import base_pars
from InterventionAlgorithms import NHS_2025_lambdamuForVaccCohortsOnly, NHS_2025_lambdamu, NHS_Vacc, GlobalScreeningParameters



# Extract parameter values from GlobalScreeningParameters and interventions themselves, to append to naming of result-files to specify their scenarios (for sensitivity analysis)
uptk = GlobalScreeningParameters.projected_teen_vaccination_uptake 
cte = GlobalScreeningParameters.cancer_treatment_effectiveness 
psp1 = GlobalScreeningParameters.primary_screen_prob_under50 
psp2 = GlobalScreeningParameters.primary_screen_prob_50andover 
ssp = GlobalScreeningParameters.secondary_screen_prob 
tsp = GlobalScreeningParameters.third_screen_prob 
colpo_prob = GlobalScreeningParameters.colpo_prob
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

def twosig(value):
    #We specify all these numerical scenario parameters by the first two significant figures as 2 digits on their own - this is sufficient to tell apart all cases we are considering in our sensitivity analysis
    s = str(value).replace('.', '').replace('-', '') #convert to string, remove decimal point (and any minus sign, though this is irrelevant for this, included for completeness of this function)
    s = s.lstrip('0') #remove all leading zeros
    return s[:2] #print two significant digits (or 1 if there is only 1 digit left)
    

    #{calibration_code} uniquely identifies the calibrated set of parameters to use for this sensitivity analysis
calibration_code = f"_{twosig(cte)}_{twosig(psp1)}_{twosig(psp2)}_{twosig(ssp)}_{twosig(tsp)}_{twosig(colpo_prob)}_{twosig(ablate_prob)}_{twosig(gtp)}_{twosig(ohr_rel_imm)}_{twosig(cytology_cin_sensitivity)}"
    #{scenario_code} uniquely identifies the full scenario being modelled (i.e. the modelling parameters relevant pre-2026 that needs a calibratiob for it, and the assumed later vaccine uptake)
scenario_code = f"_{twosig(uptk)}{calibration_code}"



# Define seeds we are generating the results over (these raw seeds will be offset by parameterisation number, so that we get a much wider set of seeds when iterating over different parameters from the calibration )
seeds = [0,1,2,3]

# Define the algorithm we are using in this script and add to basepars
alg_name = "5_15" + scenario_code
screening_interventions = NHS_2025_lambdamuForVaccCohortsOnly.get_interventions(l=1, m=3)

adapted_pars = copy.deepcopy(base_pars)
adapted_pars['interventions'] =  screening_interventions + NHS_Vacc.vaccinations

# Define the name of the calibration relevant for this
cal_filename = f"project_modelling/calibration_results/FinalPooledCalibration{calibration_code}"



        
        