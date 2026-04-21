"""
This script evaluates an ensemble-calibration, by computing GOF against seen (and then unseen) real-world observed data

Either run [2019,2020,2021,2022]  as the years, with the START datafile, or   [2012,2013,2014,2015,2016,2017,2018] as the years with the END datafile to evaluate GOF against seen or unseen data.

The genotype distribution for cancerous cases is hardcoded from aquired data, but can also be changed.

To run this code, we need to change the name of the dummy calibration with which we set up a sim object to then have its parameters changed, to a file that does not yet exist.
"""

## Imports
import matplotlib.pyplot as plt

import numpy as np

import sciris as sc

import pandas as pd

from tqdm import tqdm

import pickle

from InterventionAlgorithms import screeningAlgorithms, NHS_Vacc
from basePars import base_pars
import hpvsim as hpv

from scenario_configuration import cal_filename  #import the name of this calibration according to the parameterisation of the model we are calibrating to


#Whether to plot for start years (2012-2018 + 2011 cancerous genotype data) or end years (2019-2022)
plot_start = True 


#Name for the dummy calibration
cal_name = "project_modelling/calibration_results/dummy_21Apr26_1"

#Simulation start and end dates
start = 1980
end = 2025
base_pars['start'] = start
base_pars['end'] = end

#hardcoded data for cancerous genotype distribution, 2011
cancerous_genotype_dist_year=2011
cancerous_genotype_dist_TRUEDATA = [0.677,0.188,0.163,0.059] #hpv16, hpv18, hi5, ohr

#Number of repeat runs for each parameter set
repeat_runs = 4

#Load in final calibration
with open(cal_filename, 'rb') as file:
    loadeddata = pickle.load(file)

final_cal_data = loadeddata['final_cal_data'][:2] #this is a list of parameter-tuples, so has a fixed ordering
par_labels = loadeddata['par_labels']



if __name__=="__main__":   
    #Configure sim template for use in the rest of this script
    pars = base_pars
    pars['interventions'] =  screeningAlgorithms.get_interventions(v=5,u=5,a=5)  + NHS_Vacc.vaccinations
    sim = hpv.Sim(pars) 
    
    #Set up the dummy calibration - to end up having an object with which we can easily change the parameters of a HPVsim instance
    calib_pars = dict( #to make the dummy calibration work, doing a cal with just beta. to make sure this script works as desired, make sure the dummy calibration is for a varaible which will be overriden with our hardcoded parameters - picking beta here as i think any calibration I do will involve beta, in which case it will be always overriden
            beta=[0.05, 0.00, 0.20], 
        )
    calib = hpv.Calibration(
        sim,
        calib_pars=calib_pars,
        datafiles= ["project_modelling/real_world_data/mesherHPVCancerDist.csv"],
        total_trials=1,
        n_workers=1,
        keep_db=True,
        name=cal_name
    )
    calib.calibrate(die=False)
    sim_pars = calib.trial_pars_to_sim_pars() # Returns best parameters from calibration in a format ready for updating other sims' pars


    #Set up analyzer for recording sim results corresponding to accumulated cancer-by-age-and-year data
    az_cancers = hpv.age_results(
        result_args=sc.objdict(
            cancers=sc.objdict(
                datafile="project_modelling/real_world_data/new_cervical_cancer_cases_ENGSCALED1P19TOUK_START.csv" if plot_start else "project_modelling/real_world_data/new_cervical_cancer_cases_ENGSCALED1P19TOUK_END.csv"
            ),
        )
    )
    years =  [2012,2013,2014,2015,2016,2017,2018] if plot_start else [2019,2020,2021,2022]     #Corresponds to above files


    

    #Define accumulating variables - each simulation run adds an element to each of these lists
    all_cancer_results = [] #stores analyser results for cancers by age by year
    gofs_cancers_normalised = [] #storing goodness of fits against number of cancers by year and age bin
    total_gofs_normalised = [] #storing overall goodness of fits (i.e. if {plot_start}, adds up GOF against both cancers and cancerous genotype distribution, else this is the same array as {gof_cancers_normalised})

    if plot_start:
        #genotype distribution
        all_cancerous_genotype_results_16 = [] 
        all_cancerous_genotype_results_18 = [] 
        all_cancerous_genotype_results_hi5 = [] 
        all_cancerous_genotype_results_ohr = [] 
        gofs_cancerous_genotypes_normalised = [] #storing goodness of fits against cancerous genotype distribution


    #Iterate over each parameterisation in our final calibration
    for i in tqdm(range(len(final_cal_data))):
        print(f"--- Generating results for parameterisation {i} of {len(final_cal_data)} ---")

        #Update sim parameters with ith parametersation in final calibration
        cal_params = final_cal_data[i][0]

        sim_pars['genotype_pars']['hpv16']['cin_fn']['k'] = cal_params[par_labels.index('params_hpv16_cin_fn_k')]
        sim_pars['genotype_pars']['hpv16']['dur_cin']['par1'] = cal_params[par_labels.index('params_hpv16_dur_cin_par1')]

        sim_pars['genotype_pars']['hpv18']['cin_fn']['k'] = cal_params[par_labels.index('params_hpv18_cin_fn_k')]
        sim_pars['genotype_pars']['hpv18']['dur_cin']['par1'] = cal_params[par_labels.index('params_hpv18_dur_cin_par1')]

        sim_pars['genotype_pars']['hi5']['cin_fn']['k'] = cal_params[par_labels.index('params_hi5_cin_fn_k')]
        sim_pars['genotype_pars']['hi5']['dur_cin']['par1'] = cal_params[par_labels.index('params_hi5_dur_cin_par1')]
        sim_pars['genotype_pars']['hi5']['rel_beta'] = cal_params[par_labels.index('params_hi5_rel_beta')]

        sim_pars['genotype_pars']['ohr']['cin_fn']['k'] = cal_params[par_labels.index('params_ohr_cin_fn_k')]
        sim_pars['genotype_pars']['ohr']['dur_cin']['par1'] = cal_params[par_labels.index('params_ohr_dur_cin_par1')]
        sim_pars['genotype_pars']['ohr']['rel_beta'] = cal_params[par_labels.index('params_ohr_rel_beta')]

        sim_pars['beta'] = cal_params[par_labels.index('params_beta')]
        sim_pars['f_cross_layer'] = cal_params[par_labels.index('params_f_cross_layer')]
        sim_pars['m_cross_layer'] = cal_params[par_labels.index('params_m_cross_layer')]

        sim = hpv.Sim(base_pars, analyzers=[az_cancers])
        sim.update_pars(sim_pars)

        #Run sim as a multisim, with required number of repeats
        msim = hpv.MultiSim(sim, n_runs=repeat_runs)
        msim.run()

        #Extract results for each sim in the multisim
        for msim_sim in msim.sims:
            an = msim_sim.get_analyzer() #If I am to add mutliple analysers, this would need to change slightly
            all_cancer_results.append(an.results)

            an.result_args['cancers'].weights=np.ones(an.result_args['cancers'].data['value'].shape)
            this_gof_cancers_normalised = an.compute_mismatch('cancers')/(len(years)*17)
            gofs_cancers_normalised.append(this_gof_cancers_normalised)

            if plot_start:
                a,b,c,d = msim_sim.results['cancerous_genotype_dist'][:,(cancerous_genotype_dist_year-start)]
                all_cancerous_genotype_results_16.append(a)
                all_cancerous_genotype_results_18.append(b)
                all_cancerous_genotype_results_hi5.append(c)
                all_cancerous_genotype_results_ohr.append(d)

                gof = hpv.misc.compute_gof(cancerous_genotype_dist_TRUEDATA, [a,b,c,d])
                this_gof_cancerous_genotypes_normalised = gof.sum()/4
                gofs_cancerous_genotypes_normalised.append(this_gof_cancerous_genotypes_normalised)

                total_gofs_normalised.append(this_gof_cancers_normalised+this_gof_cancerous_genotypes_normalised)
            else:
                total_gofs_normalised.append(this_gof_cancers_normalised)
        
    all_cancer_results = np.array(all_cancer_results) #Converts cancer analyser results to array shape (n_runs, n_timepoints) (n_runs=len(param_list)*repeat_runs)
    

    print("Completed result generation. Now plotting...")

    #Set up axes on which things are to be plotted
    
    if plot_start:
        figGeno,axGeno = plt.subplots(1,1) #If 2011 is included in the years we are evaluating against, we also evaluate against our hardcoded genotype distribution
    
    figGOF_normalised, axGOF_normalised = plt.subplots(1,1) #overall normalised GOF (i.e. divided by number of datapoints in the calculation to get approximate error per datapoint)

            
    #Plot cancer results
    axsCanc = [] #to store yearly plots of cancers by age
    for _ in range(len(years)):
        _,a = plt.subplots(1,1)
        axsCanc.append(a)

    cancers_bins = [ 0., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90.] #17 elements here
    cancers = {} #cancers[year][bin_index] = [list of values across parametersations and repeat runs for bin cancers_bins[bin_index]]

    for year in years:
        year_data = []
        for bin_index in range(len(cancers_bins)):
            bin_data = []
            for r in range(all_cancer_results.shape[0]):
                bin_data.append(all_cancer_results[r]['cancers'][year][bin_index])
            year_data.append(bin_data)
        cancers[year] = year_data
    
    #plot by year
    box = 0 #starting box, offset for blank space at the start
    for year in years:
        ax = axsCanc[box]; box+=1

        cancer_year_data = cancers[year]
        bp = ax.boxplot(cancer_year_data, tick_labels=cancers_bins if i==1 else [' ']*len(cancers_bins),
                                            patch_artist=True,
            positions = np.arange(1,1+len(cancers_bins)))
        
        #Overlay plot from datafile
        datafilesource = f"project_modelling/real_world_data/new_cervical_cancer_cases_ENGSCALED1P19TOUK_{year}.csv" 
        df = pd.read_csv(datafilesource)
        true_data = df['value']
        true_data.tolist()
        for j in range(1, len(true_data)+1):
            ax.hlines(true_data[j-1],j-0.5, j+0.5, color='red' )

        ax.set_title(f"Cancers by age bracket, {year}")
            


    #Plot cancerous genotypes results
    if plot_start:
        axGenoPlot = axGeno.boxplot([all_cancerous_genotype_results_16, all_cancerous_genotype_results_18, all_cancerous_genotype_results_hi5, all_cancerous_genotype_results_ohr], 
                                                patch_artist=True,
                    tick_labels=["HPV 16", "HPV 18", "HI5", "OHR"] if i==1 else [' ',' ',' ',' '],
                positions = np.array([1,2,3,4]))
        for j in range(1,5):
            axGeno.hlines(cancerous_genotype_dist_TRUEDATA[j-1],j-0.5, j+0.5, color='red' )

        axGeno.set_title(f"Cancerous Genotype Distribution, {cancerous_genotype_dist_year}")
        


    #--- Plot normalised GOF data ---#
    if plot_start:
        axGOFplot = axGOF_normalised.boxplot([gofs_cancers_normalised, gofs_cancerous_genotypes_normalised], 
                                                patch_artist=True,
                                            tick_labels=["cancers", "cancerous genotype distribution"],
                                            positions = np.array([1,2]))
    else:
        axGOFplot = axGOF_normalised.boxplot([gofs_cancers_normalised], 
                                                patch_artist=True,
                                            tick_labels=["cancers"],
                                            positions = np.array([1]))
    axGOF_normalised.set_title("GOF (relative to data sizes)")
    axGOF_normalised.set_ylim(0,0.25)

    print("Info on total GOF:")
    print(f"\t min:{min(total_gofs_normalised)} \n\tLQ:{np.percentile(total_gofs_normalised, 25)} \n\tmedian:{np.median(total_gofs_normalised)} \n\tUQ:{np.percentile(total_gofs_normalised, 75)} \n\t max:{max(total_gofs_normalised)}")
    print(f"(recall, this is for final pooled calibration {cal_filename})")





    plt.show()
