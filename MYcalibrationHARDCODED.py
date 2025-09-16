# Import HPVsim
import hpvsim as hpv
import numpy as np

import numpy as np

import pickle

import os

import sciris as sc

import matplotlib.pyplot as plt

#Name for the dummy calibration
cal_name = "CalibrationRawResults\\d2Cal_16Sep25_132"

#TODO: need to add genotype tracking stuff here tooooooo for the genotype data file

#TODO: analyzers have a 'compute mismatch' function which should allow me to compute mismatch in here! try it out!
    #^ i think because the analysier has a buildin datafile, mismatch may be added to the sim??

    #TODO: confrim both the above are correct wth a run of a recent ca;

#TODO: seperately report GOFs against data that has not been used in the calibration, just by swapping out the cancers analyser for one with later data!

#data for cancerous genotype distribution
cancerous_genotype_dist_year=2011
cancerous_genotype_dist_TRUEDATA = [0.679,0.19,0.161]

#Number of repeat runs for each parameter set
repeat_runs = 80#20

#Define the hardcoded parameters I will be using
sophie_married_matrix = [        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [10, 0, 0, 0.08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [15, 0, 0, 0.08, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [20, 0, 0, 0, 0, 0.6, 2, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0],        [25, 0, 0, 0, 0, 0.6, 1, 2, 0.4, 0.1, 0, 0, 0, 0, 0, 0, 0],        [30, 0, 0, 0, 0, 0.5, 0.5, 2, 1, 0.5, 0.1, 0, 0, 0, 0, 0, 0],        [35, 0, 0, 0, 0, 1, 0.5, 1, 2, 1, 0.5, 0.2, 0, 0, 0, 0, 0],        [40, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0],        [45, 0, 0, 0, 0, 0.1, 1, 2, 2, 2, 1, 0.5, 0.2, 0.08, 0, 0, 0],        [50, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 2, 2, 0.5, 0.2, 0.05, 0, 0],        [55, 0, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 3, 2, 1, 0.3, 0.1, 0.1],        [60, 0, 0, 0, 0, 0, 0, 0.1, 0.5, 1, 2, 3, 3, 2, 0.5, 0.3, 0.1],        [65, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 2, 3, 3, 2, 1, 0.2],        [70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 3, 2, 1],        [75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3],    ]
sophie_married_matrix = np.array(sophie_married_matrix)
sophie_casual_matrix = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],        [5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],        [10,0,0,0,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0],        [15,0,0,1,2,3,2,1,0.5,0,0,0,0,0,0,0,0],        [20,0,0,0.15,2,3,2,2,1,0.15,0,0,0,0,0,0,0],        [25,0,0,0.15,0.25,1,2,2,1,1,0,0,0,0,0,0,0],        [30,0,0,0,0,0.5,0.5,2,1,0.15,0,0,0,0,0,0,0],        [35,0,0,0,0,1,0.5,1,2,1,0.5,0,0,0,0,0,0],        [40,0,0,0,0,1,1,1,1,1,0.5,0.25,0,0,0,0,0],        [45,0,0,0,0,0.15,1,2,2,2,1,0.5,0.2,0.1,0,0,0],        [50,0,0,0,0,0,0.15,1,2,3,2,2,0.5,0.2,0.05,0,0],        [55,0,0,0,0,0,0,0.15,1,2,3,3,2,1,0.25,0.1,0.1],        [60,0,0,0,0,0,0,0.15,0.15,1,2,3,3,2,0.5,0.25,0.1],        [65,0,0,0,0,0,0,0,0,0,1,1,2,2,1,0.5,0],        [70,0,0,0,0,0,0,0,0,0,0,0,0,0.8,1,0.7,0.5],        [75,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.25]    ]
sophie_casual_matrix = np.array(sophie_casual_matrix)

start = 1980; end =2017
base_pars = dict(n_agents= 100e3,
                start=start, end=end, dt=0.25, 
                location='united kingdom', 
                verbose=-1,
                debut=dict(f=dict(dist='normal', par1=16.0, par2=3.1), m=dict(dist='normal', par1=16.0, par2=4.1)),
                mixing = {'m':sophie_married_matrix,
                          'c':sophie_casual_matrix},
                condoms = dict(m=0.17, c=0.50), #condom usage in (m)arried and (c)asual relationships
                genotypes     = ['hpv16', 'hpv18', 'hi5'],
                init_hpv_prev = {
                    'age_brackets'  : np.array([  12,   25,   34,   44,  54,   64, 150]),
                    'm'             : np.array([ 0.0, 0.25, 0.15,   0.075,  0.06,   0.06, 0.03]),
                    'f'             : np.array([ 0.0,0.25, 0.15,   0.075,  0.06,   0.06, 0.03])
                },
                init_hpv_dist = {
                    'hpv16': 0.023,
                    'hpv18': 0.009,
                    'hi5':  0.022 #HPV 33 is not listed as one of the top 10 most prevalent in general population in (), so we can assume its prevalence is at most 0.004 - so not adding this to the sum
                })

param_list = []
#d2Cal_15Sep25_XPS_EB_5 pars (best 5)
param_list.append({ #best, GOF=18.73
    "beta":0.07494,
    "f_cross_layer":0.0642,
    "hi5_cin_fn_k":0.1154,
    "hi5_dur_cin_par1":3.117,
    "hi5_rel_beta":0.7839,
    "hpv16_cin_fn_k":0.001433,
    "hpv16_dur_cin_par1":4.5429,
    "hpv18_cin_fn_k":0.32812,
    "hpv18_dur_cin_par1":2.8408,
    "m_cross_layer":0.07521,
})
"""
param_list.append({ #2nd best, GOF=19.54
    "beta":0.34148,
    "f_cross_layer":0.1116,
    "hi5_cin_fn_k":0.04146,
    "hi5_dur_cin_par1":4.9345,
    "hi5_rel_beta":0.6599,
    "hpv16_cin_fn_k":0.10139,
    "hpv16_dur_cin_par1":2.2458,
    "hpv18_cin_fn_k":0.3061,
    "hpv18_dur_cin_par1":1.8204,
    "m_cross_layer":0.58533,
})
param_list.append({ #3rd best, GOF=19.98
    "beta":0.1954,
    "f_cross_layer":0.1405,
    "hi5_cin_fn_k":0.01748,
    "hi5_dur_cin_par1":4.8042,
    "hi5_rel_beta":0.6885,
    "hpv16_cin_fn_k":0.092008,
    "hpv16_dur_cin_par1":2.8175,
    "hpv18_cin_fn_k":0.23367,
    "hpv18_dur_cin_par1":1.8061,
    "m_cross_layer":0.4554,
})
param_list.append({ #4th best, GOF=20.09
    "beta":0.2347,
    "f_cross_layer":0.1922,
    "hi5_cin_fn_k":0.01974,
    "hi5_dur_cin_par1":3.5182,
    "hi5_rel_beta":0.7973,
    "hpv16_cin_fn_k":0.05673,
    "hpv16_dur_cin_par1":3.4901,
    "hpv18_cin_fn_k":0.25921,
    "hpv18_dur_cin_par1":1.9170,
    "m_cross_layer":0.21027,
})
"""

"""
#d2Cal_8Sep25_6 pars (best 4)
param_list.append({"hpv16_cin_fn_k": 0.246296,  # No1 pars from d2Cal_8Sep25_6, fit 20.19 
                "hpv16_dur_cin_par1": 5.19969,
                "hpv18_cin_fn_k":0.268885,
                "hpv18_dur_cin_par1":4.00459,
                "beta":0.04204,
                "f_cross_layer":0.000168,
                "m_cross_layer":0.001859,
            })
param_list.append({"hpv16_cin_fn_k": 0.21203,  # No2 pars from d2Cal_8Sep25_6, fit 20.77
                "hpv16_dur_cin_par1":5.3747,
                "hpv18_cin_fn_k":0.20886,
                "hpv18_dur_cin_par1":5.52436,
                "beta":0.06974,
                "f_cross_layer":0.00007572,
                "m_cross_layer":0.0005274,
            })
param_list.append({"hpv16_cin_fn_k": 0.21586,  # No3 pars from d2Cal_8Sep25_6, fit 21.25
                "hpv16_dur_cin_par1":4.8178,
                "hpv18_cin_fn_k":0.21695,
                "hpv18_dur_cin_par1":5.36992,
                "beta":0.05142,
                "f_cross_layer":0.0006683,
                "m_cross_layer":0.001476,
            })
param_list.append({"hpv16_cin_fn_k": 0.2541,  # No4 pars from d2Cal_8Sep25_6, fit 21.39
                "hpv16_dur_cin_par1":5.2199,
                "hpv18_cin_fn_k":0.27596,
                "hpv18_dur_cin_par1":4.7121,
                "beta":0.04471,
                "f_cross_layer":0.0003373,
                "m_cross_layer":0.0009538,
            })
"""





if __name__=="__main__":   

    #---Define interventions to add to sim ---#
    #TODO: do this and add interventions to sim for cal and sim(s) for results generation

    #---Set up a dummy calibration for me to be able to subsequently reset particular parmaeters---#
    sim = hpv.Sim(base_pars) 
    
    #Set up the dummy calibration
    calib_pars = dict( #to make the dummy calibration work, doing a cal with just beta. to make sure this script works as desired, make sure the dummy calibration is for a varaible which will be overriden with our hardcoded parameters - picking beta here as i think any calibration I do will involve beta, in which case it will be always overriden
            beta=[0.05, 0.00, 0.20], 
        )
    calib = hpv.Calibration(
        sim,
        calib_pars=calib_pars,
        datafiles= [#"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2012to27.csv",
                    "data\\d3\\mesherHPVCancerDist.csv",
                    ],
        total_trials=1,
        n_workers=1,
        keep_db=True,
        name=cal_name
    )
    calib.calibrate(die=False)


    #---Set up analyzer(s) for recording sim results corresponding to accumulated data---#
    az_cancers = hpv.age_results(
        result_args=sc.objdict(
            cancers=sc.objdict(
                datafile="C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2012to27.csv",
            ),
        )
    )


    #---Accumulate results over all parameter-sets and repeat runs of the sets---#
    sim_pars = calib.trial_pars_to_sim_pars() # Returns best parameters from calibration in a format ready for updating other sims' pars
    i=0
    all_cancer_results = []

    all_cancerous_genotype_results_16 = [] 
    all_cancerous_genotype_results_18 = [] 
    all_cancerous_genotype_results_hi5 = [] 

    gofs_cancers = []
    gofs_cancerous_genotypes = []

    for params in param_list:
        i+=1; print(f"Accumulating results for parameter set {i} of {len(param_list)}.")

        #Update sim_pars
        sim_pars['genotype_pars']['hpv16']['cin_fn']['k'] = params['hpv16_cin_fn_k']
        sim_pars['genotype_pars']['hpv16']['dur_cin']['par1'] = params['hpv16_dur_cin_par1']
        sim_pars['genotype_pars']['hpv18']['cin_fn']['k'] = params['hpv18_cin_fn_k']
        sim_pars['genotype_pars']['hpv18']['dur_cin']['par1'] = params['hpv18_dur_cin_par1']
        sim_pars['genotype_pars']['hi5']['cin_fn']['k'] = params['hi5_cin_fn_k']
        sim_pars['genotype_pars']['hi5']['dur_cin']['par1'] = params['hi5_dur_cin_par1']
        sim_pars['beta'] = params['beta']
        sim_pars['f_cross_layer'] = params['f_cross_layer']
        sim_pars['m_cross_layer'] = params['m_cross_layer']

        sim_pars['genotype_pars']['hi5']['rel_beta'] = params['hi5_rel_beta'] #NOTE: not sure about the inclusion of this in the cals!

        sim = hpv.Sim(base_pars, analyzers=[az_cancers])
        sim.update_pars(sim_pars)

        #Run multisim with {params}
        msim = hpv.MultiSim(sim, n_runs=repeat_runs)
        msim.run()

        #Extract results for each sim in the multisim
        for s in msim.sims:
            an = s.get_analyzer() #NOTE: if I am to add mutliple analyzers, I will need to change this slightly
            all_cancer_results.append(an.results)   

            an.result_args['cancers'].weights=np.ones(an.result_args['cancers'].data['value'].shape)
            gofs_cancers.append(an.compute_mismatch('cancers'))

            a,b,c = s.results['cancerous_genotype_dist'][:,(cancerous_genotype_dist_year-start)]
            all_cancerous_genotype_results_16.append(a)
            all_cancerous_genotype_results_18.append(b)
            all_cancerous_genotype_results_hi5.append(c)

            gof = hpv.misc.compute_gof(cancerous_genotype_dist_TRUEDATA, [a,b,c])
            gofs_cancerous_genotypes.append(gof.sum())
            
    # Convert results to array: shape (n_runs, n_timepoints) (n_runs=len(param_list)*repeat_runs)
    all_cancer_results = np.array(all_cancer_results)


    #---Plot cancer results ---#
    cancers_bins = [ 0., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80., 85., 90.]
    cancers = {} #cancers[year][bin_index] = [list of values across sims for bin cancers_bins[bin_index]]

    years = [2012,2013,2014,2015,2016,2017]

    for year in years:
        year_data = []
        for bin_index in range(len(cancers_bins)):
            bin_data = []
            for i in range(all_cancer_results.shape[0]):
                bin_data.append(all_cancer_results[i]['cancers'][year][bin_index])
            year_data.append(bin_data)
        cancers[year] = year_data
        

    #plot it by year
    box = 0
    rows = 3
    cols = 2
    fig, axs = plt.subplots(rows, cols)
    for year in years:
        row = box//cols; col = box%cols; box+=1
        ax = axs[row,col]

        cancer_year_data = cancers[year]
        ax.boxplot(cancer_year_data, tick_labels=cancers_bins)
        
        #Overlay plot from datafile
        datafilesource = f"C:\\Users\\fabia\\Documents\\Uni\\DPhil\\SophieHPV\\HPVsim_FABIAN\\data\\d2\\new_cervical_cancer_cases_{year}.csv"
        import pandas as pd
        df = pd.read_csv(datafilesource)
        true_data = df['value']
        true_data.tolist()
        for i in range(1, len(true_data)+1):
            ax.plot(i, true_data[i-1], 'ro')

        ax.set_title(f"Cancers by age bracket, {year}")
        
    plt.show()


    #--- Plot 'genotypes of cancers' results ---#
    fig,ax = plt.subplots(1,1)
    print([all_cancerous_genotype_results_16, all_cancerous_genotype_results_18, all_cancerous_genotype_results_hi5])
    ax.boxplot([all_cancerous_genotype_results_16, all_cancerous_genotype_results_18, all_cancerous_genotype_results_hi5],
                tick_labels=["HPV 16", "HPV 18", "HI5"])
    for i in range(1,4):
        ax.plot(i, cancerous_genotype_dist_TRUEDATA[i-1], 'ro')
    ax.set_title(f"Cancerous Genotype Distribution, {cancerous_genotype_dist_year}")
    plt.show()


    #--- Plot GOF data ---#
    fig, ax = plt.subplots(1,1)
    ax.boxplot([gofs_cancers, gofs_cancerous_genotypes], 
               tick_labels=["cancers", "cancerous genotype distribution"])
    ax.set_title("GOFs")
    plt.show()
    
    