# Import HPVsim
import hpvsim as hpv
import numpy as np

import numpy as np

import pickle

#Importing interventions from specific models of NHS England interventions
from InterventionAlgorithms import NHS_Vacc as vacc
#TODO: implement each era of NHS england screening algorithm in seperate files and test them and import them



#Define sexual mixing matrices according to [TODO: find sources sophie used to inform the mixing matrices she made] 
sophie_married_matrix = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [10, 0, 0, 0.08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [15, 0, 0, 0.08, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [20, 0, 0, 0, 0, 0.6, 2, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0],
    [25, 0, 0, 0, 0, 0.6, 1, 2, 0.4, 0.1, 0, 0, 0, 0, 0, 0, 0],
    [30, 0, 0, 0, 0, 0.5, 0.5, 2, 1, 0.5, 0.1, 0, 0, 0, 0, 0, 0],
    [35, 0, 0, 0, 0, 1, 0.5, 1, 2, 1, 0.5, 0.2, 0, 0, 0, 0, 0],
    [40, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0],
    [45, 0, 0, 0, 0, 0.1, 1, 2, 2, 2, 1, 0.5, 0.2, 0.08, 0, 0, 0],
    [50, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 2, 2, 0.5, 0.2, 0.05, 0, 0],
    [55, 0, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 3, 2, 1, 0.3, 0.1, 0.1],
    [60, 0, 0, 0, 0, 0, 0, 0.1, 0.5, 1, 2, 3, 3, 2, 0.5, 0.3, 0.1],
    [65, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 2, 3, 3, 2, 1, 0.2],
    [70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 3, 2, 1],
    [75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3],
]
sophie_married_matrix = np.array(sophie_married_matrix)
sophie_casual_matrix = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
       [5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
       [10,0,0,0,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0],
       [15,0,0,1,2,3,2,1,0.5,0,0,0,0,0,0,0,0],
       [20,0,0,0.15,2,3,2,2,1,0.15,0,0,0,0,0,0,0],
       [25,0,0,0.15,0.25,1,2,2,1,1,0,0,0,0,0,0,0],
       [30,0,0,0,0,0.5,0.5,2,1,0.15,0,0,0,0,0,0,0],
       [35,0,0,0,0,1,0.5,1,2,1,0.5,0,0,0,0,0,0],
       [40,0,0,0,0,1,1,1,1,1,0.5,0.25,0,0,0,0,0],
       [45,0,0,0,0,0.15,1,2,2,2,1,0.5,0.2,0.1,0,0,0],
       [50,0,0,0,0,0,0.15,1,2,3,2,2,0.5,0.2,0.05,0,0],
       [55,0,0,0,0,0,0,0.15,1,2,3,3,2,1,0.25,0.1,0.1],
       [60,0,0,0,0,0,0,0.15,0.15,1,2,3,3,2,0.5,0.25,0.1],
       [65,0,0,0,0,0,0,0,0,0,1,1,2,2,1,0.5,0],
       [70,0,0,0,0,0,0,0,0,0,0,0,0,0.8,1,0.7,0.5],
       [75,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.25]
]
sophie_casual_matrix = np.array(sophie_casual_matrix)

if __name__=="__main__":
    #---SET UP SIMULATION TO CALIBRATE---#
    pars = dict(n_agents= 100e3,
                start=1980, end=2022, dt=0.25, #start=1980, end=2020, dt=0.25, 
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
                }, #(note, this measure will be rescaled to a prob distribution by hpvsim.utils.choose_w)

                interventions = vacc.vaccinations + [],

                #TODO: get analysers which can test my intiial considitons are what i want them to be, that is, I am setting the parameters correctly

                )


    #---Implemeting screening algorithm(s) to match historic+projected screening ---#
    #TODO: once my NHS_Screening_Pathway.py code is all good, I can pretty much just copy it over here

    sim = hpv.Sim(pars) 

    #---SET UP CALIBRATION---#

    # Configure a simulation with some parameters
    

    # Specify some parameters to adjust during calibration.
    # The parameters in the calib_pars dictionary don't vary by genotype,
    # whereas those in the genotype_pars dictionary do. Both kinds are
    # given in the order [best, lower_bound, upper_bound].

    #I have added new extended ranges where the best parameters from the [B] cals are near to the sides of the ranges. old ranges commented to the right  
    calib_pars = dict(
            beta=[0.25,0.00,0.50],#[0.05, 0.00, 0.20], #still happy
            f_cross_layer= [0.15, 0, 1], #always close to 0 but as i think this cant be negaative, still happy
            m_cross_layer= [0.25, 0, 1], #pretty close to 0, not quite aas much aas f, but still i think i cnt extend the range so happy (for 3, this is also super super close to 0)
            
        )

    genotype_pars = dict(
        hpv16=dict(
            cin_fn=dict(k=[0.5, 0.0, 1.0]),#(k=[0.5, 0.2, 1.0]), haappy now they can be around 0.1 - 3 agrees
            dur_cin=dict(par1=[6, 1, 12])#dict(par1=[6, 4, 12]) perhaaps push this even lower to 1 min, it seems aaraound 2-3.5 but still. 3 cnt seem to settle here, maybe doesnt maatter thaat much
        ),
        hpv18=dict(
            cin_fn=dict(k=[0.5, 0.0, 1.0]),#dict(k=[0.5, 0.2, 1.0]), haappy now these cn be v small for both caals, it does seem to waant to be quite small
            dur_cin=dict(par1=[6, 1, 12])#dict(par1=[6, 4, 12]) these seem to waant to be big in 2 and 3, but i wamnder if th9is is because of asymmtery with range of 16, so first lets make them match
        ),
        hi5=dict(
            cin_fn=dict(k=[0.5, 0.0, 1.0]),
            dur_cin=dict(par1=[6, 1, 12]),
            rel_beta=[0.9,0,1]
        )
    )

    # List the datafiles that contain data that we wish to compare the model to:
    datafiles=[
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d1/new_cervical_cancer_cases_condensed.csv"]
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d1/new_cervical_cancer_cases_2017.csv"]
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2017.csv"] #for thiss data, i am assuming 'cancers' is the number of newly detected cancers in a year; but I suppose the same assumption is the case if I were to use 'cancer_incidence'
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_UNIFIED.csv"] #for this data, i am assuming 'cancers' is the number of newly detected cancers in a year; but I suppose the same assumption is the case if I were to use 'cancer_incidence'
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2012to27.csv"]
        
        #"data\\d2\\new_cervical_cancer_cases_2012to27.csv",
        "data\\d3\\new_cervical_cancer_cases_ENGSCALED1P19TOUK.csv",
        "data\\d3\\mesherHPVCancerDist.csv", 
        ]

    # List extra results that we don't have data on, but wish to include in the
    # calibration object so we can plot them.
    results_to_plot = ['cancer_incidence', 'asr_cancer_incidence']

    # Create the calibration object, run it, and plot the results
    calib = hpv.Calibration(
        sim,
        calib_pars=calib_pars,
        genotype_pars=genotype_pars,
        extra_sim_result_keys=results_to_plot,
        datafiles=datafiles,

        total_trials=5000, #50
        n_workers=12, #5

        keep_db=True,
        name="CalibrationRawResults\\d2Cal_17Sep25_XPS_F1VACCONLY_1"

    )

    #---PERFORM CALIBRATION---#
    calib.calibrate(die=False,
                    #plots=["learning_curve", "timeline"],                 
                                   # detailed_contplot=[('f_cross_layer','m_cross_layer')], #can't do contour for everything; uses a lot of memory. Likely can't do it for a very large calibration either!
    )

    #---ASSESS CALIBRATION QUALITY---#
    #Plot goodness of fit compared to observed data
    calib.plot(res_to_plot=50)
    calib.plot(res_to_plot=10)
    calib.plot(res_to_plot=5)
    calib.plot(res_to_plot=1)

    #Additional plotting to assess how the calibration went
    #calib.plot_learning_curve()
    #calib.plot_timeline()

    #---SAVE CALIBRATED PARAMETERS---#
    print(calib.df.head(10))

    