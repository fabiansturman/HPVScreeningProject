# Import HPVsim
import hpvsim as hpv
import numpy as np

import numpy as np

import pickle




#NOTE: if results from a calibration aren't satisfyingly tight, may need to check over, say, calibration.py in Sophie's code to see...
#... if her other chanes to either the sim parameters, or the genotype_pars dictionary, improve things!

#TODO: definitely worth going over that table at the end of sophie's paper, to find justifications for any parameter choices she has not already justified, and then manually set sim parameters here to match the UK as Sophie has done!

#TODO: also deffo need to get gentoype data into here!

#TODO: worth calibrating sev_dist, which can be seen as an immunocompromise factor controlling how quickly we progress thorugh the disease

#TODO: once I have Sophie's BSP, look into her appendices for the values for the mixing matrix for unmarried individuals. The married one I have, but I am not sure where she has gotten this info - is it specific to england (as in, found though lit review or calibration?) or just a default from within HPVsim? I can at least deduce whether it is a HPVsim default/loaded in as part of teh demographic info!

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
    pars = dict(n_agents= 100e3,#100e3,#40e4,#20e4, #20e3, 
                start=1980, end=2017, dt=0.25, #start=1980, end=2020, dt=0.25, 
                location='united kingdom', 
                verbose=-1,
                debut=dict(f=dict(dist='normal', par1=16.0, par2=3.1), m=dict(dist='normal', par1=16.0, par2=4.1)),
                mixing = {'m':sophie_married_matrix,
                          'c':sophie_casual_matrix},
                condoms = dict(m=0.01, c=0.2),  #condom usage in (m)arried and (c)asual relationships
                #NOTE: I think having too many young cancers is probably caused by the debut age being too low <-or it could just be the whole data of that being in shambles, now we have nowhere near enough younger cancers
                    #^or because of the mixing matrices - looks like perhaps adding sophies mixing matrices helps things?
                )



    # Implement Current NHS Algorithm (): defining inteventions to screen, triage, assign treatment, and administer treatment
    
    #TODO: gradually add intervnetions in, to match the real NHS screening programme. The stuff below is a crude approximation of the programme but definitley not the programme itself.

    start_year = 2015
    primary_screen_prob = 0.2
    triage_screen_prob = 0.9
    ablate_prob = 0.9
    prob = 0.6
    screen_eligible = lambda sim: np.isnan(sim.people.date_screened) | (sim.t > (sim.people.date_screened + 3 / sim['dt']))

    #Algorithm 7. HPV DNA as the primary screening test, followed by cytology triage, followed by colposcopy and treatment.
    hpv_primary7  = hpv.routine_screening(eligibility=screen_eligible, start_year=start_year, prob=primary_screen_prob, product='hpv', label='hpv_primary7') #HPV Testing
    # Send HPV+ women* for cytology
    to_cytology = lambda sim: sim.get_intervention('hpv_primary7').outcomes['positive']
    cytology7   = hpv.routine_triage(eligibility=to_cytology, prob=triage_screen_prob, product='lbc', annual_prob=False, label='cytology7')
    #Send ASCUS and abnormal cytology results for colpo.
    to_colpo7 = lambda sim:list(set(sim.get_intervention('cytology7').outcomes['abnormal'].tolist() + sim.get_intervention('cytology7').outcomes['ascus'].tolist())) #Define who's eligible for colpo
    colpo7 = hpv.routine_triage(eligibility=to_colpo7, prob=triage_screen_prob, product='colposcopy', annual_prob=False, label='colposcopy') #Send people to colposcopy
    #After colpo, treat HSILs with Ablation
    hsils7 = lambda sim: sim.get_intervention('colposcopy').outcomes['hsil'] # Define who's eligible for ablation
    ablation7  = hpv.treat_num(eligibility=hsils7, prob=ablate_prob, product='ablation', label='ablation') # Administer ablation
    
    
    
    #sim = hpv.Sim(pars, interventions = [hpv_primary7, cytology7, colpo7, ablation7])

    sim = hpv.Sim(pars) #<-simulation without interventions

    #---SET UP CALIBRATION---#

    # Configure a simulation with some parameters
    

    # Specify some parameters to adjust during calibration.
    # The parameters in the calib_pars dictionary don't vary by genotype,
    # whereas those in the genotype_pars dictionary do. Both kinds are
    # given in the order [best, lower_bound, upper_bound].

    #I have added new extended ranges where the best parameters from the [B] cals are near to the sides of the ranges. old ranges commented to the right  
    calib_pars = dict(
            beta=[0.05, 0.00, 0.20], #still happy
            f_cross_layer= [0.05, 0, 2], #always close to 0 but as i think this cant be negaative, still happy
            m_cross_layer= [0.05, 0, 2], #pretty close to 0, not quite aas much aas f, but still i think i cnt extend the range so happy (for 3, this is also super super close to 0)
            #init_hpv_prev
        )

    genotype_pars = dict(
        hpv16=dict(
            cin_fn=dict(k=[0.5, 0.0, 1.0]),#(k=[0.5, 0.2, 1.0]), haappy now they can be around 0.1 - 3 agrees
            dur_cin=dict(par1=[6, 1, 12])#dict(par1=[6, 4, 12]) perhaaps push this even lower to 1 min, it seems aaraound 2-3.5 but still. 3 cnt seem to settle here, maybe doesnt maatter thaat much
        ),
        hpv18=dict(
            cin_fn=dict(k=[0.5, 0.0, 1.0]),#dict(k=[0.5, 0.2, 1.0]), haappy now these cn be v small for both caals, it does seem to waant to be quite small
            dur_cin=dict(par1=[6, 1, 12])#dict(par1=[6, 4, 12]) these seem to waant to be big in 2 and 3, but i wamnder if th9is is because of asymmtery with range of 16, so first lets make them match
        )
        #TODO: dont i want to do the sme with hi5 then??
    )

    # List the datafiles that contain data that we wish to compare the model to:
    datafiles=[
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d1/new_cervical_cancer_cases_condensed.csv"]
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d1/new_cervical_cancer_cases_2017.csv"]
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2017.csv"] #for this data, i am assuming 'cancers' is the number of newly detected cancers in a year; but I suppose the same assumption is the case if I were to use 'cancer_incidence'
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_UNIFIED.csv"] #for this data, i am assuming 'cancers' is the number of newly detected cancers in a year; but I suppose the same assumption is the case if I were to use 'cancer_incidence'
        #"C:/Users/fabia/Documents/Uni/DPhil/SophieHPV/HPVsim_FABIAN/data/d2/new_cervical_cancer_cases_2012to27.csv"]
        
        "data\\d2\\new_cervical_cancer_cases_2012to27.csv"]

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
        n_workers=5,#12, #2

        keep_db=True,
        name="CalibrationRawResults\\d2Cal_12Sep25_XPS_D6"
  
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

    