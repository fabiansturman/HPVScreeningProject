# Import HPVsim
import hpvsim as hpv
import numpy as np

import numpy as np

from scenario_configuration import calibration_code  #import the basic name of this calibration according to the parameterisation of the model we are calibrating to


#Importing interventions from specific models of NHS England interventions
from InterventionAlgorithms import screeningAlgorithms, NHS_Vacc
from basePars import base_pars

#NOTE: this calibration is using the NHS alg of 2025+ as an approximation of the full history of NHS screening - note that my F caloibrations show that this works when trainign to all the data, so I may as well go and train on just the start of teh data right now!

if __name__=="__main__":
    #---SET UP SIMULATION TO CALIBRATE---#
    pars = base_pars
    pars['interventions'] =  screeningAlgorithms.get_interventions(v=5,u=5,a=5)  + NHS_Vacc.vaccinations
    pars['end'] = 2025
    
    sim = hpv.Sim(pars) 
    
    #sim.run()
    #quit()

    #---SET UP CALIBRATION---#

    # Configure a simulation with some parameters
    
    #I have added new extended ranges where the best parameters from the [B] cals are near to the sides of the ranges. old ranges commented to the right  
    calib_pars = dict(
            beta=[0.25,0.00,0.50],#[0.05, 0.00, 0.20], #still happy
            f_cross_layer= [0.15, 0, 1], #always close to 0 but as i think this cant be negaative, still happy
            m_cross_layer= [0.25, 0, 1], #pretty close to 0, not quite aas much aas f, but still i think i cnt extend the range so happy (for 3, this is also super super close to 0)
            
        )

    genotype_pars = dict(
        hpv16=dict(
            cin_fn=dict(k=[0.25, 0.0, 1.0]),#(k=[0.5, 0.2, 1.0]), haappy now they can be around 0.1 - 3 agrees
            dur_cin=dict(par1=[5, 1, 12])#dict(par1=[6, 4, 12]) perhaaps push this even lower to 1 min, it seems aaraound 2-3.5 but still. 3 cnt seem to settle here, maybe doesnt maatter thaat much
        ),
        hpv18=dict(
            cin_fn=dict(k=[0.25, 0.0, 1.0]),#dict(k=[0.5, 0.2, 1.0]), haappy now these cn be v small for both caals, it does seem to waant to be quite small
            dur_cin=dict(par1=[5, 1, 12])#dict(par1=[6, 4, 12]) these seem to waant to be big in 2 and 3, but i wamnder if th9is is because of asymmtery with range of 16, so first lets make them match
        ),
        hi5=dict(
            cin_fn=dict(k=[0.25, 0.0, 1.0]),
            dur_cin=dict(par1=[5, 1, 12]),
            rel_beta=[0.9,0,1]
        ),
        ohr=dict(
            cin_fn=dict(k=[0.25, 0.0, 1.0]),
            dur_cin=dict(par1=[5, 1, 12]),
            rel_beta=[0.9,0,1]
        )
    )

    # List the datafiles that contain data that we wish to compare the model to:
    datafiles=["project_modelling/real_world_data/new_cervical_cancer_cases_ENGSCALED1P19TOUK_START.csv",
               "project_modelling/real_world_data/mesherHPVCancerDist.csv"]

    """
    Current cals that are running that I am yet to process into final cals (with their sibling cals)

    -> Apr12_2__8_68_76_55_55_9_9__96_7_13 : COMPELTE this is a cal (on MB) using Sophie's mixing matrices and base sensitivity analysis case otherwise
    -> Apr12_4__8_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (done on BMRC, transferred to MB) using Sophie's mixing matrices and base sensitivity analysis case otherwise
    -> Apr13_1__8_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (done on BMRC, transferred to MB) using Sophie's mixing matrices and base sensitivity analysis case otherwise
        => ** summarised into FinalPooledCalibration_8_68_76_55_55_9_9__96_7_13_OLDMIXING.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:
        -show this result still in the appendix

    -> Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (on MB), using my own mixing matrices (Married dervied from natsal cohabitation data, and Casual according to sexual encounter mixing matrix and base sensitivuty analysis case otehrwise (the sexual encounter mixing matrix hopefully is a good proxy for causal - looks good intuitively as it is a diagnonal with a bit larger entropy than marriage -but  really also includes marriage but hopefully other cal fits for the doubling? maybe not hmm but we will see. I guess we could calibrate condom usage in marriage to serve as a proxy to allow us to alter the effect of the married layer to cover for doubling up?))
    -> Apr15_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (done on BMRC, transferred to MB), identical setup to the calibration directly above
    -> Apr15_2_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (on MB) identical setup to above
        => ** summarised into FinalPooledCalibration_8_68_76_55_55_9_9__96_7_13.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:
        - this is the preferred choice of mixing matrix and what all the other sensitivtiy analyses will use (as well as the base case, which is this here)

    -> Apr16_1__65_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (on BMRC, transferred to MB) doing a Sensitivity Analysis on cancer treatment effectiveness (trying 65% rather than 80%) (otherwise, it is identical to the setup in Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13)
    -> Apr17_1__65_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (on MB) the same as above
    -> Apr17_2__65_68_76_55_55_9_9__96_7_13 : COMPLETE this is a cal (on BMRC, transferred to MB) the same as above
        => ** summarised into FinalPooledCalibration_65_68_76_55_55_9_9__96_7_13.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:

    -> Apr18_1__8_68_76_36_36_8_8__96_7_13 : COMPLETE this is a cal (on MB) doing a Sensitivity Analysis that is a lower bound on Screening Algorothm Compliance (otherwise is identical to the setup in Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13)
    -> Apr18_2__8_68_76_36_36_8_8__96_7_13 : COMPLETE this is a cal (on XPS, transferred to MB) as above
    -> Apr18_4__8_68_76_36_36_8_8__96_7_13 : COMPLETE this is a cal (on BMRC, transferred to MB) as above
        => ** summarised into FinalPooledCalibration_8_68_76_36_36_8_8__96_7_13.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:
    
    -> Apr19_1__8_8_8_68_68_95_95__96_7_13 : COMPLETE this is a cal (on MB) doing a Sensitivity Analysis that is an a upper bound on Screening Algorothm Compliance (otherwise is identical to the setup in Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13)
    -> Apr20_1__8_8_8_68_68_95_95__96_7_13 : COMPLETE this is a cal (on BMRC, copied to MB) as above
    -> Apr20_2__8_8_8_68_68_95_95__96_7_13 : COMPLETE this is a cal (on MB) as above
        => ** summarised into FinalPooledCalibration_8_8_8_68_68_95_95__96_7_13.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:

    -> Apr22_3__8_68_76_55_55_9_9_1_96_7_13 : COMPLETE this is a cal (on BMRC, on MB) doing a Sensitivity Analysis considering ohr rel imm of 0.1. (otherwise is identical to the setup in Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13)
    -> Apr23_1__8_68_76_55_55_9_9_1_96_7_13 : INPROGRESS this is a cal (on XPS) as above     
    -> Apr23_3__8_68_76_55_55_9_9_1_96_7_13 : INPROGRESS this is a cal (on MB) as above

    -> Apr21_1__8_68_76_55_55_9_9__96_7_62 : COMPLETE this is a cal (on XPS, copied to MB) doing a Sensitivity Analysis on Initial HPV Prevalence among Men, seeing what will happen if we make it drop off with age like we do for women (otherwise is identical to the setup in Apr14_1_USINGNATSALANDLIT__8_68_76_55_55_9_9__96_7_13)
    -> Apr21_2__8_68_76_55_55_9_9__96_7_62 : COMPLETE this is a cal (on BMRC, copied to MB) as above
    -> Apr22_2__8_68_76_55_55_9_9__96_7_62 : COMPLETE this is a cal (on MB) as above
        => ** summarised into FinalPooledCalibration_8_68_76_55_55_9_9__96_7_62.pickle **
                ^ GOF details of best cals:
                    min:
                    LQ:
                    median:
                    UQ:
                    max:
    
    -> Apr23_4
    ->
    ->


    TODO: -> complete all remaining senstiivty analysis calibrations 
        NOTE: on reconsideration, I am NOT going to compare calibration evaluations between different sensitivity analyses unless there is a clear difference (we are talking, UQ of one set of SA GOFs is below LQ of the main cases's GOFs). Otherwise, it is a totally different problem to be thinking about what sensitivty analysis most makes sense, especially as my main case is informed by the literature with SAs currently probing other potential parameterisations of which wew are less certain of the correctness. And besides, it is very very likely that my conclusions end up being the same across all the SAs. So unless the conclusions start being quite different between the SAs or one of the SAs ends up with an evaluation GOF distribution that is MUCHH better than the main case's evaluation GOF distribution (implying i may be way off with my main case, or at least way off how good I could have been in my fit), then keep the cases as they are right now

        TODO: Compare my GOF in the end to my previous cal (i.e. getting sub 0.2 gofs in 75% at least of cases, ideally 100%). i should be there!!
    """

    # Create the calibration object, run it, and plot the results
    name =  f"project_modelling/calibration_results/Apr23_4_{calibration_code}" # "calibration_results/Feb14_1" #"calibration_results/Feb13_1"
    
    calib = hpv.Calibration(
        sim,
        calib_pars=calib_pars,
        genotype_pars=genotype_pars,
        datafiles=datafiles,

        total_trials=5000,
        n_workers=5, #to be changed according to the multiprocessing capability of the hardware being used

        keep_db=True,
        name=name
    )

    #---PERFORM CALIBRATION---#
    calib.calibrate(die=False)

    #---ASSESS CALIBRATION QUALITY---#
    #Plot goodness of fit compared to observed data
    calib.plot(res_to_plot=50)
    calib.plot(res_to_plot=10)
    calib.plot(res_to_plot=5)
    calib.plot(res_to_plot=1)

    #---OUTPUT AND SAVE CALIBRATED PARAMETERS---#
    print(calib.df.head(10))

    calib.df.to_csv(f"{name}.csv")

    #TODO: Once I have settled on a final initial setup for calibration: any changes I make to how I calibrate must be shown where relevant in the excel paramter spreadsheet!!

    #TODO: delete all practice/test calibrations from my local copy of the git repo and then change the .gitignore so that my calibration files get uploaded 

    #TODO: cal evaluation shpuld include plots just along an axis for each param with dots on the axis for the sampled parameter positions, with the size + height of teh dot according to its GOF so the goopd ones are most visible in a sea of samples, and we colouyr dots by the calibration we got them from, so we can see how cals differ in the parameters they pick, in the marginals (some persoective, if only in the marginals, is better than none)
