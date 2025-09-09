#Making my own products: my own screening, treatment, and vaccination, and testing them in a 

#This code is based on: the HPVsim tutorial, the documentation, and the .csv file data for tx, dx, and vx products

import hpvsim as hpv
import pandas as pd

if __name__ == "__main__":
    ##Making my own products##

    #Treatment products (tx) are defined by:

    #Vaccination products (vx) are defined by:
  #My Vaccination products (vx) are defined by:
    sophies_vx1_data = pd.DataFrame({'name':'Sophie\'s vaccine ',
                                    'genotype':['hpv16', 'hpv18'],
                                    'rel_imm':[0.75,0.75]
                                    })
    
    sophies_vx1 = hpv.vx(genotype_pars=sophies_vx1_data, #information from our DataFrame (basically an excel clipping) about our vaccination's effectiveness against particular genotypes (only need to specify the relevant genotypes!)
                            imm_init={'dist':'normal', 'par1':0.8, 'par2':0.02} #samples peak immunity from a distribution as defined in the list - see utils.sample for more info
                                )

    #My Vaccination products (vx) are defined by:
    sophies_vx2_data = pd.DataFrame({'name':'Sophie\'s vaccine ',
                                        'genotype':['hpv16', 'hpv18'],
                                        'rel_imm':[0.9,0.9]
                                        })
        
    sophies_vx2 = hpv.vx(genotype_pars=sophies_vx2_data, #information from our DataFrame (basically an excel clipping) about our vaccination's effectiveness against particular genotypes (only need to specify the relevant genotypes!)
                            imm_init={'dist':'normal', 'par1':0.8, 'par2':0.02} #samples peak immunity from a distribution as defined in the list - see utils.sample for more info
            )
    sophies_vx3_data = pd.DataFrame({'name':'Fabian\'s vaccine ',
                                        'genotype':['hpv16', 'hpv18'],
                                        'rel_imm':[1,1]
                                        })
        
    sophies_vx3 = hpv.vx(genotype_pars=sophies_vx3_data, #information from our DataFrame (basically an excel clipping) about our vaccination's effectiveness against particular genotypes (only need to specify the relevant genotypes!)
                            imm_init={'dist':'normal', 'par1':0.9, 'par2':0.03} #samples peak immunity from a distribution as defined in the list - see utils.sample for more info
                            )
     #IF WE WANT, WE CAN ALWAYS ADD INTERVENTIONS TO THE EXCEL FILES IN hpvsim.data AS THESE FILES ARE LOADED IN EACH TIME!!
     # (...this bypasses the need to create the classes here, though it is messier - ideally those excel files should only contain real life data!!) 

    ##Define the interventions themselves that we are using##       . Note how we assign the interventions labels which are chained together to form a treatment process/algorithm
    prob = 0.9        #(For simplicity, I am using the same probability for each intervention)
    
        ##Set up a prophylactic vaccination##
    vx1 = hpv.routine_vx(prob=prob, start_year=2000, age_range=[9,15], product=sophies_vx1)
    vx2 = hpv.routine_vx(prob=prob, start_year=2000, age_range=[9,15], product=sophies_vx2)
    vx3 = hpv.routine_vx(prob=prob, start_year=2000, age_range=[9,15], product=sophies_vx3)

        ##Create and run sims, plot data##
    pars = dict(
            start       = 1975,
            use_waning  = True,
            n_agents    = 20e3,
            n_years     = 100,
            verbose     = 0,
            rand_seed   = 1,            #This sets a non-default random seed
            genotypes   = [16,18]       #We are including the two genootypes of greatest general interest in our simulation
        )

    baseline_sim = hpv.Sim(pars, label="No Vaccinaction")
    sim2 = hpv.Sim(pars, interventions = [vx2] , label='2-Doses of the Nonavalent Gardasil9® Vaccination (9vHPV)')
    sim1 = hpv.Sim(pars, interventions = [vx1] , label='1-Dose of the Nonavalent Gardasil9® Vaccination (9vHPV)')
    sim3 = hpv.Sim(pars, interventions = [vx3] , label='3-Doses of the Nonavalent Gardasil9® Vaccination (9vHPV)'
                )

    msim = hpv.parallel(baseline_sim, sim1, sim2, sim3)
    msim.plot(style='ggplot')