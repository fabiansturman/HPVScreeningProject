#Making my own products: my own screening, treatment, and vaccination, and testing them in a 

#This code is based on: the HPVsim tutorial, the documentation, and the .csv file data for tx, dx, and vx products

import hpvsim as hpv
import pandas as pd

if __name__ == "__main__":
    rel_imm_values = list(range(50, 101))  # Convert the range 50-100 to a list of integers (0.5 to 1)
    for rel_imm_value in rel_imm_values:
        rel_imm = rel_imm_value / 100.0  # Convert back to a float between 0.5 and 1
        
    ##Making my own products##
    #My Vaccination products (vx) are defined by:
    sophies_vx_data = pd.DataFrame({'name':'Sophie\'s vaccine ',
                                    'genotype':['hpv16', 'hpv18'],
                                    'rel_imm':[rel_imm,rel_imm]
                                    })
    
    sophies_vx = hpv.vx(genotype_pars=sophies_vx_data, #information from our DataFrame (basically an excel clipping) about our vaccination's effectiveness against particular genotypes (only need to specify the relevant genotypes!)
                        imm_init={'dist':'normal', 'par1':0.8, 'par2':0.02} #samples peak immunity from a distribution as defined in the list - see utils.sample for more info
                            )
    
    ##Define the interventions themselves that we are using##       . Note how we assign the interventions labels which are chained together to form a treatment process/algorithm
    prob = 0.9        #(For simplicity, I am using the same probability for each intervention)

        ##Set up my vaccination##
    vx = hpv.routine_vx(prob=prob, start_year=2000, age_range=[9,15], product=sophies_vx)

        ##Create and run sims, plot data##
    pars = dict(
            start       = 1965,
            use_waning  = True,
            n_agents    = 20e3,
            n_years     = 80,
            verbose     = 0,
            rand_seed   = 1,            #This sets a non-default random seed
            genotypes   = [16,18]       #We are including the two genootypes of greatest general interest in our simulation
        )

    baseline_sim = hpv.Sim(pars, label="No Vaccinaction")
    sim = hpv.Sim(pars, interventions = [vx] , label='Vaccination')

    msim = hpv.parallel(baseline_sim, sim)
    msim.plot(style='ggplot');