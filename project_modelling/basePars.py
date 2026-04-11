"""
Stores the base parameters for a HPVsim simulation using the current NHS strategy.
Can override aspects of this dictionary, particularly interventions, if needed.
"""
import numpy as np
from InterventionAlgorithms import NHS_2025_lambdamu, NHS_Vacc



married_matrix = [        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [10, 0, 0, 0.08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [15, 0, 0, 0.08, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],        [20, 0, 0, 0, 0, 0.6, 2, 0.2, 0.1, 0, 0, 0, 0, 0, 0, 0, 0],        [25, 0, 0, 0, 0, 0.6, 1, 2, 0.4, 0.1, 0, 0, 0, 0, 0, 0, 0],        [30, 0, 0, 0, 0, 0.5, 0.5, 2, 1, 0.5, 0.1, 0, 0, 0, 0, 0, 0],        [35, 0, 0, 0, 0, 1, 0.5, 1, 2, 1, 0.5, 0.2, 0, 0, 0, 0, 0],        [40, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0.5, 0.3, 0.1, 0, 0, 0, 0],        [45, 0, 0, 0, 0, 0.1, 1, 2, 2, 2, 1, 0.5, 0.2, 0.08, 0, 0, 0],        [50, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 2, 2, 0.5, 0.2, 0.05, 0, 0],        [55, 0, 0, 0, 0, 0, 0, 0.1, 1, 2, 3, 3, 2, 1, 0.3, 0.1, 0.1],        [60, 0, 0, 0, 0, 0, 0, 0.1, 0.5, 1, 2, 3, 3, 2, 0.5, 0.3, 0.1],        [65, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 2, 3, 3, 2, 1, 0.2],        [70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 1, 2, 3, 3, 2, 1],        [75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 3],    ]
married_matrix = np.array(married_matrix)
casual_matrix = [[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],        [5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],        [10,0,0,0,0.2,0.1,0.05,0,0,0,0,0,0,0,0,0,0],        [15,0,0,1,2,3,2,1,0.5,0,0,0,0,0,0,0,0],        [20,0,0,0.15,2,3,2,2,1,0.15,0,0,0,0,0,0,0],        [25,0,0,0.15,0.25,1,2,2,1,1,0,0,0,0,0,0,0],        [30,0,0,0,0,0.5,0.5,2,1,0.15,0,0,0,0,0,0,0],        [35,0,0,0,0,1,0.5,1,2,1,0.5,0,0,0,0,0,0],        [40,0,0,0,0,1,1,1,1,1,0.5,0.25,0,0,0,0,0],        [45,0,0,0,0,0.15,1,2,2,2,1,0.5,0.2,0.1,0,0,0],        [50,0,0,0,0,0,0.15,1,2,3,2,2,0.5,0.2,0.05,0,0],        [55,0,0,0,0,0,0,0.15,1,2,3,3,2,1,0.25,0.1,0.1],        [60,0,0,0,0,0,0,0.15,0.15,1,2,3,3,2,0.5,0.25,0.1],        [65,0,0,0,0,0,0,0,0,0,1,1,2,2,1,0.5,0],        [70,0,0,0,0,0,0,0,0,0,0,0,0,0.8,1,0.7,0.5],        [75,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.25]    ]
casual_matrix = np.array(casual_matrix)

start = 1980
end = 2055

base_pars = dict(n_agents= 200_000, 
                start=start, end=end, dt=0.25, 
                location='united kingdom', 
                verbose=-1,
                debut=dict(f=dict(dist='normal', par1=16.0, par2=3.1), m=dict(dist='normal', par1=16.0, par2=4.1)),
                mixing = {'m':married_matrix,
                          'c':casual_matrix},
                condoms = dict(m=0.17, c=0.50), #condom usage in (m)arried and (c)asual relationships
                
                genotypes     = ['hpv16', 'hpv18', 'hi5', 'ohr'],

                init_hpv_prev = {
                    'age_brackets'  : np.array([  16,   24,   34,   44,  54,   64, 150]),
                    'f'             : np.array([ 0.0, 0.25, 0.14,   0.08,  0.06,   0.06, 0.03]),
                    'm'             : np.array([ 0.0,0.20, 0.23,   0.20,  0.18,   0.25, 0.25])
                },

                init_hpv_dist = {
                    'hpv16': 2.3,
                    'hpv18': 0.9,
                    'hi5':  2.2, #HPV 33 is not listed as one of the top 10 most prevalent in general population in (), so we can assume its prevalence is at most 0.004 - so not adding this to the sum
                    'ohr': 2.1,
                }, #(note, this measure will be rescaled to a prob distribution by hpvsim.utils.choose_w)

                #TODO: get analysers which can test my intiial considitons are what i want them to be, that is, I am setting the parameters correctly
                burnin = 20,
                )


if __name__ == "__main__":
    import hpvsim as hpv

    base_pars['verbose'] = 1

    sim = hpv.Sim(base_pars)
    sim.run()
    sim.plot()