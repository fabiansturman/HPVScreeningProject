"""
The purpose of this code is to make an intervention which
"""

import hpvsim as hpv

import numpy as np

from tqdm import tqdm

# Define a custom intervention
class CureCancer(hpv.Intervention):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def apply(self, sim):
        # Identify people with cancer
        people = sim.people
        cancer_inds = people.cancerous > 0  # people with any cancer
        # Cure them
        
        for genotype in range(people.cancerous.shape[0]):
            for pid in range(people.cancerous.shape[1]):
                people.cancerous[genotype, pid] = False

        #people.cancerous[cancer_inds] = 0
        #people.dead[cancer_inds] = False  # ensure they don't die from it

def cure_cancers(sim):
    """
    #The below DOES work to make everyone not infectious and therefore set cancers to 0, great
    for genotype in range(sim.people.infectious.shape[0]):
        for pid in range(sim.people.infectious.shape[1]):
            sim.people.infectious[genotype, pid] = False
    """
    if sim.t % 5 == 0 or True:#50 == 0:
        for pid in tqdm(range(len(sim.people))):
            sim.people[pid].cancerous=np.array([False, False, False])

    if sim.t==120 and False:
        for pid in range(len(sim.people)):
            print(sim.people[pid])
            """
            When printing the person, if we screen them using routine screening that DOES get recognised in thier internal state, good.
            """
            quit(0)

    for genotype in range(sim.people.cancerous.shape[0]):
        for pid in range(sim.people.cancerous.shape[1]):
            sim.people.cancerous[genotype, pid] = False



def eligible(sim):
    for_screening = np.isnan(sim.people.date_screened) & (sim.people.age>=24) & (sim.people.age<=60) #eligible for screening if they have not yet been screened
    
    for i in range(sim.n): #sim.n==len(sim.people)
        for_screening[i] = True
    
    regularity = 5 #deliver the treatment every {regularity} timesteps

    if sim.t % regularity == 0:
        return np.where(for_screening==True)
    else:
        return np.where(for_screening==False)



perf = hpv.treat_num(eligibility=eligible, 
                                prob=1, 
                                product=hpv.interventions.radiation(dur=dict(dist='normal', par1=16, par2=2.)), #'ablation' 'cancerPerfectTreatment'
                                label='perf') #TODO: can model a max capacity with treat_num too, but as it stands we pretend infinite capacity
        

"""
Findings:
- administering HPVsim-standard ablation to everyone every timestep does set cancers to 0, presumably as nothing is able to get to the stage of being cancer
    - ^ but standard ablation even as frequent as every 15 timesteps results in thousands of total cancers - and every 5 timesteps still gives 549 total.

- it does look like doing my cure_cancers thing works - because although the total number of cancers is the same (well, even higher!) when running it every 50 timesteps, the deaths is only around 15k.
    - ^i am now trying to run this every 5 timesteps to see what happens... hmm still 15k deaths although I definitely am doing it every 5 timesteps.
    - still, this does seem to be working somewhat...
    - i am now doing it running it every timestep to try and catch maximal cancers. To speed things up im going from 1e4 agents to 1e3 agents but that should not affect stuff crazy i hope, its just x10 fewer agents.
        ^hmm, it did it every timestep but cancer deaths is actually more.... 
        ^THIS DOESNT MAKE SENSE TO ME, SURELY THERE MUST BE A BUG.. IF I WERE TO BE CORRECTLY REMOVING ALL CANCERS EACH TIMESTEP, HOW COULD ANY CANCER DEATHS POSSIBLY HAPPEN? (particularly as radiation seems to work in eliminating cancer deaths if we do similar ). It does look like its working a bit, jsut not as much as it surely should be.
    
-nothing seems to be happening with my cancerPerfectTreatment treatment (but it definitely is loading in because when i run it outside of anaconda i get an error while when i run it inside anaconda all is fine)

-uhhh.. looks like radiation works just fine though at extending lifespan???

- in this case, it looks like I can make an intervention a bit like how they implemented radiation which allows, of the eligible people, a certain number to be cleared of cancer and a certain proportion not to be?
    ^and then, as radiation is the only intervention whose use currently even makes 'cancer treated' go up, I can make this intervention to also add to the tally of treated cancers manually, yay!

"""

#print("trying out cancerPerfectTreatment every timestep to see what happens")

# Create and run the sim
sim = hpv.Sim(rand_seed=1,
              #verbose=-1,
              start=1970, end=2010,
              n_agents=1e3,#1e4,
              #interventions=CureCancer(),
              #interventions=[cure_cancers],
              interventions=[perf],
              #interventions=hpv.interventions.radiation(),
              )
sim.run()

df = sim.to_df()

#print(df.columns)

#print(df['cancer_deaths'])

print(df['cum_cancer_treated'])

#sim.plot()
