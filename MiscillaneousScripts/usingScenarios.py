#AS scenarios don't have as much flexibility as multisims, though, the HPVsim tutorial reccomends only using Scenarios for quick investigations

import hpvsim as hpv
import pandas as pd

if __name__=="__main__":
  # Set base parameters -- these will be shared across all scenarios
  basepars = {'n_agents':10e3}

  # Configure the settings for each scenario
  scenarios = {'baseline': {
                'name':'Baseline',
                'pars': {}
                },
              'high_rel_trans': {
                'name':'High rel trans (0.75)',
                'pars': {
                    'beta': 0.75,
                    }
                },
              'low_rel_trans': {
                'name':'Low rel trans(0.25)',
                'pars': {
                    'beta': 0.25,
                    }
                },
              }

  # Run and plot the scenarios
  scens = hpv.Scenarios(basepars=basepars, scenarios=scenarios)
  scens.run()
  scens.plot();