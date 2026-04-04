
## Result Replication

Run any of the "*Gen.py" files to run a scenario (these load in parameters from CalibrationResults.py, which are lists of the usually-10 best parameters from different calibrations.).

This generates one run of a HPVsim instance, with each instance's run information saved as a pickle file in sim_logs/. Each such file contains a detailed timeseries for the simulation, storing information such as cancer/infection/cin state of each agent at each point,and  interventions administered and intervention outcomes at any given point.

These 'raw result' pickle files can be processed using agentLogging_READINGLOGS_withPairingAbility.py, which appends to (or creates) a 'processed results' pickle file, containing timeseries of the population-level quantities we report (cancers over time, cins over time, etc).

These can then be read by the analysis scripts, such as pairedTTestRig.py.

If you wish to perform alternate analysis on the agent-wise timeseries to get different population-level results, or for more detail on how simulation results are logged in the simulation-wise 'raw result' files, agentLogging_READINGLOGS_withPairingAbility.py is a good starting point.
