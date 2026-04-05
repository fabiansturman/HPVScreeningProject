## Installation

1. Make a local copy/clone of this repo.
2. Make a fresh Python virtual environment, version 3.11.15, and enter this virtual environment for remainder of all following steps.
3. Install required modules - versioning is important. See requirements.txt for a `pip freeze` of the packages installed on the Python venv with which we generated our results.
5. Inside the repo root directory (i.e. the directory within which we can access 'hpvsim', 'project_documents', 'project_modelling', etc), run `pip install -e .`. This installs our version of HPVsim into the repo, in editable mode (such that any changes you subsequently make within the 'hpvsim' directory are reflected when you run HPVsim, automatically).
6. Run the following
```
python
import hpvsim as hpv
```
  You should observe the output 'Loading HPVsim adapted for England Cervical Screening Project.'. If this appears, installation is successful. If not, consult previous steps.

## Disclaimer

This code allows replication of modelling results from a research project into Cervical Screening in England, using the HPVsim model as developed by the Institute of Disease Modelling (IDM), the Burnet Institute, and other collaborators. Its underlying model is HPVsim 2.2.4, with adaptations to align with modelling needs in England. See https://github.com/StarSimHub/HPVSim for the latest version of the base model. 



## Result Replication

Run any of the "*Gen.py" files to run a scenario (these load in parameters from CalibrationResults.py, which are lists of the usually-10 best parameters from different calibrations.).

This generates one run of a HPVsim instance, with each instance's run information saved as a pickle file in sim_logs/. Each such file contains a detailed timeseries for the simulation, storing information such as cancer/infection/cin state of each agent at each point,and  interventions administered and intervention outcomes at any given point.

These 'raw result' pickle files can be processed using agentLogging_READINGLOGS_withPairingAbility.py, which appends to (or creates) a 'processed results' pickle file, containing timeseries of the population-level quantities we report (cancers over time, cins over time, etc).

These can then be read by the analysis scripts, such as pairedTTestRig.py.

If you wish to perform alternate analysis on the agent-wise timeseries to get different population-level results, or for more detail on how simulation results are logged in the simulation-wise 'raw result' files, agentLogging_READINGLOGS_withPairingAbility.py is a good starting point.
