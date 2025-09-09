"""This script makes the figure for the Ratio screening algorithm, where a parameter
     lambda controls how much more frequently the unvaccinated have regular screenings 
     compared to the vaccinated, and then the vaccinated have less screenings to 
     make total screening burden around constant"""

import pickle
import os
import matplotlib.pyplot as plt

#Below, set the directory from which we load results and the lambdas of the results we want to load
dir = "results/AlternateScreeningRatio/"
lambdas = [0.5,1.0,1.5,2.0,2.5,3.0]

#TODO: add a horizontal line as a baseline, equivalent to lambda=1, for the current strategy. Or maybe just add it as a bar idk 
    # ^though, the baseline may not be exactly equivalent to lambda=1 if my v is wrong (and as it is static, it will be) - but i suppose the key result here is as long as i get the sorta curve shape i want im good

#Reading data (dictionary {lambda: data}, where each data is a list of dictionaries, where each dictionary is the output of a simulation run)
results_by_lambda = {l:[] for l in lambdas}

files_in_dir = os.listdir(dir)
files_in_dir = [filename for filename in files_in_dir if '.pickle'==filename[-7:]]

for l in results_by_lambda.keys():
    l_str = str(l)
    relevant_files = [filename for filename in files_in_dir if l_str==filename[-len(l_str)-7:-7]]
    
    for filename in relevant_files:
        with open(dir+filename, 'rb') as file:
            data_dict = pickle.load(file)
            results_by_lambda[l].append(data_dict)

#TODO: the script below gives the general gist, but I really want it to be a curve

#Only counting people who got cancer after timepoint T1 and before timepoint T2, calculate and plot number of (cancers, cancers that get put through colposcopy, ratio of the two), as lists indexed by dir
start_year = 2025
end_year = 2080

T1 = int((start_year-1970)/0.25)
T2 = int((end_year-1970)/0.25)
    #^inclusive bounds

    

cancer_plotdata_TRUNCATED = []
cancer_found_plotdata_TRUNCATED = []
prop_cancer_found_plotdata_TRUNCATED = []
labels = []
for dir in results_by_lambda.keys():
    cancer = []
    cancer_found = []
    prop_cancer_found = []

    for result in results_by_lambda[dir]:
        #We only count patients who first get cancer/first have a colposcopy between timepoints T1 and T2
        #TODO: this is not ideal! I should actually probs record ALL TIMEPOINTS where a person has a colposcopy with cancer, and ALL TIMEPOINTS where they have cancer, as e.g. if the cancer starts just beore T1 but persists into the interval of interest, we should still count it surely!
            #^I think this may be seriously messing with my results, perhaps!
            
        people_with_cancer_in_timerange = 0
        screened_people_with_cancer_in_timerange = 0

        for pid in set(result['with_cancer']):
            tp, scale = result['time_with_first_cancer'][pid]
            if T1 <= tp <= T2:
                people_with_cancer_in_timerange += scale

        if 'colposcopy_eligible_with_cancer' in result.keys(): #only in keys if colposcopies were administered 
            for pid in set(result['colposcopy_eligible_with_cancer']):
                tp, scale = result['time_with_first_cancer_colposcopy'][pid]
                if T1 <= tp <= T2:
                    screened_people_with_cancer_in_timerange += scale
                
        cancer.append(people_with_cancer_in_timerange)
        cancer_found.append(screened_people_with_cancer_in_timerange)
        prop_cancer_found.append(screened_people_with_cancer_in_timerange/people_with_cancer_in_timerange)

    cancer_plotdata_TRUNCATED.append(cancer)
    cancer_found_plotdata_TRUNCATED.append(cancer_found)
    prop_cancer_found_plotdata_TRUNCATED.append(prop_cancer_found)

    labels.append(f"{dir}, n={len(cancer)}")

plt.boxplot(cancer_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.xlabel("lambda")
plt.show()

plt.boxplot(cancer_found_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers found through colposcopy, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.xlabel("lambda")
plt.show()

plt.boxplot(prop_cancer_found_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers found thorugh colposcopy / cancers, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.xlabel("lambda")
plt.show()