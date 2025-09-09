import pickle
import os
import matplotlib.pyplot as plt

#track proportion of hpv cases which are diagnosed through a positive screen - #positive tests/#positives


#it may be interesting/useful to literally track at each timepoint who has had each treatment, from which i can calculate things like box plots/means for how many screens each woman gets, tracking the time it takes between onset of  cancer for a woman and it being caught by screening (in a better way than I currently have)
#yeah maybe even better would be to store lists for each agent keeping track of, at each year - cancer status (or just, state), vaccination status, and which interventions they have undergone

#Define the directories which hold the data of interest
names_by_dir = {"results/AlternateScreeningC/":"C", #NOTE: to properly test this path i really need to make sure the effectiveness of cytology is bang on! and really also the effectiveness of all the other treatments
                "results/AlternateScreeningA/":"A",
                "results/AlternateScreeningF/":"F",
                "results/AlternateScreeningD/":"D",
                "results/AlternateScreeningD2/":"D2",
                    #^I just want to see D2 absolutely sucking, to sanity check that my code is as intended. Because surely it should suck.
                        #Unless, am i just not modelling the link between number of screenings and number of treated cancers tightly enough?

                
                "results/AlternateScreeningCD/":"ASP CD",

                #"results/AlternateScreeningRatio/":"ASP R", <-need to do special thing with ratio cause of lambda
                "results/NHSPathwayResults/":"cnt.",
            #    "results/NHSPathwayResultsNOSCREEN/":"v only",
             #   "results/NHSPathwayResultsNOVACC/":"scrn only",
              #  "results/NHSPathwayResultsNOVACCNOSCREEN/":"no int",

    #            "results/lowvacc/AlternateScreeningC/":"ASP C, 0.3 vacc",
     #           "results/lowvacc/AlternateScreeningA/":"ASP A, 0.3 vacc",
      #          "results/lowvacc/AlternateScreeningF/":"ASP F, 0.3 vacc",
       #         "results/lowvacc/AlternateScreeningD/":"ASP D, 0.3 vacc",
        #        "results/lowvacc/NHSPathwayResults/":"cnt., 0.3 vacc",
                } 



#Reading data (each set of data is a list of dictionarys, where each dictionary is the output of a simulation run)
results_by_dir = {x:[] for x in names_by_dir.keys()}
for dir in names_by_dir.keys(): 
    files_in_dir = os.listdir(dir) #find all files in the relevant results directiory


    for file in files_in_dir: #for each found file, load its relevant series
        if file[-7:] == ".pickle": #TODO: probably should fix this  naming conflict, but also works fine
            with open(dir+file, 'rb') as file:
                data_dict = pickle.load(file)
                results_by_dir[dir].append(data_dict)


#Only counting people who got cancer after timepoint T1 and before timepoint T2, calculate and plot number of (cancers, cancers that get put through colposcopy, ratio of the two), as lists indexed by dir
start_year = 2026
end_year = 2080

T1 = int((start_year-1970)/0.25)
T2 = int((end_year-1970)/0.25)
    #^inclusive bounds

    

cancer_plotdata_TRUNCATED = []
cancer_found_plotdata_TRUNCATED = []
prop_cancer_found_plotdata_TRUNCATED = []
labels = []
for dir in results_by_dir.keys():
    cancer = []
    cancer_found = []
    prop_cancer_found = []

    for result in results_by_dir[dir]:
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

    labels.append(names_by_dir[dir])

plt.boxplot(cancer_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.show()

plt.boxplot(cancer_found_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers found through colposcopy, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.show()

plt.boxplot(prop_cancer_found_plotdata_TRUNCATED, labels=labels)
plt.ylabel(f"cancers found thorugh colposcopy / cancers, {int(T1*0.25+1970)} - {int(T2*0.25+1970)}")
plt.show()




"""
pickle.dump({'final_timepoint':algo.t,
                        'with_cancer':with_cancer,
                        'colposcopy_eligible_with_cancer':colposcopy_eligible_with_cancer,
                        'time_with_first_cancer': time_with_first_cancer,
                        'time_with_first_cancer_colposcopy':time_with_first_cancer_colposcopy}, file=file)
                        """