import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

from scipy import stats

#Set the following - which csv column are we plotting a time series of? note we will, by default, index everything by directory

report_type = "CIs" #"percentiles", "CIs"

column_name = "new_screens" #"cancers", "infections", "cum_screens", "new_screens", "cancer_deaths"


series_names = {"results/AlternateScreeningC/":"ASP C", #NOTE: to properly test this path i really need to make sure the effectiveness of cytology is bang on! and really also the effectiveness of all the other treatments
                "results/AlternateScreeningA/":"ASP A",
                "results/AlternateScreeningF/":"ASP F",
                "results/AlternateScreeningD/":"ASP D",
                "results/AlternateScreeningD2/":"ASP D2",
                    #^I just want to see D2 absolutely sucking, to sanity check that my code is as intended. Because surely it should suck.
                        #Unless, am i just not modelling the link between number of screenings and number of treated cancers tightly enough?

                "results/AlternateScreeningCD/":"ASP CD",

                #"results/AlternateScreeningRatio/":"ASP R", <-need to do special thing with ratio cause of lambda
                "results/NHSPathwayResults/":"Current Screening, Vaccination",
       #         "results/NHSPathwayResultsNOSCREEN/":"Vaccination, No screening",
         #       "results/NHSPathwayResultsNOVACC/":"Current Screening, No Vaccination",
              #  "results/NHSPathwayResultsNOVACCNOSCREEN/":"No Vaccination, No Screening",


  #              "results/lowvacc/AlternateScreeningC/":"ASP C, 0.3 vacc",
   #             "results/lowvacc/AlternateScreeningA/":"ASP A, 0.3 vacc",
    #            "results/lowvacc/AlternateScreeningF/":"ASP F, 0.3 vacc",
     #           "results/lowvacc/AlternateScreeningD/":"ASP D, 0.3 vacc",
           #     "results/lowvacc/NHSPathwayResults/":"Current Screening, 0.3 vacc",
                } #{folder in which to find the data for the series: series label}

#I think the 0.3vacc values are not annual so the vaccination is actually much higher than I intend??
#quite possibly they are actually. I got around 60% vaccinated full pop in 2080 it seems with v=0.8, looking at the excel, and this is not crazy high so probs correct? Though I suppose I would need to see how many of a certain age there are in 2080?
    #^ actually, I can properly verify this by running a single sim for a really long time, and then all people alive will have had a chance to get vaccinated - then if the proportion matches what I want it to be, I am doing it right!
# (also, v=0.3 was about 3/8ths number of vaccinated people as v=0.8, which is good!)

#Reading data
seriess_by_dirs = {x:[] for x in series_names.keys()}
for dir in seriess_by_dirs.keys(): #loop through all results directories
    print(f"Reading {dir}")
    files_in_dir = os.listdir(dir) #find all files in the relevant results directiory


    for file in files_in_dir: #for each found file, load its relevant series
        if file[-7:] != ".pickle":
            df = pd.read_excel(dir+file)
            series = list(df[column_name])
            seriess_by_dirs[dir].append(series)


           


"""
#Do a basic plot of the data
colours = ["blue", "orange", "green", "aqua", "gray", "purple"]
xs = np.arange(1970,2081,1)

i=0
for dir in seriess_by_dirs.keys():
    j=0
    for series in seriess_by_dirs[dir]:
        plt.plot(xs,series, label=series_names[dir] if j==0 else None, color=colours[i])
        j+=1
    i+=1


plt.legend()
plt.show()
plt.cla() 
"""

#Do a shaded plot of the central 50% of the data
print("Plotting...")
colours = ["blue", "orange", "green", "aqua", "gray", "purple", "magenta", "black", "red"]
xs = np.arange(1970,2081,1)

i=0
for dir in seriess_by_dirs.keys():
    seriess = seriess_by_dirs[dir]

    if report_type == "CIs":
        #mean and 95% CI on either side
        averages = np.mean(seriess,axis=0)
        sems = stats.sem(seriess, axis=0)
        margins_of_error = sems * stats.t.ppf((1+0.95) / 2, len(seriess) - 1) #seriess is square, so no need to do the multiplication factor differently across timepoints
        lowers = averages-margins_of_error 
        uppers = averages+margins_of_error
    elif report_type == "percentiles":
        #LQ, Median, UQ
        lowers = np.percentile(seriess, 25, axis=0)
        averages = np.percentile(seriess, 50, axis=0)#
        uppers = np.percentile(seriess, 75, axis=0)
    
    
    plt.plot(xs, averages, label=f"{series_names[dir]}, n={len(seriess)}", color=colours[i])
  #  plt.fill_between(xs, lowers, uppers, color=colours[i], alpha=0.5)

    i+=1

plt.vlines(2008, 0,6400, color="red", linestyles="dashed")

plt.vlines(2020, 0,6400, color="red", linestyles="dashed")

plt.xlim((2000,2081))

plt.xlabel("Year")
plt.ylabel(f"{column_name}")

plt.legend()
plt.show()