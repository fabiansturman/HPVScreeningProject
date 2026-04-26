import numpy as np
import sciris as sc
import pandas as pd

from colorScheme import natureColorScheme

import matplotlib.pyplot as plt


colours = [natureColorScheme['reds'][2], natureColorScheme['blues'][2]]
hatches = ['/////','\\\\\\\\\\']
alpha = 0.2

fig, ax = plt.subplots(1,1)

plt.hlines(90, 0,11, colors='green', linestyles = 'dashed', label='NHS Target Vaccination Level: 90%')
plt.hlines(80, 0,11, colors='black', linestyles = 'dashed', label='Middle Vaccination Scenario 80%')

plt.hlines(60, 0, 11, colors='red', linestyles = 'dashed', label='Low Vaccination Scenario: 60% ')


boxplot = ax.bar([0.35,1.35,2.35,3.35,4.35,],
                     [90.4,88.8,89.1,88.9,88.9],
                     width=0.5,
                     color=colours[0],
                     alpha=alpha,
                     edgecolor = 'none',
                     hatch=hatches[0])
boxplot = ax.bar([0.35,1.35,2.35,3.35,4.35,],
                     [90.4,88.8,89.1,88.9,88.9],
                     width=0.5,
                     color='none',
                     edgecolor = colours[0],
                     hatch=hatches[0],
                     label='Girls')

boxplot = ax.bar([5,6,7,8,9,10 ],
                     [88.5, 83.2, 76.7,75.5,75.3,71.7 ],
                     width = 0.3, 
                     align='edge',
                     color=colours[0],
                     alpha=alpha,
                     edgecolor='none',
                     hatch=hatches[0])
boxplot = ax.bar([5,6,7,8,9, 10],
                     [88.5, 83.2, 76.7, 75.5,75.3,71.7],
                     width = 0.3, 
                     align='edge',
                     color='none',
                     edgecolor=colours[0],
                     hatch=hatches[0])


boxplot = ax.bar([ 5.4, 6.4, 7.4, 8.4, 9.4,10.4],
                     [81.5, 78.6, 71.2, 70.5, 70.5 ,67.0 ],
                     width = 0.3, 
                     align='edge',
                     color=colours[1],
                     alpha=alpha,
                     edgecolor = 'none',
                     hatch=hatches[1])

boxplot = ax.bar([ 5.4, 6.4, 7.4, 8.4, 9.4, 10.4],
                     [81.5, 78.6, 71.2, 70.5, 70.5 ,67.0],
                     width = 0.3, 
                     align='edge',
                     color='none',
                     edgecolor = colours[1],
                     hatch=hatches[1],
                     label='Boys')

plt.xticks([0.35,1.35,2.35,3.35,4.35,5.35,6.35,7.35,8.35,9.35, 10.35], ['2014', '2015', '2016','2017', '2018', '2019', '2020', '2021', '2022','2023', '2024'])
plt.yticks([0,10,20,30,40,50,60,70,80,90])

def round_str_to_str(n):
    if int(n[-1])>=5:
        return str(int(str(d[i])[:2])+1)
    else:
        return str(d[i])[:2]

d = [90.4,88.8,89.1,88.9,88.9]
for i in range(5):
    if d[i]!=88.9:
        plt.text(0.27+i, d[i]+1.5, round_str_to_str(str(d[i])))
    else:
        plt.text(0.27+i, d[i]+2, round_str_to_str(str(d[i])))

d = [88.5, 83.2, 76.7, 75.5, 75.3, 71.7 ]
for i in range(6):
    if i==0:
        plt.text(5.07+i, d[i]+2, round_str_to_str(str(d[i])))
    elif i==2:
        plt.text(5.07+i, d[i]+4, round_str_to_str(str(d[i])))
    else:
        plt.text(5.07+i, d[i]+1.5, round_str_to_str(str(d[i])))

d = [81.5, 78.6, 71.2,  70.5, 70.5 ,67.0 ]
for i in range(6):
    if i!=1:
        plt.text(5.47+i, d[i]+1.5, round_str_to_str(str(d[i])))
    else:
        plt.text(5.47+i, d[i]+2.5, round_str_to_str(str(d[i])))

"""
d = [90.4,88.8,89.1,88.9,88.9]
for i in range(5):
    plt.text(0.27+i, d[i]+1.5, str(d[i])[:2])

d = [88.5, 83.2, 76.7, 74.1, 72.9 ]
for i in range(5):
    plt.text(5.07+i, d[i]+1.5, str(d[i])[:2])

d = [81.5, 78.6, 71.2, 68.5, 67.7 ]
for i in range(5):
    plt.text(5.47+i, d[i]+1.5, str(d[i])[:2])
"""

plt.xlabel('Calander year in which students began academic Year 8')
plt.ylabel('Rate of vaccination uptake (%) in cohort, between Years 8-10')

#plt.legend()

plt.show()

