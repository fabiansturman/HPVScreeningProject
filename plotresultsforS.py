"""This script makes the figure for the S screening algorithm, where parameters
     lambda, my controls how much more frequently the unvaccinated/vaccinated have regular screenings 
     after a negative HPV test"""

import pickle
import os
import matplotlib.pyplot as plt

#Below, set the directory from which we load results and the lambdas of the results we want to load
dir = "results/AlternateScreeningS/"

#TODO: complete this to plot the desired quantity(ies) over a grid with colour represetning what we want. 
#       ^Perhaps have green shades for when we achieve our 2040 cervical cancer elimination goal? 