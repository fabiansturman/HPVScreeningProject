import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load the CSV file
file_path = '/Users/sophiekeenan/Desktop/Calibration/sensdata.csv'
data = pd.read_csv(file_path)

# Remove commas and convert the data to floats
data_cleaned = data.drop(columns=data.columns[0])  # Drop the first column
data_cleaned = data_cleaned.apply(lambda x: x.str.replace(',', '').astype(float))

# Preparing X, Y, Z data for the 3D plot, swapping X and Y
y = np.arange(data_cleaned.shape[1])  # Use the number of columns for Y
x = np.arange(data_cleaned.shape[0])  # Use the number of rows for X
X, Y = np.meshgrid(y, x)  # Swap x and y in meshgrid
Z = data_cleaned.values

# Create a 3D surface plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Generate the surface plot
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

# Add a color bar which maps values to colors
fig.colorbar(surf, shrink=0.5, aspect=5)

# Update Labels to reflect the swap
ax.set_xlabel('Vaccine Efficacy (2-Dose Schedule)')  # Now represents the original Y data
ax.set_ylabel('Vaccine Uptake(% of Eligible Population)')# Now represents the original X data
ax.set_zlabel('Total HPV Infections (million)')
ax.set_title('3D Surface Plot Heatmap with Swapped Axes')

# Rotate the graph 90 degrees clockwise about the Z-axis
ax.view_init(elev=10, azim=-20)  # Adjust the azimuth to -90 degrees

plt.show()
