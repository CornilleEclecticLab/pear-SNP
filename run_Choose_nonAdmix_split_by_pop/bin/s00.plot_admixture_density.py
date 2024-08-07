#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File     : s00.plot_admix_density.py
# @Version  : 1.0.0
# @Author   : NIE Yuqi
# @Email    : nieyuqi.cn@gmail.com
# @Time(CET): 2023/11/27 15:15:19
# @Description:
#     

import datetime
import sys
import textwrap
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

start_time = datetime.datetime.now()
print(f'{" Start ":=^79}')



version = "1.0.0"
script_basename = os.path.splitext(os.path.basename(sys.argv[0]))[0]
script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
bin_dir = script_path
work_dir = os.path.dirname(script_path)
input_dir = os.path.join(work_dir,'input')
output_dir = os.path.join(work_dir,'output',script_basename)
sub_script_dir = os.path.join(bin_dir,script_basename)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to wrap text to 79 characters
def wrap(text, width=79):
    return textwrap.fill(text, width=width, subsequent_indent=' '*4)


Qfile="data/K.15.MajorCluster.Q"
matrix = np.loadtxt(Qfile)

# Calculate the maximum coefficient values for each individual
max_coefficients = np.max(matrix, axis=1)

# Create a kernel density estimate (KDE) for the maximum coefficient values
kde = gaussian_kde(max_coefficients)

# Generate x values for the density plot
x_values = np.linspace(0, 1, 100)

# Add a vertical line at x=0.8
plt.axvline(x=0.8, color='red', linestyle='--', label='Threshold')


# Plot the density
plt.plot(x_values, kde(x_values), label='Density')
plt.title('Density Plot of Maximum Membership Coefficient Values')
plt.xlabel('Maximum Membership Coefficient Value')
plt.ylabel('Number of Individuals')
plt.legend()
plt.show()



end_time = datetime.datetime.now()
print('')
print(' END '.center(79,'='))
print(str(end_time-start_time).center(79))