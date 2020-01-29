# -*- coding: utf-8 -*-
"""
FUNCTION: To look for anomalies, or trends in the blade generation data
@author: kjohri
"""
# import modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Select fluid being analysed
fluid = 'air'
folders = sorted(os.listdir("/home/kjohri/Documents/bladeGeneration/"+fluid))

# Set column names based on data present in the meangen.out file per blade row
cols_rotor = ['nBlades','beta1','beta2','Vx','M_in','M_exit','rho_exit','P_exit','Pt_in','Pt_exit','Pt_rel_in','T_exit','Tt_exit','r_tip','b_in','c_ax','AR','pitch']
cols_stator = ['nBlades','alpha1','alpha2','Vx','M_in','M_exit','rho_exit','P_exit','Pt_in','Pt_exit','Pt_rel_in','T_exit','Tt_exit','r_tip','b_in','c_ax','AR','pitch']

# Set index of the dataframe to that of the foldernames and create dataframe
fols = np.array(folders,dtype = int)
df_rotor = pd.DataFrame(columns = cols_rotor, index = fols)
df_stator = pd.DataFrame(columns = cols_stator, index = fols)

# Extract values from meangen.out per folder
folder_counter = 0
for folder in folders:  
    os.chdir("/home/kjohri/Documents/bladeGeneration/air/"+folder)
    f = open('meandesign.out')
    lines = f.readlines()
    f.close()
    ii = 0
    stator_value = []    
    rotor_value = []
    for line in lines:
        # Extract blade numbers
        if ii == 0:
            stator_value.append(int(line.split()[-1]))
        elif (ii > 0) and (ii < 2):
            rotor_value.append(int(line.split()[-1]))
        # Extract flow angles and other parameters of the stator
        elif (ii > 6) and (ii < 8):
            stator_value.append(float(line.split()[-2]))
            stator_value.append(float(line.split()[-1]))
        elif (ii > 7) and (ii < 23):
            stator_value.append(float(line.split()[-1]))
        # Extract flow angles and other parameters of the rotor
        elif (ii > 26) and (ii < 28):
            rotor_value.append(float(line.split()[-2]))
            rotor_value.append(float(line.split()[-1]))
        elif (ii > 27) and (ii < 43):
            rotor_value.append(float(line.split()[-1]))
        ii = ii+ 1
    # Append values to the dataframes
    df_stator.loc[int(folder)] = stator_value
    df_rotor.loc[int(folder)] = rotor_value
    
    folder_counter = folder_counter + 1

#%% Plotting 
    
df_rotor.reset_index().plot(x = 'index', y = 'beta1', kind ='line')
df_stator.reset_index().plot(x = 'index', y = 'alpha1', kind ='bar')

# NOTE: More needs to be looked into the plotting capabilities.

