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
folders = sorted(os.listdir("/home/kjohri/Documents/codeThesis/"+fluid))

# Set column names based on data present in the meangen.out file per blade row
cols_rotor = ['phi','psi','nBlades','beta1' ,'beta2' ,'Vx','M_in','M_exit','P_in','P_exit','Tt_in','Tt_exit','Ht_in','Ht_exit','rho_in','rho_exit','Pt_in','Pt_exit','Pt_in_rel','Pt_exit_rel','r_tip','b_in','b_out','sigma_ax','AR','c_ax','g_ax','H1/Rm','TETc','Rm']
cols_stator= ['phi','psi','nBlades','alpha1','alpha2','Vx','M_in','M_exit','P_in','P_exit','Tt_in','Tt_exit','Ht_in','Ht_exit','rho_in','rho_exit','Pt_in','Pt_exit','r_tip','b_in','b_out','sigma_ax','AR','c_ax','g_ax','H1/Rm']

# Set index of the dataframe to that of the foldernames and create dataframe
fols = np.array(folders,dtype = int)
df_rotor = pd.DataFrame(columns = cols_rotor, index = fols)
df_stator = pd.DataFrame(columns = cols_stator, index = fols)

# Extract values from meangen.out per folder
folder_counter = 0
for folder in folders:  
    os.chdir("/home/kjohri/Documents/codeThesis/air/"+folder)
    f = open('meandesign.out')
    lines = f.readlines()
    f.close()
    ii = 0
    Nblades = []
    stator_angles = []
    rotor_angles = []
    stator_values = []
    rotor_values = []
    # Store phi
    stator_values.append(float(folder[0:4])/1000)
    stator_values.append(float(folder[4:])/1000)
    #Store psi
    rotor_values.append(float(folder[0:4])/1000)
    rotor_values.append(float(folder[4:])/1000)
    # Read meangen.out files
    for line in lines:
        # Extract blade numbers
        if ii == 0:
            stator_values.append(int(line.split()[-1]))
        elif (ii > 0) and (ii < 2):
            rotor_values.append(int(line.split()[-1]))            
        elif (ii > 6) and (ii < 8):
            stator_values.append(float(line.split()[-2]))
            stator_values.append(float(line.split()[-1]))
        elif (ii > 7) and (ii < 29):
            stator_values.append(float(line.split()[-1]))
        # Extract flow angles and other parameters of the rotor
        elif (ii > 32) and (ii < 34):
            rotor_values.append(float(line.split()[-2]))
            rotor_values.append(float(line.split()[-1]))
        elif (ii > 33) and (ii < 59):
            rotor_values.append(float(line.split()[-1]))
        ii = ii+ 1
    # Append values to the dataframes
    df_stator.loc[int(folder)] = stator_values
    df_rotor.loc[int(folder)] = rotor_values
    
    folder_counter = folder_counter + 1

#%% Plotting 
#    
df_rotor.reset_index().plot(x = 'index', y = 'sigma_ax', kind ='bar',title = 'Zcr = 0.95')
df_rotor.reset_index().plot(x = 'index', y = 'M_exit', kind ='bar')
df_rotor.reset_index().plot(x = 'index', y = 'M_in', kind ='bar')
#df_stator.reset_index().plot(x = 'index', y = 'sigma_ax', kind ='bar')

# NOTE: More needs to be looked into the plotting capabilities.

