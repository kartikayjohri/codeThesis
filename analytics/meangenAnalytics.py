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
        if (ii > 0) and (ii < 2):
            rotor_value.append(int(line.split()[-1]))
        # Extract flow angles and other parameters of the stator
        if (ii > 6) and (ii < 8):
            stator_value.append(float(line.split()[-2]))
            stator_value.append(float(line.split()[-1]))
        if (ii > 7) and (ii < 23):
            stator_value.append(float(line.split()[-1]))
        # Extract flow angles and other parameters of the rotor
        if (ii > 26) and (ii < 28):
            rotor_value.append(float(line.split()[-2]))
            rotor_value.append(float(line.split()[-1]))
        if (ii > 27) and (ii < 43):
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

#%%

import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar
from math import sqrt,cos,tan,atan,pi,sin

#% Define function to find optimum Tr for Pr = 1
# x0 is the initial guess of a factor   
def TcritEval(x0,Fluid,zInput):
    # Z at which Tcrit is to be found for the fluid
    fluid = Fluid   
    Tcrit = CP.PropsSI('TCRIT',fluid)
    Pcrit = CP.PropsSI('PCRIT',fluid)
    Tinput = x0*Tcrit
    Z = CP.PropsSI('Z','P',Pcrit,'T',Tinput,fluid)
    error = abs((1- Z/zInput))
    return error

FLUID = 'air'
Z = 1

PHI = 0.8
PSI = 1.3
Rinput = 0.5
pressure_ratio = 2
############################ SAME FROM THIS POINT ######################
# Call TcritEval function to determine operating conditions
bnds = ((0.5, 5))
res = minimize_scalar(TcritEval, bounds=bnds,args=(FLUID,Z,), tol = 1e-3, method='bounded')
Tr = res.x
# Critical temperatures and pressures for the fluid
Tcrit = CP.PropsSI('Tcrit',FLUID)     # in K
Pcrit = CP.PropsSI('Pcrit',FLUID)/1e5 # in bars
# Inlet conditions
inletP    = Pcrit
inletT    = Tr*Tcrit
# Gas properties
Rgas  = CP.PropsSI('gas_constant',FLUID)/CP.PropsSI('M',FLUID)
Gamma = CP.PropsSI('Cpmass','T',Tr*Tcrit,'P',Pcrit*1e5,FLUID)/CP.PropsSI('Cvmass','T',Tr*Tcrit,'P',Pcrit*1e5,FLUID)
Cp = CP.PropsSI('Cpmass','T',Tr*Tcrit,'P',Pcrit*1e5,FLUID)    
# Velocity triangles
phi = PHI
psi = PSI
R   = Rinput
# Machine efficiency
eta = 0.85
# Pressure ratio over stage
PR = pressure_ratio
# Operating conditions and sizing
massFlow  = 25.0
RPM       = 9000
##############################################################################

# Inlet conditions
inletRho = CP.PropsSI('D','P',inletP*1e5,'T',inletT,FLUID)
inletH0 = CP.PropsSI('Hmass','P',inletP*1e5,'T',inletT,FLUID)
# Outlet conditions
outletP = inletP/PR
outletT = inletT/PR**(1-1/Gamma)
outletRho = CP.PropsSI('D','P',outletP*1e5,'T',outletT,FLUID)
outletH0is = CP.PropsSI('Hmass','P',outletP*1e5,'T',outletT,FLUID)
# Enthalpy drop
DeltaH0 = abs(outletH0is-inletH0)*eta
DeltaH0is = abs(outletH0is-inletH0)
# Density and volumetric ratio
DR = PR**-Gamma
VR = 1/DR # VR defined as V1/V2
# Velocity triangle
alpha1 = atan(((psi/2)+1-R)/phi) # radians
beta2 = -alpha1 # radians
beta1 = atan(tan(alpha1)-(1/phi)) # radians
alpha2 = -beta1 # radians
# Blade height to mean radius 
H_r_des = 0.2
# Aspect ratio 
AR = 2

def optimumMassFlow(x0,phi,psi,R,DeltaH0,massFlow,RPM,inletRho,H_r_des):
    
    U = sqrt(DeltaH0/psi)
    radiusDes = U/(RPM*pi/30)
    H1 = x0*massFlow/(inletRho*phi*U*2*pi*radiusDes)
    H_r_calc = H1/radiusDes
    error = abs(1 - H_r_calc/H_r_des)
    return error

bnds1 = ((1E-2, 100))
res1 = minimize_scalar(optimumMassFlow, bounds=bnds1,args=(phi,psi,R,DeltaH0,massFlow,RPM,inletRho,H_r_des,), tol = 1e-5, method='bounded')

massFlow_opt = massFlow*res1.x
U = sqrt(DeltaH0/psi)
radiusDes = U/(RPM*pi/30)
H1_new = massFlow_opt/(inletRho*phi*U*2*pi*radiusDes)
H2 = massFlow_opt/(outletRho*phi*U*2*pi*radiusDes)
H_bar = (H1_new+H2)*0.5
c_ax = H_bar/AR

H_r_opt = H1_new/radiusDes

# Velocities
Vm = (RPM*pi/30)*radiusDes*phi
V1 = Vm/cos(alpha1)
V2 = Vm/cos(alpha2)
W1 = Vm/cos(beta1)
W2 = Vm/cos(beta2)





#%%
import numpy as np

phi_range     = np.linspace(0.5,1.3,3) 
psi_range     = np.linspace(0.8,2.5,3)
rr            = 0.5

phi,psi = np.meshgrid(phi_range,psi_range)
# Velocity triangle
alpha1 = np.arctan(((psi/2)+1-rr)/phi) # radians
beta2 = -alpha1 # radians
beta1 = np.arctan(np.tan(alpha1)-(1/phi)) # radians
alpha2 = -beta1 # radians

# Solidity
Z = 0.85
sigma = 2*(1/Z)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2




