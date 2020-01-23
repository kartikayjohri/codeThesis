#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 18:56:31 2020

@author: kjohri
"""
#%% Import modules

import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar
from math import sqrt,cos,tan,atan,pi,sin

#%% Define functions 

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

#% Calculate opitimum mass flow
def optimumMassFlow(x0,phi,psi,R,DeltaH0,massFlow,RPM,inletRho,H_r_des):
    
    U = sqrt(DeltaH0/psi)
    radiusDes = U/(RPM*pi/30)
    H1 = x0*massFlow/(inletRho*phi*U*2*pi*radiusDes)
    H_r_calc = H1/radiusDes
    error_mflow = abs(1 - H_r_calc/H_r_des)
    return error_mflow

#% Zweifel Criteria
def Zweifel(Zcr,alpha1,alpha2):
    # Solidity
    sigma = 2*(1/Zcr)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2
    return sigma

# Determine preliminary inputs for meangen
def meangenInputs(PHI,PSI,Rinput,Z,pressure_ratio,massFlow,RPM,H_r_des,AR):
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
    massFlow  = massFlow
    RPM       = RPM
##############################################################################    
    # Inlet conditions
    inletRho = CP.PropsSI('D','P',inletP*1e5,'T',inletT,FLUID)
    inletH0  = CP.PropsSI('Hmass','P',inletP*1e5,'T',inletT,FLUID)
    # Outlet conditions
    outletP    = inletP/PR
    outletT    = inletT/PR**(1-1/Gamma)
    outletRho  = CP.PropsSI('D','P',outletP*1e5,'T',outletT,FLUID)
    outletH0is = CP.PropsSI('Hmass','P',outletP*1e5,'T',outletT,FLUID)
    # Enthalpy drop
    DeltaH0   = abs(outletH0is-inletH0)*eta
    DeltaH0is = abs(outletH0is-inletH0)
    # Density and volumetric ratio
    DR = PR**-Gamma
    VR = 1/DR # VR defined as V1/V2
    # Velocity triangle
    alpha1 = atan(((psi/2)+1-R)/phi) # radians
    beta2  = -alpha1 # radians
    beta1  = atan(tan(alpha1)-(1/phi)) # radians
    alpha2 = -beta1 # radians
    
    # Iterate over massflow to determine the optimum mass flow for requirements
    bnds1 = ((1E-3, 1000))
    res1 = minimize_scalar(optimumMassFlow, bounds=bnds1,args=(phi,psi,R,DeltaH0,massFlow,RPM,inletRho,H_r_des,), tol = 1e-5, method='bounded')
    
    # Calculate U, radiusDes
    U = sqrt(DeltaH0/psi)
    radiusDes = U/(RPM*pi/30)
    # Determine values with optimum mass flow
    massFlow_opt = massFlow*res1.x
    # Rotor blade height at inlet, outlet and average
    H1 = massFlow_opt/(inletRho*phi*U*2*pi*radiusDes)
    H2 = massFlow_opt/(outletRho*phi*U*2*pi*radiusDes)
    H_bar = (H1+H2)*0.5
    # Calculate axial chord
    c_ax = H_bar/AR
    H_r_opt = H1/radiusDes
    # Velocities
    Vm = (RPM*pi/30)*radiusDes*phi
    V1 = Vm/cos(alpha1)
    V2 = Vm/cos(alpha2)
    W1 = Vm/cos(beta1)
    W2 = Vm/cos(beta2)
    
    return alpha1,alpha2,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_opt,c_ax,RPM,DeltaH0,eta

#%% Write inputs and call functions 
# Fluid
FLUID = 'air'
# Compressibility
Z = 1
# Stage coefficients
PHI = 0.8
PSI = 1.3
Rinput = 0.5
# Design pressure ratio
pressure_ratio = 2
# Machine design
RPM = 9000
massFlow = 25
# Geometrical parameters
# Blade height to design radius
H_r_des = 0.2
# Aspect ratio
AR = 2
# Zweifel coefficient
Zcr = 0.9

# Call functions and assign values
alpha1,alpha2,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta = meangenInputs(PHI,PSI,Rinput,Z,pressure_ratio,massFlow,RPM,H_r_des,AR)
solidity = Zweifel(Zcr,alpha1,alpha2)

from inputMeangen_old import inputMeangen
inputMeangen(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta)
