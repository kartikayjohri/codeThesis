"""
FUNCTION: To perform preliminary calculations for Meangen and then optimise 
          inputs for meangen to ensure similarity across design space.
@author: Kartikay Johri
"""

# Import modules
import os
import time
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize
from math import sqrt,cos,tan,atan,pi
from scipy.optimize import minimize_scalar

# Setting the present directory to be sure
os.chdir('/home/kjohri/Documents/codeThesis/preliminaryCalculations/')
from inputMeangen_new import inputMeangen_new



# Measure run time
start_time = time.time()
#%% Define functions 

#% Define function to find optimum Tr for Pr = 1
def tcritEval(x0,Fluid,zInput):      
    Tcrit = CP.PropsSI('TCRIT',Fluid)
    Pcrit = CP.PropsSI('PCRIT',Fluid)
    Tinput = x0*Tcrit
    Z = CP.PropsSI('Z','P',Pcrit,'T',Tinput,Fluid)
    error = abs((1- Z/zInput))
    return error

#% Calculate opitimum mass flow
def optimumMassFlow(x0,phi,psi,R,deltaH0,massFlow,RPM,inletRho,outletRho,Hr_des):
    # Calculate U, radiusDes
    U = sqrt(deltaH0/psi)
    radiusDes = U/(RPM*pi/30)
    # Rotor blade height at inlet, outlet and average
    H1_stator = x0*massFlow/(inletRho*phi*U*2*pi*radiusDes)
    H2_rotor = x0*massFlow/(outletRho*phi*U*2*pi*radiusDes)
    H2_stator = (H1_stator+H2_rotor)*0.5
    Hr_calc = H2_stator/radiusDes
    error_mflow = abs(1 - Hr_calc/Hr_des)
    return error_mflow

# Determine preliminary inputs for meangen
def meangenInputs(PHI,PSI,Rinput,Z,pressure_ratio,massFlow,RPM,Hr_des,AR):
    # Call tcritEval function to determine operating conditions
    bnds = ((0.5, 5))
    res = minimize_scalar(tcritEval, bounds=bnds,args=(FLUID,Z,), tol = 1e-3, method='bounded')
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

    # Stage coefficients
    phi = PHI
    psi = PSI
    R   = Rinput
    # Velocity triangle
    alpha1 = atan(((psi/2)+1-R)/phi)    # radians
    beta1  = atan(tan(alpha1)-(1/phi))  # radians
    alpha2 = -beta1                     # radians
    beta2  = -alpha1                    # radians

    # Pressure ratio over stage
    PR = pressure_ratio    
    # Machine efficiency
    eta = 0.85

    # Inlet thermodynamic conditions
    inletRho = CP.PropsSI('D','P',inletP*1e5,'T',inletT,FLUID)
    inletH0  = CP.PropsSI('Hmass','P',inletP*1e5,'T',inletT,FLUID)
    # Outlet thermodynamic conditions
    outletP    = inletP/PR
    outletT    = inletT/PR**(1-1/Gamma)
    outletRho  = CP.PropsSI('D','P',outletP*1e5,'T',outletT,FLUID)
    outletH0is = CP.PropsSI('Hmass','P',outletP*1e5,'T',outletT,FLUID)
    # Enthalpy drop across stage
    deltaH0is = abs(inletH0-outletH0is)
    deltaH0   = abs(inletH0-outletH0is)*eta

    # Iterate over massflow to determine the optimum mass flow for requirements
    bnds1 = ((1E-3, 1000))
    res1 = minimize_scalar(optimumMassFlow, bounds=bnds1,args=(phi,psi,R,deltaH0,massFlow,RPM,inletRho,outletRho,Hr_des,), tol = 1e-5, method='bounded')
    
    # Calculate U, radiusDes
    U = sqrt(deltaH0/psi)
    radiusDes = U/(RPM*pi/30)
    # Determine values with optimum mass flow
    massFlow_opt = massFlow*res1.x
    # Rotor blade height at inlet, outlet and average
    H1_stator = massFlow_opt/(inletRho*phi*U*2*pi*radiusDes)
    H2_rotor = massFlow_opt/(outletRho*phi*U*2*pi*radiusDes)
    H2_stator = (H1_stator+H2_rotor)*0.5
    H1_rotor = H2_stator
    H_bar_stator = (H1_stator+H2_stator)*0.5
    H_bar_rotor = (H1_rotor+H2_rotor)*0.5
    
    # Calculate axial chord
    c_ax_stator = H_bar_stator/AR[0]
    c_ax_rotor = H_bar_rotor/AR[1]
    c_ax = [c_ax_stator,c_ax_rotor]
    # Velocities
    Vm = (RPM*pi/30)*radiusDes*phi
    V1 = Vm/cos(alpha1)
    V2 = Vm/cos(alpha2)
    W1 = Vm/cos(beta1)
    W2 = Vm/cos(beta2)
    
    return PHI,PSI,Rinput,massFlow_opt,RPM,Hr_des,AR,Rgas,Gamma,inletP,inletT,c_ax,deltaH0,deltaH0is,eta,alpha1,alpha2

#% Zweifel Criteria
def Zweifel(Zcr,alpha1,alpha2):
    sigma = 2*(1/Zcr)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2
    return sigma


def optimiseMflowMeangenInput(x0,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,deltaH0is,eta,solidity,TETc,Hr_des):
    
    inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,x0*massFlow,c_ax,RPM,deltaH0is,eta,solidity,TETc)      
    # Run Meangen
    os.system('./runMeangen.sh')        
    
    f = open('meandesign.out')
    lines = f.readlines()
    f.close()
    ii = 0    
    for line in lines:
        # Extract H/R for stator or leading edge of the rotor    
        if (ii > 27) and (ii < 29):
            H_r1_meangen = float(line.split()[-1])
        # Extract solidity
        elif (ii > 51) and (ii < 53):
            Solidity_meangen = float(line.split()[-1])
        # Extract Aspect Ratio of the rotor row
        elif (ii > 52) and (ii < 54):
            AR_meangen = float(line.split()[-1])
        # Extract H/R for trailing edge of the rotor    
        elif (ii > 55) and (ii < 57):
            H_r2_meangen = float(line.split()[-1])
        elif (ii > 56) and (ii < 58):
            TETs_meangen = float(line.split()[-1])
        elif (ii > 57) and (ii < 59):
            radiusDes_meangen = float(line.split()[-1])        
        ii = ii+ 1
    
    Hr_calculated = (H_r1_meangen+H_r2_meangen)*0.5
    error_Hr = abs((1- Hr_calculated/Hr_des))
    return error_Hr

def optimiseARMeangenInput(x0,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,deltaH0is,eta,solidity,TETc,AR):
    
    inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,x0*np.array(c_ax),RPM,deltaH0is,eta,solidity,TETc)      
    # Run Meangen
    os.system('./runMeangen.sh')        
    
    f = open('meandesign.out')
    lines = f.readlines()
    f.close()
    ii = 0    
    for line in lines:
        # Extract Aspect Ratio of the rotor row
        if (ii > 24) and (ii < 26):
            AR_meangen_stator = float(line.split()[-1])
        elif (ii > 52) and (ii < 54):
            AR_meangen_rotor = float(line.split()[-1])
        ii = ii+ 1

    error_AR = abs((1- AR_meangen_stator/AR[0])) + abs((1- AR_meangen_rotor/AR[1]))
    print(error_AR,AR_meangen_stator,AR_meangen_rotor)
    return error_AR
    
#%% Write Inputs 

# Fluid
FLUID = 'air'
# Compressibility
Z = 1.0
# Stage coefficients
phi_range = np.linspace(0.5,1.2,3)
psi_range = np.linspace(0.8,2,3)
Rinput = 0.5
# Design pressure ratio
pressure_ratio = 2
# Machine design
RPM = 9000
massFlow = 25
## Geometrical parameters
# Blade height to design radius
Hr_des = 0.2
# Aspect ratio
AR = [1.2,1.2]
# TE thickness to pitch
TETs = 0.02
# Zweifel coefficient
#Zcr = 0.95
Zcr_tab = np.ones((3,3))
Zcr_tab[:,0] = np.linspace(1.1,0.7,3)
Zcr_tab[:,-1] = np.linspace(0.9,0.8,3)
for row in range(0,np.shape(Zcr_tab)[0]):
    Zcr_tab[row,:] = np.linspace(Zcr_tab[row,0],Zcr_tab[row,-1],3)
# Create path to generate folders
pathFluid = '/home/kjohri/Documents/codeThesis/'+FLUID

#%%
# Check if path exists
if os.path.isdir(pathFluid) == False:
    os.mkdir(pathFluid)

# Start iterating over Smith chart values
counter_phi = 0
for PHI in phi_range:
    counter_psi = 0
    for PSI in psi_range:
        # Preliminary calculations
        PHI,PSI,Rinput,massFlow_opt,RPM,Hr_des,AR,Rgas,Gamma,inletP,inletT,c_ax,deltaH0,deltaH0is,eta,alpha1,alpha2 = meangenInputs(PHI,PSI,Rinput,Z,pressure_ratio,massFlow,RPM,Hr_des,AR)
        # Calculate solidity to input TE thickness to chord 
        Zcr = Zcr_tab[counter_psi][counter_phi]        
        solidity = Zweifel(Zcr,alpha1,alpha2)
        TETc = (TETs/solidity)
        # Optimise meangen.in by iterating over massflow for H_rmean
        bnds = ((0.5, 5))
        res2 = minimize_scalar(optimiseMflowMeangenInput,bounds=bnds,args=(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_opt,c_ax,RPM,deltaH0is,eta,solidity,TETc,Hr_des), tol = 1e-5, method='bounded')
        massFlow_coeff = res2.x
        # Optimise meangen.in by iterating over massflow for H_rmean
        bnds_AR = ((0.01, 100),(0.01, 100))
        x0 = [1,1]
        res3 = minimize(optimiseARMeangenInput,x0,args=(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_coeff*massFlow_opt,c_ax,RPM,deltaH0is,eta,solidity,TETc,AR),method = 'Nelder-Mead',bounds=bnds_AR, tol = 1e-5)
        AR_coeff = res3.x
        # Write meangen.in file with optimum mass flow such that H_rmean calculated = H_rmean desired
        inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_coeff*massFlow_opt,AR_coeff*np.array(c_ax),RPM,deltaH0is,eta,solidity,TETc)      
        # Run Meangen
        os.system('./runMeangen.sh')        
        # Read meangen output values to finalise the input for 
        f = open('meandesign.out')
        lines = f.readlines()
        f.close()
        ii = 0    
        for line in lines:
            if (ii > 24) and (ii < 26):
                AR_meangen_stator = float(line.split()[-1])
            # Extract H/R for stator or leading edge of the rotor    
            elif (ii > 27) and (ii < 29):
                H_r1_meangen = float(line.split()[-1])
            # Extract solidity
            elif (ii > 51) and (ii < 53):
                Solidity_meangen = float(line.split()[-1])
            # Extract Aspect Ratio of the rotor row
            elif (ii > 52) and (ii < 54):
                AR_meangen_rotor = float(line.split()[-1])
            # Extract H/R for trailing edge of the rotor    
            elif (ii > 55) and (ii < 57):
                H_r2_meangen = float(line.split()[-1])
            elif (ii > 56) and (ii < 58):
                TETc_meangen = float(line.split()[-1])
            elif (ii > 57) and (ii < 59):
                radiusDes_meangen = float(line.split()[-1])        
            ii = ii+ 1
        
        Hr_calculated_final = (H_r1_meangen+H_r2_meangen)*0.5
        H2_stator_meangen = Hr_calculated_final*radiusDes_meangen
        H1_rotor_meangen = H2_stator_meangen
        H1_stator_meangen = H_r1_meangen*radiusDes_meangen        
        H2_rotor_meangen = H_r2_meangen*radiusDes_meangen        
        H_bar_stator_meangen = (H1_stator_meangen+H2_stator_meangen)*0.5
        H_bar_rotor_meangen = (H1_rotor_meangen+H2_rotor_meangen)*0.5
   
        # Calculate the optimum TETc such that TETs calculated = TETs desired                
        TETs_meangen = TETc_meangen*Solidity_meangen    
    
        # Name folders
        foldername = '{:.3f}'.format(PHI).replace('.', '')+'{:.3f}'.format(PSI).replace('.', '')
        # Make folder 
        os.mkdir(pathFluid+'/'+str(foldername))
        # Make Db folder for meshing
        os.mkdir(pathFluid+'/'+str(foldername)+'/Db')
        # Move files        
        os.system('mv *.in '+pathFluid+'/'+str(foldername))
        # Up the counter 
        counter_psi = counter_psi + 1
    counter_phi = counter_phi + 1
#%% Execution time
print("--- %s seconds ---" % (time.time() -start_time))


#%%%%

#import numpy as np
#import matplotlib.pyplot as plt
#
## Zweifel across psi
#Zcr_dist = np.ones((51,18))
#Zcr_dist[:,:] = np.linspace(1.1,0.90,51).reshape(51,1)
### Zweifel across phi and psi
#Zcr_dist1 = np.ones((51,18))
#Zcr_dist1[:,0] = np.linspace(1.1,0.7,51)
#Zcr_dist1[:,-1] = np.linspace(0.9,0.8,51)
#for row in range(0,np.shape(Zcr_dist)[0]):
#    Zcr_dist1[row,:] = np.linspace(Zcr_dist1[row,0],Zcr_dist1[row,-1],18)
## Fixed Zweifel
#Zcr = 0.95
#
#phi_range = np.linspace(0.5,1.2,18)
#psi_range = np.linspace(0.8,2,51)
#phi,psi = np.meshgrid(phi_range,psi_range)
#R = 0.5
#
#alpha1 = np.arctan(((psi/2)+1-R)/phi)
#beta1  = np.arctan(np.tan(alpha1)-(1/phi))  # radians
#alpha2 = -beta1
#
#sigma = 2*(1/Zcr)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2
#sigma_dist = 2*(1/Zcr_dist)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2
#sigma_dist1 = 2*(1/Zcr_dist1)*(np.tan(alpha1)+np.tan(alpha2))*np.cos(alpha1)**2
#
#plt.figure()
#plt.contourf(phi,psi,sigma,20)
#plt.colorbar()
#
#plt.figure()
#plt.contourf(phi,psi,sigma_dist,20)
#plt.colorbar()
#
#plt.figure()
#plt.contourf(phi,psi,sigma_dist1,20)
#plt.colorbar()























