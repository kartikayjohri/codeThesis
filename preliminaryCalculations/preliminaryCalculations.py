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
from math import sqrt,cos,tan,atan,pi
from scipy.optimize import minimize_scalar
from inputMeangen_new import inputMeangen_new

# Measure run time
start_time = time.time()
#%% Define functions 

#% Define function to find optimum Tr for Pr = 1
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
    # Velocity triangle
    alpha1 = atan(((psi/2)+1-R)/phi) # radians
    beta1  = atan(tan(alpha1)-(1/phi)) # radians
    alpha2 = -beta1 # radians
    beta2  = -alpha1 # radians

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


def optimiseMeangenInput(x0,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta,solidity,TETc,H_r_des,AR,TETs):
    
    inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,x0*massFlow,c_ax,RPM,DeltaH0,eta,solidity,TETc)      
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
    
    H_r_calculated = (H_r1_meangen+H_r2_meangen)*0.5
    error_Hr = abs((1- H_r_calculated/H_r_des))

#    c_ax_new = H_r_calculated*radiusDes_meangen/AR
#    TETc_new = TETs/solidity
#    
#    inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,x0*massFlow,c_ax_new,RPM,DeltaH0,eta,solidity,TETc_new)      
    return error_Hr
    
#%% Write Inputs 

# Fluid
FLUID = 'air'
# Compressibility
Z = 1.0
# Stage coefficients
phi_range = np.linspace(0.5,1.3,3)
psi_range = np.linspace(0.8,2.5,3)
Rinput = 0.5
# Design pressure ratio
pressure_ratio = 2
# Machine design
RPM = 9000
massFlow = 25
## Geometrical parameters
# Blade height to design radius
H_r_des = 0.2
# Aspect ratio
AR = 1.2
# TE thickness to pitch
TETs = 0.02
# Zweifel coefficient
Zcr = 0.9
# Create path to generate folders
pathFluid = '/home/kjohri/Documents/codeThesis/'+FLUID

# Check if path exists
if os.path.isdir(pathFluid) == False:
    os.mkdir(pathFluid)

# Start iterating over Smith chart values
for PHI in phi_range:
    for PSI in psi_range:
        # Preliminary calculations
        alpha1,alpha2,PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta = meangenInputs(PHI,PSI,Rinput,Z,pressure_ratio,massFlow,RPM,H_r_des,AR)
        # Calculate solidity to input TE thickness to chord 
        solidity = Zweifel(Zcr,alpha1,alpha2)
        TETc = (TETs/solidity)
        # Optimise meangen.in by iterating over massflow for H_rmean
        bnds = ((0.5, 5))
        res = minimize_scalar(optimiseMeangenInput,bounds=bnds,args=(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta,solidity,TETc,H_r_des,AR,TETs), tol = 1e-5, method='bounded')
        massFlow_coeff = res.x
        # Write meangen.in file with optimum mass flow such that H_rmean calculated = H_rmean desired
        inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_coeff*massFlow,c_ax,RPM,DeltaH0,eta,solidity,TETc)      
        # Run Meangen
        os.system('./runMeangen.sh')        
        # Read meangen output values to finalise the input for 
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
                TETc_meangen = float(line.split()[-1])
            elif (ii > 57) and (ii < 59):
                radiusDes_meangen = float(line.split()[-1])        
            ii = ii+ 1
        
        H_r_calculated = (H_r1_meangen+H_r2_meangen)*0.5
        # Calculate the optimum axial chord such that AR calculated = AR desired        
        c_ax_new = H_r_calculated*radiusDes_meangen/AR
        # Calculate the optimum TETc such that TETs calculated = TETs desired                
        TETs_meangen = TETc_meangen*Solidity_meangen
        # Write meangen.in with final values 
        inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow_coeff*massFlow,c_ax_new,RPM,DeltaH0,eta,solidity,TETc)      
        
        # Name folders
        foldername = '{:.3f}'.format(PHI).replace('.', '')+'{:.3f}'.format(PSI).replace('.', '')
        # Make folder 
        os.mkdir(pathFluid+'/'+str(foldername))
        # Make Db folder for meshing
        os.mkdir(pathFluid+'/'+str(foldername)+'/Db')
        # Move files        
        os.system('mv *.in '+pathFluid+'/'+str(foldername))

#%% Execution time
print("--- %s seconds ---" % (time.time() -start_time))