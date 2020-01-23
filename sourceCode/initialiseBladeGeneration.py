"""
FUNCTION: Initiate Blade Generation
@author: kjohri
"""

import os
import time
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize_scalar

# Measure time to run code

start_time = time.time()
#%% Define function to find optimum Tr for Pr = 1

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

#%% Define function to write meangen file

def inputMeangen(PHI,PSI,Rinput,FLUID,Z,pressure_ratio):
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
    PR = 1/pressure_ratio
    # Operating conditions and sizing
    massFlow  = 25.0
    RPM       = 3000
    # Work over turbine
    W_turb = eta*Cp*(Tr*Tcrit)*(1-(PR)**((float(Gamma)-1)/float(Gamma)))/1000 # kJ/kg
    
    # Write meangen
    file = open("meangen.in","w") 
    file.write("T \n")                  # TURBO_TYP,"C" FOR A COMPRESSOR,"T" FOR A TURBINE
    file.write("AXI \n")                # FLO_TYP FOR AXIAL OR MIXED FLOW MACHINE 
    file.write(str(Rgas)+" "+str(Gamma)+"\n")     # GAS PROPERTOES, RGAS, GAMMA 
    file.write(str(inletP)+" "+str(inletT)+"\n")  # POIN,  TOIN 
    file.write("1 \n")                  # NUMBER OF STAGES IN THE MACHINE 
    file.write("M \n")                  # CHOICE OF DESIGN POINT RADIUS, HUB, MID or TIP
    file.write( str(RPM)+"\n")               # ROTATION SPEED, RPM 
    file.write(str(massFlow)+"\n")           # MASS FLOW RATE, FLOWIN. 
    file.write("A \n")                  # INTYPE, TO CHOOSE THE METHOD OF DEFINING THE VELOCITY TRIANGLES
    file.write(str(R)+" "+str(phi)+" "+str(psi)+"\n")  # REACTION, FLOW COEFF., LOADING COEFF.
    file.write("B \n")                  # RADTYPE, TO CHOOSE THE DESIGN POINT RADIUS
    file.write(str(W_turb)+"  \n")      # THE DESIGN POINT RADIUS 
    file.write("0.050  0.040 \n")       # BLADE AXIAL CHORDS IN METRES.
    file.write("0.250  0.500 \n")       # ROW GAP  AND STAGE GAP 
    file.write("0.000  0.000 \n")       # BLOCKAGE FACTORS, FBLOCK_LE,  FBLOCK_TE 
    file.write(str(eta)+"  \n")             # GUESS OF THE STAGE ISENTROPIC EFFICIENCY
    file.write("1.000  1.000 \n")       # ESTIMATE OF THE FIRST AND SECOND ROW DEVIATION ANGLES
    file.write("-2.000  -2.000 \n")     # FIRST AND SECOND ROW INCIDENCE ANGLES
    file.write("1.000 \n")              # BLADE TWIST OPTION, FRAC_TWIST
    file.write("N \n")                  # BLADE ROTATION OPTION , Y or N
    file.write("92.000  88.000 \n")     # QO ANGLES AT LE  AND TE OF ROW 1 
    file.write("88.000  92.000 \n")     # QO ANGLES AT LE  AND TE OF ROW 2 
    file.write("N \n")      #DO YOU WANT TO CHANGE THE ANGLES FOR THIS STAGE ? "Y" or "N"
    file.write("Y \n")      #IFSAME_ALL, SET = "Y" TO REPEAT THE LAST STAGE INPUT TYPE AND VELOCITY TRIANGLES, SET = "C" TO CHANGE INPUT TYPE.
    file.write("Y \n")      # STATOR No.  1 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE
    file.write("Y \n")      # ROTOR No.   1 SET ANSTK = "Y" TO USE THE SAME  BLADE SECTIONS AS THE LAST STAGE 
    file.close() 
    return

#%% Initiate loop

# Input
phi_range     = np.linspace(0.5,1.3,3) 
psi_range     = np.linspace(0.8,2.5,3)
Rinput        = 0.5
Z             = 1.0
pressure_ratio= 4
FLUID         = 'air'
pathFluid     = '/home/kjohri/Documents/bladeGeneration/'+FLUID
# Create base folder based on fluid
os.mkdir(pathFluid)

for phi in phi_range:
    for psi in psi_range:
        inputMeangen(phi,psi,Rinput,FLUID,Z,pressure_ratio)
        os.system('./fortranBash.sh')
        # Name folders
        foldername = '{:.3f}'.format(phi).replace('.', '')+'{:.3f}'.format(psi).replace('.', '')
        # Make folder 
        os.mkdir(pathFluid+'/'+str(foldername))
        # Make Db folder for meshing
        os.mkdir(pathFluid+'/'+str(foldername)+'/Db')
        # Move files        
        os.system('mv *.out *.dat *.curve '+pathFluid+'/'+str(foldername))
        
#%% Execution time        
print("--- %s seconds ---" % (time.time() - start_time))
