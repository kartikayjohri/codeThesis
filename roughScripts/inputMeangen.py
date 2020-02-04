"""
FUNCTION: WRITE MEANGEN INPUT FILE
@author: kjohri
"""

#def inputMeangen(PHI,PSI,Rinput):
# Gas properties
Rgas  = str(287.15)
Gamma = str(1.40)
# Operating conditions and sizing
inletP    = str(8.0)
inletT    = str(600.0)
massFlow  = str(25.0)
RPM       = str(3000)
# Velocity triangles
phi = str(0.5)#PHI)
psi = str(4.0)#PSI)
R   = str(0.5)#Rinput)
# Write meangen
file = open("meangen.in","w") 
file.write("T \n")                  # TURBO_TYP,"C" FOR A COMPRESSOR,"T" FOR A TURBINE
file.write("AXI \n")                # FLO_TYP FOR AXIAL OR MIXED FLOW MACHINE 
file.write(Rgas+" "+Gamma+"\n")     # GAS PROPERTOES, RGAS, GAMMA 
file.write(inletP+" "+inletT+"\n")  # POIN,  TOIN 
file.write("1 \n")                  # NUMBER OF STAGES IN THE MACHINE 
file.write("M \n")                  # CHOICE OF DESIGN POINT RADIUS, HUB, MID or TIP
file.write( RPM+"\n")               # ROTATION SPEED, RPM 
file.write(massFlow+"\n")           # MASS FLOW RATE, FLOWIN. 
file.write("A \n")                  # INTYPE, TO CHOOSE THE METHOD OF DEFINING THE VELOCITY TRIANGLES
file.write(R+" "+phi+" "+psi+"\n")  # REACTION, FLOW COEFF., LOADING COEFF.
file.write("A \n")                  # RADTYPE, TO CHOOSE THE DESIGN POINT RADIUS
file.write("0.400  \n")             # THE DESIGN POINT RADIUS 
file.write("0.050  0.040 \n")       # BLADE AXIAL CHORDS IN METRES.
file.write("0.250  0.500 \n")       # ROW GAP  AND STAGE GAP 
file.write("0.000  0.000 \n")       # BLOCKAGE FACTORS, FBLOCK_LE,  FBLOCK_TE 
file.write("0.850  \n")             # GUESS OF THE STAGE ISENTROPIC EFFICIENCY
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
#return

import subprocess
subprocess.call("testBash.sh")
