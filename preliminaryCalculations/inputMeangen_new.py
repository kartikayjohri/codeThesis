"""
FUNCTION: WRITE MEANGEN INPUT FILE
@author: kjohri
"""

def inputMeangen_new(PHI,PSI,Rinput,Rgas,Gamma,inletP,inletT,massFlow,c_ax,RPM,DeltaH0,eta,TETc):
    # Write meangen
    file = open("meangen.in","w") 
    file.write("T \n")                  # TURBO_TYP,"C" FOR A COMPRESSOR,"T" FOR A TURBINE
    file.write("AXI \n")                # FLO_TYP FOR AXIAL OR MIXED FLOW MACHINE 
    file.write(str(Rgas)+" "+str(Gamma)+"\n")      # GAS PROPERTOES, RGAS, GAMMA 
    file.write(str(inletP)+" "+ str(inletT)+"\n")  # POIN,  TOIN 
    file.write("1 \n")                  # NUMBER OF STAGES IN THE MACHINE 
    file.write("M \n")                  # CHOICE OF DESIGN POINT RADIUS, HUB, MID or TIP
    file.write(str(RPM)+"\n")           # ROTATION SPEED, RPM 
    file.write(str(massFlow)+"\n")      # MASS FLOW RATE, FLOWIN. 
    file.write("A \n")                  # INTYPE, TO CHOOSE THE METHOD OF DEFINING THE VELOCITY TRIANGLES
    file.write(str(Rinput)+" "+str(PHI)+" "+str(PSI)+"\n")  # REACTION, FLOW COEFF., LOADING COEFF.
    file.write("B \n")                  # RADTYPE, TO CHOOSE THE ENTHALPY DROP
    file.write(str(DeltaH0/1000)+"\n")  # ENTHALPY DROP IN KJ/kg
    file.write("0.050"+str(c_ax)+"\n")  # BLADE AXIAL CHORDS IN METRES.
    file.write("0.250  0.500 \n")       # ROW GAP  AND STAGE GAP 
    file.write("1.000  1.000 \n")       # EXPONENT TO MOVE BLADE LOADING TK_TYP_S AND TK_TYP_R
    file.write("0.040  0.040 \n")       # LE THICKNESS TO CHORD RATIO  
    file.write("0.040"+str(TETc)+"\n")  # TE THICKNESS TO CHORD RATIO  
    file.write("2.000  2.000 \n")       # THICKNESS DISTRIBUTION  
    file.write("0.000  0.000 \n")       # BLOCKAGE FACTORS, FBLOCK_LE,  FBLOCK_TE 
    file.write(str(eta)+"\n")           # GUESS OF THE STAGE ISENTROPIC EFFICIENCY
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