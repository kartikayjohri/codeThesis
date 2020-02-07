# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 19:40:29 2020

@author: kjohri
"""

"""
FUNCTION: Initialise blade meshing
@author: kjohri
"""
import os
import numpy as np

import matplotlib.pyplot as plt

fluid = 'air'

# Extract pitch and chord from meangen.out per folder
folders = sorted(os.listdir("/home/kjohri/Documents/bladeGeneration/"+fluid))
number = 2
for number in range(2,3):
    folder = folders[number]
    os.chdir("/home/kjohri/Documents/bladeGeneration/air/"+folder)
    f = open('meandesign.out')
    lines = f.readlines()
    f.close()
    ii = 0
    stator_value = []    
    rotor_value = []
    for line in lines:
        # Extract flow angles and other parameters of the stator
        if (ii > 39) and (ii < 41):
            chord = (float(line.split()[-1]))
        if (ii > 41) and (ii < 43):
            pitch = (float(line.split()[-1]))
        ii = ii+ 1
    
    # Extract all coordinates from profile
    f = open('rotor_profile.curve')
    lines = f.readlines()
    f.close()
    
    # Set empty array
    x_coords = []
    y_coords = []
    
    ii = 0
    for line in lines:
        if (ii > 0) and (ii <+ len(lines)-1):
            x_coords.append(float(line.split()[0])*chord)
            y_coords.append(float(line.split()[1])*chord)
        ii = ii+ 1 
    
    # Segregate blade coordinates in upper and lower
    x_blade_u = x_coords[int(len(x_coords)/2):]
    x_blade_u = x_blade_u[::-1]
    x_blade_l = x_coords[0:int(len(x_coords)/2)]
    
    y_blade_u = np.array(y_coords[int(len(y_coords)/2):])
    y_blade_u = y_blade_u[::-1]
    y_blade_l = np.array(y_coords[0:int(len(y_coords)/2)])
    
    # Blade thickness,camber
    t_blade = abs(y_blade_u-y_blade_l)
    camber_blade = y_blade_u - (t_blade/2) # Not the best method because coords are determined based on tangent
    t_blade_l = abs(camber_blade - y_blade_l)
    t_blade_u = abs(camber_blade - y_blade_u)
    
    # Inflow and outflow bc
    bc_in = 0.25*chord
    bc_out = 1.25*chord 
    
    # Reference point based on thickness
    x_ref = x_blade_u[np.argmax(t_blade)]
    y_ref = camber_blade[np.argmax(t_blade)]
    
    # Inlet coordinates
    x_in = [-bc_in,-bc_in]
    y_in = [y_ref-pitch/2,y_ref+pitch/2]
    
    # Upper periodic gap - inlet
    x_in_u = [-bc_in,x_blade_u[np.argmax(t_blade)]]
    y_in_u = [y_in[-1],y_in[-1]]
    
    # Lower periodic gap - inlet
    x_in_l = [-bc_in,x_blade_u[np.argmax(t_blade)]]
    y_in_l = [y_in[0],y_in[0]]
    
    # Determining boundary coordinates
    bc_u = [y_in_u[-1],y_ref+camber_blade[-1]+pitch/2]
    bc_l = [y_in_l[-1],y_ref+camber_blade[-1]-pitch/2]
    bc_x = [x_blade_u[np.argmax(t_blade)],chord]
    
    # Upper periodic gap - outlet
    x_out_u = [bc_x[-1],bc_out]
    y_out_u = [bc_u[-1],bc_u[-1]]
    
    # Lower periodic gap - outlet
    x_out_l = [bc_x[-1],bc_out]
    y_out_l = [bc_l[-1],bc_l[-1]]
    
    # Outlet coordinates
    x_out = [bc_out,bc_out]
    y_out = [bc_u[-1],bc_l[-1]]
    
    
    #% Plotting
    plt.figure(number)
    plt.plot(x_coords,y_coords)
    plt.plot(x_coords,[pitch+i for i in y_coords])
    plt.plot(x_in,y_in)
    plt.plot(x_in_u,y_in_u)
    plt.plot(x_in_l,y_in_l)
    plt.plot(bc_x,bc_u)
    plt.plot(bc_x,bc_l)
    plt.plot(x_out_u,y_out_u)
    plt.plot(x_out_l,y_out_l)
    plt.plot(x_out,y_out)
    plt.xlabel('Meridional flow direction')
    plt.ylabel('Tangential flow direction')
    plt.grid()
    plt.gca().set_aspect("equal")
    plt.plot(x_blade_u,camber_blade,'--')
    plt.title('$\phi$='+str(int(folder[0:4])/1000)+','+'$\psi$='+str(int(folder[4:])/1000))


#%% Generating Geometry file
os.chdir("/home/kjohri/Documents/bladeGeneration/air/"+folder+'/Db')

# Create and write in new file
f = open('geometry.rotor', 'w')
template = '{0:.8f} {1:.8f}'

f.write('Number of surfaces \n')
f.write('10 \n')

surface_coords = [y_blade_u,y_blade_l]
for i in range(0,len(surface_coords)):
    f.write('BLADE \n')
    f.write("       ' S '\n")
    f.write('       dim         np\n')
    f.write('         2        '+str(len(x_blade_u))+'\n')
    f.write('         x          y\n')
    for ii in range(0,len(x_blade_u)):
        f.write(template.format(x_blade_u[ii], surface_coords[i][ii]))
        f.write('\n')

f.write('INFLOW \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_in))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_in)):
    f.write(template.format(x_in[i], y_in[i]))
    f.write('\n')

f.write('PERIODIC IN L \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_in_l))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_in_l)):
    f.write(template.format(x_in_l[i], y_in_l[i]))
    f.write('\n')    

f.write('PERIODIC BLD L \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(bc_x))+'\n')
f.write('         x          y\n')
for i in range(0,len(bc_x)):
    f.write(template.format(bc_x[i], bc_l[i]))
    f.write('\n')

f.write('PERIODIC OUT L - Lower \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_out_l))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_out_l)):
    f.write(template.format(x_out_l[i], y_out_l[i]))
    f.write('\n')

f.write('PERIODIC IN U \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_in_u))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_in_u)):
    f.write(template.format(x_in_u[i], y_in_u[i]))
    f.write('\n')

f.write('PERIODIC BLD U \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(bc_x))+'\n')
f.write('         x          y\n')
for i in range(0,len(bc_x)):
    f.write(template.format(bc_x[i], bc_u[i]))
    f.write('\n')

f.write('PERIODIC OUT U - Lower \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_out_u))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_out_u)):
    f.write(template.format(x_out_u[i], y_out_u[i]))   
    f.write('\n')

f.write('OUTFLOW \n')
f.write("       ' S '\n")
f.write('       dim         np\n')
f.write('         2        '+str(len(x_out))+'\n')
f.write('         x          y\n')
for i in range(0,len(x_out)):
    if i == 0:
        f.write(template.format(x_out[i], y_out[i]))
        f.write('\n')
    else:
        f.write(template.format(x_out[i], y_out[i]))

f.close()
