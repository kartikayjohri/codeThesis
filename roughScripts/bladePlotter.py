import numpy as np
import matplotlib.pyplot as plt



def PLOT(filename,Fluid):
    path = "/home/kjohri/Documents/bladeGeneration/"+Fluid
    f = open(filename,'r')
    lines = f.readlines()
    x_2D_st = []
    ySS_2D_st = []
    yPS_2D_st = []
    x_2D_rot = []
    ySS_2D_rot = []
    yPS_2D_rot = []
    
    ii = 0
    for line in lines:
        if (ii >= 320) and (ii <+ 520):
            x_2D_st.append(float(line.split()[0]))
            ySS_2D_st.append(float(line.split()[1]))
            yPS_2D_st.append(float(line.split()[2]))
        elif (ii >= 1271) and (ii <+ 1471):
            x_2D_rot.append(float(line.split()[0]))
            ySS_2D_rot.append(float(line.split()[1]))
            yPS_2D_rot.append(float(line.split()[2]))
        ii += 1
    f.close()
    
    x_st = np.concatenate((np.asarray(x_2D_st), np.asarray(x_2D_st[::-1])))
    y_st = np.concatenate((np.asarray(yPS_2D_st), np.asarray(ySS_2D_st[::-1])))
    x_rot = np.concatenate((np.asarray(x_2D_rot), np.asarray(x_2D_rot[::-1])))
    y_rot = np.concatenate((np.asarray(yPS_2D_rot), np.asarray(ySS_2D_rot[::-1])))
    
    f = open('stator_profile.curve', 'w')
    f.write('# Profile ' + str(3) + ' at 50%')
    f.write('\n')
    for ii in range(len(x_st)):
        f.write('    ' + str(x_st[ii]) + '    ' + str(y_st[ii]) + '    ' + str(0.0))
        f.write('\n')
    f.write('\n')
    f.write('\n')
    f.close()
    
    f = open('rotor_profile.curve', 'w')
    f.write('# Profile ' + str(3) + ' at 50%')
    f.write('\n')
    for ii in range(len(x_rot)):
        f.write('    ' + str(x_rot[ii]) + '    ' + str(y_rot[ii]) + '    ' + str(0.0))
        f.write('\n')
    f.write('\n')
    f.write('\n')
    f.close()
    
      z
#    fig_st = plt.figure()
#    ax_st = fig_st.add_subplot(111)
#    fig_rot = plt.figure()
#    ax_rot = fig_rot.add_subplot(111)
#    ax_st.set_title('stator')
#    ax_st.plot(x_st, y_st, 'r')
#    ax_rot.set_title('rotor')
#    ax_rot.plot(x_rot, y_rot, 'b')
#    plt.show()
#    plt.hold(True)
    return x_rot, y_rot, x_st, y_st
# %%
from itertools import cycle
    
filenames = ['blade_profiles.dat']
#filenames = ['blade_profiles0515.dat','blade_profiles05275.dat','blade_profiles0540.dat']
#filenames = ['blade_profiles1315.dat','blade_profiles13275.dat','blade_profiles1340.dat']
#filenames = ['blade_profiles0915.dat','blade_profiles09275.dat','blade_profiles0940.dat']
lines = ["-","--","-.",":"]
linecycler = cycle(lines)
plt.figure()
for filename in filenames:
    x_rot = PLOT(filename)[0]
    y_rot = PLOT(filename)[1]    
    
    plt.title('rotor')
    plt.plot(x_rot, y_rot,next(linecycler))
    plt.show()
    plt.hold(True)    
plt.figure()
for filename in filenames:
    x_st = PLOT(filename)[2]
    y_st = PLOT(filename)[3]
    plt.title('stator')
    plt.plot(x_st,y_st,next(linecycler))
    
#%%
    
f = open("/home/kjohri/Documents/bladeGeneration/He/05000800/rotor_profile.curve",'r')
lines = f.readlines()

x = []
y = []

N = 0
for line in lines:
    if N > 0:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))
    N = N+1

import matplotlib.pyplot as plt
plt.plot(x,y)
