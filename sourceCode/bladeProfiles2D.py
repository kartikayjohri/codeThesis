# Import modules
import numpy as np

# Make function to generate coordinate files
def Coords(filename):
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
    template = '{0:.6f} {1:.6f}'
    f.write('# Profile ' + str(2) + ' at 50%')
    f.write('\n')
    for ii in range(len(x_st)):
        f.write(template.format(x_st[ii], y_st[ii]))
        f.write('\n')
    f.write('\n')
    f.close()
    
    f = open('rotor_profile.curve', 'w')
    f.write('# Profile ' + str(2) + ' at 50%')
    f.write('\n')
    for ii in range(len(x_rot)):
        f.write(template.format(x_rot[ii], y_rot[ii]))
        f.write('\n')
    f.write('\n')
    f.close()
      
    return x_rot, y_rot, x_st, y_st

# Execute function
Coords('blade_profiles.dat')