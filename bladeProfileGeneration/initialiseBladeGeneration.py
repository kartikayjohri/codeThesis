"""
FUNCTION: Initiate Blade Profile Generation
@author: kjohri
"""

import os
import time

# Measure time to run code
start_time = time.time()

# Specify fluid name 
FLUID         = 'air'
pathFluid     = '/home/kjohri/Documents/codeThesis/'+FLUID
folders = sorted(os.listdir(pathFluid))

for folder in folders:
    # Move file from case folder to bladeProfileGeneration folder     
    os.system('mv '+pathFluid+'/'+folder+'/meangen.in /home/kjohri/Documents/codeThesis/bladeProfileGeneration')
    # Change directory    
    os.chdir('/home/kjohri/Documents/codeThesis/bladeProfileGeneration')    
    # Execute meangen and stagen and bladeProfile.py    
    os.system('./fortranBash.sh')    
    # Move files from current directory to case folder       
    os.system('mv *.in *.out *.dat *.curve '+pathFluid+'/'+folder)
    
#%% Execution time        
print("--- %s seconds ---" % (time.time() - start_time))
