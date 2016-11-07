import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import time

starttime = time.clock()

#The program goes through all the dat files in the folder and writes the relevant data into a new file which can then be used for analysis

filenum = 0

bin = input("\nBin: ")

name = str(bin) + 'Combined'

f = open(name, 'w')

for filename in glob.glob('*.dat'):

    filenum += 1
	
    print(filenum)

	#The file in the following if statement was corrupt, remove if looking at new files
	
    if filename == 'Beam3_dm_D20140307T171132.dat':
        continue
    
    x = np.loadtxt(filename,
    dtype = {'names': ('Time', 'DM', 'SNR', 'Bin'), 'formats': ('f8', 'f8', 'f8', 'int8')}, delimiter = ', ')
    
    r = 0

    for row in x:
        if x[r]['Bin'] == bin and x[r]['DM'] == 3.0:
            f.write(str(row))
            f.write('\n')
        r += 1
        
f.close()
    
endtime = time.clock()

totaltime = endtime - starttime

print('The program took', totaltime, 'seconds to run')














