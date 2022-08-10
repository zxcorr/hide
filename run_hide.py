# HIDE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HIDE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with HIDE.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on August 1, 2018
author: Lucas Olivari

Last update: July, 2022
authors: Joao Alberto, Carlos Otobone
'''

import numpy as np
import os

####################################################################
#### THIS IS THE SCRIPT THAT WILL RUN HIDE FOR EACH HORN
####################################################################

# ==================================================================
# CHOOSE BINGO MODEL
# ==================================================================

bingo_model = 0 # either 0, 1 or 2

# ==================================================================
# CHOOSE DESTINATION AND WORKING PATHS
# ==================================================================

import hide

destination_path = os.path.join(hide.__path__[0], '/hide/config/') # change to your destination (the place where your hide package is located within your python repository)
working_path = os.getcwd() + "/hide/config/"

# ==================================================================
# HORNS AZIMUTH AND ALTITUDE (ELEVATION)
# ==================================================================

az_in = np.loadtxt('azimuth.txt')  # one hor    n -- one azimuth [degree]       
al_in = np.loadtxt('altitude.txt') # one horn -- one altitude [degree]

# ==================================================================
# SETTING THE CONFIG FILES FOR EACH HORN
# ==================================================================

dfile_short = 'bingo_horn'

initial_horn = 0
final_horn = 1 #az_in.size

for i in range(initial_horn, final_horn):
    destination = open(working_path + 'bingo_horn_' + str(i) + '.py', 'w')
    source = open(working_path + 'bingo.py', 'r')
    for line in source:
        if line == 'coordinate_file_fmt\n':
            destination.write('coordinate_file_fmt = "coord_bingo_' + 
                             str(i) + '_' + '%s.txt"' + '\n') 
                             # coordinates file name for each horn
        
        elif line == 'params_file_fmt\n':
            destination.write('params_file_fmt = "params_bingo_' + 
                             str(i) + '_' + '{}.txt"' + '\n') 
                             # params file name for each horn
        
        elif line == 'mode\n':
            destination.write('mode = ' + str(i) + '\n')  
                             # horn suffix             

        elif line == 'azimuth_pointing\n':
            if az_in.size == 1:
                  destination.write('azimuth_pointing = ' + 
                             str(float(az_in)) + '\n')
                             # azimuth pointing (this assumes a drift scan!)
            else:
                destination.write('azimuth_pointing = ' + 
                             str(az_in[i]) + '\n')
                             # azimuth pointing (this assumes a drift scan!)

        elif line == 'altitude_start_pos\n':
            if al_in.size == 1:
                destination.write('altitude_start_pos = ' +
                             str(float(al_in)) + '\n')
                             # altitude (elevation) pointing (this assumes a drift scan!)
            else:
                destination.write('altitude_start_pos = ' +
                             str(al_in[i]) + '\n') 
                             # altitude (elevation) pointing (this assumes a drift scan!)

        elif line == 'gain_path\n':
            destination.write('gain_path = "data/gain_template_fake_bingo_model_{}_'.format(bingo_model) + 
                              str(i) + '.dat"' + '\n') 
                             # gain template used for each horn

        elif line == 'background_path\n':
            destination.write('background_path = "data/background_template_fake_bingo_model_{}_'.format(bingo_model) +
                              str(i) + '.dat"' + '\n') 
                             # background template used for each horn

        elif line == 'noise_path\n':
            destination.write('noise_path = "data/noise_template_fake_bingo_model_{}_'.format(bingo_model) +
                             str(i) + '.dat"' + '\n') 
                             # noise template used for each horn

        elif line == 'rfi_path\n':
            destination.write('rfi_path = "data/gain_template_fake_bingo_model_{}_'.format(bingo_model) +
                             str(i) + '.dat"' + '\n') 
                             # rfi template used for each horn

        else:
            destination.write(line)	
    source.close()
    destination.close()

# ==================================================================
# SETTING AND RUNNING HIDE
# ==================================================================

for i in range(initial_horn, final_horn):
    print("\nExecuting horn {0}\n".format(i))
    os.system('cp ' + working_path + 'bingo_horn_' + str(i) + '.py' + ' ' + destination_path)
    os.system('hide hide.config.' + dfile_short + '_' + str(i)) # run hide


