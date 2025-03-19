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

Update: May, 2024
authors: Alessandro Marins, Thiago Pena

Update: March, 2025
authors: Luiza Ponte
'''

import numpy as np
import os
import re
import configparser 


script_dir = os.path.dirname(os.path.abspath(__file__))


ini_path = os.path.join(script_dir, "hide.ini")


if not os.path.exists(ini_path):
    raise FileNotFoundError(f".ini file not found in {ini_path}")


config = configparser.ConfigParser()
config.read(ini_path)

# ==================================================================
# CHOOSE BINGO MODEL
# ==================================================================
bingo_model = config.getint('RunHide', 'bingo_model')

# ==================================================================
# CHOOSE DESTINATION AND WORKING PATHS
# ==================================================================

destination_path = config.get('RunHide', 'destination_path')
working_path = config.get('RunHide', 'working_path')

# ==================================================================
# HORNS AZIMUTH AND ALTITUDE (ELEVATION)
# ==================================================================

az_in_file = config.get('RunHide', 'az_in_file')
al_in_file = config.get('RunHide', 'al_in_file')
az_in = np.loadtxt(az_in_file)  # one horn -- one azimuth [degree]
al_in = np.loadtxt(al_in_file)  # one horn -- one altitude [degree]


# ==================================================================
# SETTING THE CONFIG FILES FOR EACH HORN
# ==================================================================
if len(os.sys.argv) > 1:
    base_file = os.sys.argv[1]  # base file, as "bingo.py"
    dfile_short = os.sys.argv[2]  # new files fmt, as "bingo_horn_{}"
    initial_horn = int(os.sys.argv[3])
    final_horn = int(os.sys.argv[4])
else:
    base_file = "bingo.py"
    dfile_short = 'bingo_horn'
    initial_horn = 1
    final_horn = 140

dfile_short_py = dfile_short + '.py'

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

        elif line == 'seed\n':
            destination.write('seed = seed0 + {}\n'.format(i))

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
    os.system('hide hide.config.' + dfile_short + '_' + str(i))  # run hide
