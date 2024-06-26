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
Created on August, 2018

author: Lucas Olivari
'''

import numpy as np
import os

#######################################################################
#### THIS SCRIPT SIMULATES A FAKE BINGO SPECTROMETER (GAIN, BACKGROUND,
#### AND NOISE) BY FOLLOWING MODEL 1 (BASED ON BLEIEN, SEE DOCUMENTATION)
####
#### FOR THE FUTURE: ADD RFI  
#######################################################################

# ==================================================================
# CHOOSE DESTINATION PATHS
# ==================================================================

destination_path = '/scratch/bingo/thiago.pena/extra/venv/lib/python3.8/site-packages/hide-0.1.0-py3.8.egg/hide/data/' # change to your destination (the place where your hide package is located within your python repository)

# ==================================================================
# READ IN BLEIEN INFO
# ==================================================================

rfi_in = np.loadtxt("rfi_template_7m_FFT_ADU.dat") # for completeness               

# ==================================================================
# BINGO PARAMETERS
# ==================================================================

number_horns = 140 # number of horns to be simulated
freq_min = 960. # minimun frequency, in MHz
freq_bin = 10. # fequency bin, in MHz
n_channels = 30 # number of channels

frequencies = freq_min + freq_bin * np.arange(n_channels)

# ==================================================================
# GAIN = 1 // RFI = 0 // BACKGROUND = 0 // NOISE = 0
# ==================================================================

# Loop on the horns
for i in range(0, number_horns):

    bingo_fake_gain = np.full(n_channels, 1)
    bingo_fake_background = np.zeros(n_channels)

    output_gain = np.zeros((n_channels, 2))
    output_gain[:, 0] = frequencies
    output_gain[:, 1] = bingo_fake_gain

    output_background = np.zeros((n_channels, 2))
    output_background[:, 0] = frequencies
    output_background[:, 1] = bingo_fake_background

    output_noise = np.zeros((n_channels, 2))
    output_noise[:, 0] = frequencies
    output_noise[:, 1] = bingo_fake_gain

    output_rfi = np.zeros((n_channels, 3))
    output_rfi[:, 0] = frequencies 
    output_rfi[:, 1] = np.mean(rfi_in[:, 1])
    output_rfi[:, 2] = np.mean(rfi_in[:, 2])

    dfile = "gain_template_fake_bingo_model_0_" + str(i) + ".dat"
    np.savetxt(dfile, output_gain)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "background_template_fake_bingo_model_0_" + str(i) + ".dat"
    np.savetxt(dfile, output_background)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "noise_template_fake_bingo_model_0_" + str(i) + ".dat"
    np.savetxt(dfile, output_noise)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "rfi_template_fake_bingo_model_0_" + str(i) + ".dat"
    np.savetxt(dfile, output_rfi)
    os.system('cp ' + dfile + ' ' + destination_path)
