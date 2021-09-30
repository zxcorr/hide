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
#### AND NOISE) BY FOLLOWING MODEL 2 (BASED ON BLEIEN, SEE DOCUMENTATION)
####
#### FOR THE FUTURE: ADD RFI  
#######################################################################

# ==================================================================
# CHOOSE DESTINATION PATHS
# ==================================================================

destination_path = '/usr/local/lib/python2.7/dist-packages/hide-0.1.0-py2.7.egg/hide/data/' # change to your destination (the place where your hide package is located within your python repository)

# ==================================================================
# READ IN BLEIEN INFO
# ==================================================================

gain_in = np.loadtxt("gain_template_7m_FFT_phase_switch_ADU_K.dat")
background_in = np.loadtxt("background_template_7m_FFT_ADU.dat")
noise_in = np.loadtxt("noise_template_7m_FFT_ADU.dat")
rfi_in = np.loadtxt("rfi_template_7m_FFT_ADU.dat") # for completeness

amplitude_gain = np.mean(gain_in[:, 1]) # mean Bleien gain
amplitude_background = np.mean(background_in[:, 1]) # mean Bleien background

# ==================================================================
# BINGO PARAMETERS
# ==================================================================

number_horns = 1 # number of horns to be simulated
freq_min = 960. # minimun frequency, in MHz
freq_bin = 1. # fequency bin, in MHz
n_channels = 300 # number of channels

frequencies = freq_min + freq_bin * np.arange(n_channels)

# ==================================================================
# BINGO CALIBRATION PARAMETERS
# ==================================================================

calibration_error = 5. # in percentage
delta_nu_osc = 50. # in MHz
amplitude_osc = 0.15 # dimensionless

# ==================================================================
# CALCULATION
# ==================================================================

# Loop on the horns
for i in range(0, number_horns):

# GAIN MODEL
    sigma_gain = (calibration_error / 100.)
    amplitude_gain_random = np.random.normal(loc=0., scale=sigma_gain, size=n_channels) # draw the amplitude from a Gaussian distribution with mean = 0 and standard deviation = sigma_gain
    bingo_fake_gain = amplitude_gain * (1. + amplitude_gain_random + amplitude_osc * np.sin((np.pi * frequencies)/delta_nu_osc)) # BINGO Fake Gain model 2 (Kelvin to ADU)
    
# BASELINE (BACKGROUND) MODEL
    #amplitude_osc_background = (amplitude_osc_percentage / 100.)
    bingo_fake_background = amplitude_background * (1. + amplitude_osc * np.sin((np.pi * frequencies)/delta_nu_osc)) # BINGO Fake Amplitude model 2, in ADU
     
# OUTPUT
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

    dfile = "gain_template_fake_bingo_model_2_" + str(i) + ".dat"
    np.savetxt(dfile, output_gain)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "background_template_fake_bingo_model_2_" + str(i) + ".dat"
    np.savetxt(dfile, output_background)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "noise_template_fake_bingo_model_2_" + str(i) + ".dat"
    np.savetxt(dfile, output_noise)
    os.system('cp ' + dfile + ' ' + destination_path)

    dfile = "rfi_template_fake_bingo_model_2_" + str(i) + ".dat"
    np.savetxt(dfile, output_rfi)
    os.system('cp ' + dfile + ' ' + destination_path)


