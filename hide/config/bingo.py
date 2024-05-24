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

Modified on 2022
authors: Joao Alberto, Carlos Otobone
'''

from __future__ import print_function, division, absolute_import, unicode_literals

from ivy.plugin.parallel_plugin_collection import ParallelPluginCollection

import os
import numpy as np
####################################################################
#### THIS IS THE CONFIGURATION FILE THAT WILL SIMULATE BINGO
####################################################################

# ==================================================================
# LIST OF PLUGINGS THAT HIDE (IVY) WILL RUN  -- THEY ARE LOCATED IN
# THE PLUGINS DIRECTORY -- FEEL FREE TO ADD OR REMOVE PLUGINS
# ==================================================================
plugins = ["hide.plugins.initialize",
           "hide.plugins.load_beam_profile", #
           "hide.plugins.scanning_strategy",
           "hide.plugins.write_coords",
           "hide.plugins.write_calibration",
           ParallelPluginCollection([
                                    "hide.plugins.qu_opt_coord_transform",
                                     ParallelPluginCollection([
                                                            "hide.plugins.astro_signal",#
                                                            "hide.plugins.earth_signal",
                                                            "hide.plugins.combine_signals", #MODIFIQUEI, 11/2023
                                                             ],
                                                            "hide.plugins.map_frequency_plugin",#
                                                            "hide.plugins.reduce_frequency_plugin",# 
#                                                             parallel=False
                                                            ),
                                    #"hide.plugins.apply_gain",
                                    "hide.plugins.add_background",           
                                    "hide.plugins.background_noise",			# FOR WHITE NOISE
                                    "hide.plugins.write_tod_phaseswitch",
                                    "hide.plugins.clean_up",
                                     ],
                                    "hide.plugins.map_strategy_plugin",
                                    ),
            "hide.plugins.write_params",
            "ivy.plugin.show_summary_stats"
          ]

# ==================================================================
# GENERAL
# ==================================================================

seed0 = 2700

verbose = True
cpu_count = 1
backend = "sequential"
speed_of_light = 299792458          # [m/s]

script_filename = os.path.realpath(__file__)

# ==================================================================
# OUTPUT
# ==================================================================
output_path = "/scratch/bingo/thiago.pena/extra/beam_test/240412/seedtest{}/".format(seed0)  # path to output folder
overwrite = True
file_fmt = "bingo_tod_horn_{mode}_{date}.h5"     # tod file name 
coordinate_file_fmt
                   # it will be written by run_hide.py
params_file_fmt
                   # it will be written by run_hide.py
mode
                   # it will be written by run_hide.py

### TO DISCOVER
calibration_file_fmt = 'CALIBRATION_{date}.txt'
polarizations = ['PXX']
###

# ==================================================================
# TELESCOPE
# ==================================================================

#-----------------
# BINGO
#-----------------
telescope_latitude = -7.0
telescope_longitude = -38.0

telescope_elevation = 0.0        # altitude

# ==================================================================
# BEAM
# ==================================================================
#beam_profile_provider = "hide.beam.beamz"
beam_profile_provider = "hide.beam.gaussian_fwhm"  
beam_frequency_min = 960.       # minimum frequency: [MHz]
beam_frequency_max = 1260.      # maximum frequency: [MHz]            #last point discarted
beam_frequency_pixscale = 10.    # pixel scale (frequency/beam)
beam_nside = 256                 # healpix NSIDE -- must be the same as that of the input sky maps
beam_response = 1               # beam response [0..1]
beam_elevation = 5          # elevation [degree] -- to calculate the beam area (not physical, just to find the pixels for calculation)
beam_azimut = 5               # azimuth [degree] -- same as above

# Gaussian (gaussian)
dish_diameter = 40.             # effective diameter of the dish [m]
# Gaussian FWHM (gaussian_fwhm)
fwhm_0 = 0.011                  # FWHM for the minimum frequency
# Zernike Beam (beamz)
zernike_coefficients_file_name = "DoubleRectangular_hornX_displacementY.fits"
interpolation_scheme = "nearest" # [nearest, linear, cubic]


# ==================================================================
# SCANNING STRATEGY
# ==================================================================
scanning_strategy_provider = "hide.strategy.drift_scan" 
strategy_start = "2020-03-01-00:00:00"     # survey start time. Format YYYY-mm-dd-HH:MM:SS
strategy_end   = "2020-03-01-23:59:59"     # survey end time. Format YYYY-mm-dd-HH:MM:SS
strategy_step_size = 1                    # size of step in [sec]
time_range = 60*60                        # time range per file [sec]
coord_step_size = 1                        # step size in the coords file

# -------------------
# DRIFT SCAN
# -------------------
azimuth_pointing
                # it will be written by run_hide.py
altitude_start_pos
                # it will be written by run_hide.py

### NOT IMPORTANT FOR BINGO!
alt_delta = 0.                             # change in altitude per day (drift scan) [degree]
altitude_max_pos = 90.0                    # max position in altitude direction [degree]
###

# ==================================================================
# ASTRO (SIGNAL)
# ==================================================================
astro_signal_provider = "hide.astro.hi_sky"    # it will read the SKY maps
astro_signal_file_name = "sky_960mhz1260mhz_nch30_mk_fullsky_nonoise_nobeam.fits" # maps (n_channels vs n_pixels) file name, located in the data/sky directory
astro_signal_freq_file_name = "oldfreqs_bingo.fits" # frequency (n_channels) file name, located in the data/sky directory

cache_astro_signals = True         # flag if loaded signals per frequency should be kept in memory

# ==================================================================
# EARTH
# ==================================================================
earth_signal_provider = "hide.earth.constant"
# -------------------
# c o n s t a n t 
# -------------------
earth_signal_flux = 0              # flux of constant earth signal


# ==================================================================
# BACKGROUND -- AKERET ET AL. (2017) MODEL
# ==================================================================
elevation_model = [0., 0., 0.] # chose [1, -1., 1.] for model 1 and [1.26321397e+10, -1.71282810e+10, 2.79280833e+10] for model 2

# ==================================================================
# NOISE -- HARPER ET AL. (2018) MODEL
# ==================================================================
load_noise_template = False
temp_sys = 70                   # system temperature, in K
year_adj = 0.0523424      # 1/sqrt(365), adjusting noise amplitude from 1d to 1y
temp_sys  *= year_adj
delta_nu = 10*1e6                   # channel width, in Hz
color_alpha = 0                  # 1 /f alpha parameter
color_fknee = 1                 # 1/ f knee frequency, in Hz
color_beta = 1                 # 1 / f beta parameter (0.001 - 1)
sample_freq = 1.0                 # telescope sample rate, in Hz

# ==================================================================
# POST PROCESSING
# ==================================================================
instrument = "hide.spectrometer.fake_bingo_spectrometer"

seed

gain_path
         # it will be written by run_hide.py
background_path
         # it will be written by run_hide.py
noise_path
         # it will be written by run_hide.py
rfi_path
         # it will be written by run_hide.py

# ==================================================================
# RFI -- AKERET ET AL. (2017) MODEL -- NOT BEING USED FOR BINGO AT 
# THE MOMENT
# ==================================================================
load_rfi_template = False
rfideltat = 5                           # Width in time for RFI [units of pixels]
rfideltaf = .5                          # Width in frequency for RFI [units of pixels]
rfiexponent = 2                         # Exponential model (1) or Gaussian model (2) for RFI
rfienhance = 1.7                        # Enhance fraction covered by RFI
rfiday = (6.0, 22.0)                    # Beginning and end of RFI day
rfidamping = 0.1                        # Damping factor of RFI during the RFI night 
