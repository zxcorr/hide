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

Update: May, 2024
authors: Alessandro Marins, Thiago Pena

Update: November, 2024
authors: Alessandro Marins, Luiza Ponte

Update: March, 2025
authors: Luiza Ponte

Update: October, 2025
authors: Nicolli Soares
'''

from __future__ import print_function, division, absolute_import, unicode_literals
from ivy.plugin.parallel_plugin_collection import ParallelPluginCollection
import os
import numpy as np
import configparser


script_dir = os.path.dirname(os.path.abspath(__file__))

ini_path = "/data/NSOARES/hide_seek/hide/hide.ini"


if not os.path.exists(ini_path):
    raise FileNotFoundError(f".ini file not found at: {ini_path}")


config = configparser.ConfigParser()
config.read(ini_path)

# ==================================================================
# LIST OF PLUGINGS THAT HIDE (IVY) WILL RUN  -- THEY ARE LOCATED IN
# THE PLUGINS DIRECTORY -- FEEL FREE TO ADD OR REMOVE PLUGINS
# ==================================================================
plugins = ["hide.plugins.initialize",
           "hide.plugins.load_beam_profile", 
           "hide.plugins.scanning_strategy",
           "hide.plugins.write_coords",
           "hide.plugins.write_calibration",
           ParallelPluginCollection([
                                    "hide.plugins.qu_opt_coord_transform",
                                     ParallelPluginCollection([
                                                            "hide.plugins.astro_signal",
                                                            "hide.plugins.earth_signal",
                                                            "hide.plugins.combine_signals", 
                                                             ],
                                                            "hide.plugins.map_frequency_plugin",
                                                            "hide.plugins.reduce_frequency_plugin",#
#                                                             parallel=False
                                                            ),
                                    #"hide.plugins.apply_gain",
                                    #"hide.plugins.add_background",           
                                    #"hide.plugins.background_noise",			
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
seed0 = config.getint('General', 'seed0')
verbose = config.getboolean('General', 'verbose')
cpu_count = config.getint('General', 'cpu_count')
backend = config.get('General', 'backend')
speed_of_light = config.getfloat('General', 'speed_of_light')
script_filename = os.path.realpath(__file__)

# ==================================================================
# OUTPUT
# ==================================================================
output_path = config.get('Output', 'output_path')
overwrite = True
file_fmt = config.get('Output', 'file_fmt')
coordinate_file_fmt
params_file_fmt
mode
calibration_file_fmt = config.get('Output', 'calibration_file_fmt')
polarizations = config.get('Output', 'polarizations').split(',')

# ==================================================================
# TELESCOPE
# ==================================================================
telescope_latitude = config.getfloat('Telescope', 'telescope_latitude')
telescope_longitude = config.getfloat('Telescope', 'telescope_longitude')
telescope_elevation = config.getfloat('Telescope', 'telescope_elevation')

# ==================================================================
# BEAM
# ==================================================================
beam_profile_provider = config.get('Beam', 'beam_profile_provider')
beam_frequency_min = config.getfloat('Beam', 'beam_frequency_min')
beam_frequency_max = config.getfloat('Beam', 'beam_frequency_max')
beam_number_channels = config.getint('Beam', 'beam_number_channels')
beam_decimals = config.getint('Beam', 'beam_decimals')
beam_nside = config.getint('Beam', 'beam_nside')
beam_response = config.getfloat('Beam', 'beam_response')
beam_elevation = config.getfloat('Beam', 'beam_elevation')
beam_azimut = config.getfloat('Beam', 'beam_azimut')
dish_diameter = config.getfloat('Beam', 'dish_diameter')
fwhm_0 = config.getfloat('Beam', 'fwhm_0')
zernike_coefficients_file_name = config.get('Beam', 'zernike_coefficients_file_name')
interpolation_scheme = config.get('Beam', 'interpolation_scheme')

# ==================================================================
# SCANNING STRATEGY
# ==================================================================
scanning_strategy_provider = config.get('ScanningStrategy', 'scanning_strategy_provider')
strategy_start = config.get('ScanningStrategy', 'strategy_start')
strategy_end = config.get('ScanningStrategy', 'strategy_end')
strategy_step_size = config.getfloat('ScanningStrategy', 'strategy_step_size')
time_range = config.getint('ScanningStrategy', 'time_range')
coord_step_size = config.getfloat('ScanningStrategy', 'coord_step_size')
azimuth_pointing
altitude_start_pos
alt_delta = config.getfloat('ScanningStrategy', 'alt_delta')
altitude_max_pos = config.getfloat('ScanningStrategy', 'altitude_max_pos')

# ==================================================================
# ASTRO (SIGNAL)
# ==================================================================
astro_signal_provider = config.get('AstroSignal', 'astro_signal_provider')
astro_signal_file_name = config.get('AstroSignal', 'astro_signal_file_name')
astro_signal_freq_file_name = config.get('AstroSignal', 'astro_signal_freq_file_name')
cache_astro_signals = config.getboolean('AstroSignal', 'cache_astro_signals')

# ==================================================================
# EARTH
# ==================================================================
earth_signal_provider = config.get('EarthSignal', 'earth_signal_provider')
earth_signal_flux = config.getfloat('EarthSignal', 'earth_signal_flux')

# ==================================================================
# BACKGROUND
# ==================================================================
elevation_model = [float(x) for x in config.get('Background', 'elevation_model').split(',')]

# ==================================================================
# NOISE
# ==================================================================
load_noise_template = config.getboolean('Noise', 'load_noise_template')
temp_sys = config.getfloat('Noise', 'temp_sys')
delta_nu = config.getfloat('Noise', 'delta_nu')
color_alpha = config.getfloat('Noise', 'color_alpha')
color_fknee = config.getfloat('Noise', 'color_fknee')
color_beta = config.getfloat('Noise', 'color_beta')
sample_freq = config.getfloat('Noise', 'sample_freq')

# ==================================================================
# POST PROCESSING
# ==================================================================
instrument = config.get('PostProcessing', 'instrument')
seed
gain_path
background_path
noise_path
rfi_path

# ==================================================================
# RFI
# ==================================================================
load_rfi_template = config.getboolean('RFI', 'load_rfi_template')
rfideltat = config.getfloat('RFI', 'rfideltat')
rfideltaf = config.getfloat('RFI', 'rfideltaf')
rfiexponent = config.getfloat('RFI', 'rfiexponent')
rfienhance = config.getfloat('RFI', 'rfienhance')
rfiday = [float(x) for x in config.get('RFI', 'rfiday').split(',')]
#rfiday = tuple(map(float, config.get('RFI', 'rfiday').split(',')))
rfidamping = config.getfloat('RFI', 'rfidamping')
