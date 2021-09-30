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

from __future__ import print_function, division, absolute_import, unicode_literals

import glob
import os
from pkg_resources import resource_filename

import healpy as hp
import numpy as np
from astropy.io import fits as pyfits

import hide

HISKY_FILE_PATH = "/hide/data/sky/"

def load_signal(ctx):
    """
    Returns an IMSKY map dependent on the frequency.
    
    :param params: The ctx instance with the paramterization
    :returns signal: The astro signal
    """

    path = os.getcwd() # it will search for the sky data in your working directory!

    hisky_maps_file = ctx.params.astro_signal_file_name
    hisky_freq_file = ctx.params.astro_signal_freq_file_name

    gsm_maps = pyfits.getdata(path + HISKY_FILE_PATH + hisky_maps_file)

    gsm_frequencies = pyfits.getdata(path + HISKY_FILE_PATH + hisky_freq_file)

    assert ctx.frequency >= gsm_frequencies[0], "Frequency (%s) outside available frequencies (%s - %s)"%(ctx.frequency, 
                                                                                                          gsm_frequencies[0], 
                                                                                                          gsm_frequencies[-1])
    assert ctx.frequency <= gsm_frequencies[-1], "Frequency (%s) outside available frequencies (%s - %s)"%(ctx.frequency, 
                                                                                                           gsm_frequencies[0], 
                                                                                                           gsm_frequencies[-1])
    
    
    for i, frequency in enumerate(gsm_frequencies):
        if ctx.frequency == frequency:
            break

    gsm_map = gsm_maps[i, :]        # == [i][:]

    return gsm_map

 
