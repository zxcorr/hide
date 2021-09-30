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
Created on Dec 8, 2014
Modified on August, 2018 by Lucas Olivari
Modified on September, 2020 by Carlos Otobone and Joao Alberto


authors: jakeret, Lucas Olivari
'''
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np

from ivy.plugin.base_plugin import BasePlugin
from hide.utils.signal import noisegen, noise_amplitude, thermal_noise_tod, color_noise_tod
import importlib

class Plugin(BasePlugin):
    """
    Adds background noise to the time ordered data
    """

    def __call__(self):
                
        params = self.ctx.params
        size = self.ctx.tod_vx.shape
        if params.load_noise_template:
            mod = importlib.import_module(self.ctx.params.instrument)
            freq = self.ctx.frequencies
            #wn_scale, cn_amp, cn_beta = mod.get_noise_params(freq)
            wn_scale = mod.get_noise_params(self.ctx.params.noise_path, freq)
            #params.white_noise_scale = wn_scale
            #params.color_noise_amp = cn_amp
            #params.color_noise_beta = cn_beta
        
        else:
            mod = importlib.import_module(self.ctx.params.instrument)
            freq = self.ctx.frequencies
            wn_scale = 1.0

        #noise = get_noise(params.white_noise_scale, params.color_noise_amp, params.color_noise_beta, size)

        noise = get_noise_lucas(wn_scale, params.color_alpha, params.color_fknee, params.color_beta, params.sample_freq, params.temp_sys, params.delta_nu, size)
        
        self.ctx.tod_vx += noise
        #TODO: no noise for Y-polarization
#         self.ctx.tod_vy += noise
        
    def __str__(self):
        return "Add background noise"
    
def get_noise(scale, alpha, beta, size):
    wnoise = np.random.normal(scale=np.atleast_1d(scale).reshape(-1,1),
                              size=size)
    # only create colored noise if amplitude is greater than zero
    if np.any(alpha > 0):
        rnoise = noisegen(beta, size) * np.atleast_1d(alpha).reshape(-1,1)
        return wnoise + rnoise
    else:
        return wnoise

def get_noise_lucas(scale_adu, alpha, fknee, beta, sfreq, tempsys, deltanu, size):

    n_nu = size[0]
    samples = size[1]
    
    wnoise = np.zeros((n_nu, samples))

    for i in range(0, n_nu):
        wnoise[i, :] = thermal_noise_tod(sfreq, size)

    wnoise_norm = np.std(wnoise) # We must normalize our TOD

    wnoise = np.mean(scale_adu) * noise_amplitude(deltanu, tempsys) * (wnoise / wnoise_norm)
    
    rnoise = np.mean(scale_adu) * noise_amplitude(deltanu, tempsys) * (color_noise_tod(alpha, fknee, beta, deltanu, sfreq, size) / wnoise_norm)
    
    # it always simulate colored noise
    return wnoise + rnoise
