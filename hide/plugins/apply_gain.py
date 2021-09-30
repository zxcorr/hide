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
Created on Feb 27, 2015

author: jakeret
'''
from __future__ import print_function, division, absolute_import, unicode_literals

from ivy.plugin.base_plugin import BasePlugin

import importlib

class Plugin(BasePlugin):
    """
    Transform the temperature based (Kelvin) TOD into ADU by applying a
    spectrometer specific gain
    """

    def __call__(self):
        #load module
        mod = importlib.import_module(self.ctx.params.instrument)
        
        GAIN_PATH = ctx.params.gain_path
        
        #delegate loading of gain per frequency
        gain = mod.get_gain(GAIN_PATH, self.ctx.frequencies)
        
        #apply gain to TOD
        #TODO: no gain for Y-polarization
#         self.ctx.gain = gain
        self.ctx.tod_vx *= gain.reshape(-1, 1)
    
    def __str__(self):
        return "Applying gain to TOD"
