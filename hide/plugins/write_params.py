# SEEK is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# SEEK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with SEEK.  If not, see <http://www.gnu.org/licenses/>.


'''
Created on July 30, 2021

author: Joao Alberto, Carlos Otobone
'''
from __future__ import print_function, division, absolute_import, unicode_literals

import os
from collections import namedtuple

from ivy.plugin.base_plugin import BasePlugin

import numpy as np

def register_params(ctx):
    date = ctx.params.strategy_start.split("-")
    
    file_fmt = ctx.params.params_file_fmt
    file_date_fmt = date[0] + date[1] + date[2]

    file_name = file_fmt.format(file_date_fmt)
    file_path = os.path.join(ctx.params.output_path, date[0], date[1], date[2], file_name)
    
    config_file = open(ctx.params.script_filename, "r")
    p_file = open(file_path, "w")
    
    lines = config_file.readlines()
    for line in lines:
    	p_file.write(line)
        
	config_file.close()
    p_file.close()
	



class Plugin(BasePlugin):
    """
    Writes a txt file containing the ctx.params information.
    """

    def __call__(self):
        register_params(self.ctx)
        
    def __str__(self):
        return "Register Params"
