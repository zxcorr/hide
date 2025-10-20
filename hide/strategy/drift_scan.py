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
Created on Jan 15, 2015

author: jakeret

Update: October, 2025
author: Nicolli Soares
'''
from __future__ import print_function, division, absolute_import, unicode_literals

from datetime import timedelta
import numpy as np

from hide.strategy import CoordSpec
from hide.utils import sidereal


def load_strategy(ctx):
	"""
	Creates a scanning strategy that uses drift mode i.e. the 
	telescope stares at the same position for 24 hours and then changes the altiude by a certain angle
	
	:param ctx: The ctx instance with the paramterization
	
	:returns strategy: A list of CoordSpec with the scanning strategy
	"""
	
	start = ctx.strategy_start
	end   = ctx.strategy_end
	
	diff = end-start 
	
	
	lat_lon = sidereal.LatLon(np.radians(ctx.params.telescope_latitude), 
							  np.radians(ctx.params.telescope_longitude))
	
	alt_delta = np.radians(ctx.params.alt_delta)
	az_pos = np.radians(ctx.params.azimuth_pointing)
	alt_pos = np.radians(ctx.params.altitude_start_pos)
	alt_max = np.radians(ctx.params.altitude_max_pos)
	
	strategy = []

	sec = 0.0

	step = float(ctx.params.strategy_step_size)
	step_us = int(step * 10**6) # step in microseconds
	
	duration = int(diff.total_seconds())
	duration_us = int(duration * 10**6) # duration in microseconds
	n_intervals = duration_us // step_us
	delta_time = timedelta(microseconds = step_us)
	sec_us_per_day = (24*60*60)*10**6
	
	for i in range(n_intervals+1):
		
		current_time = start + i*delta_time
		sec_us = i*step_us
		sec = sec_us/(10**6)
		
		gst   = sidereal.SiderealTime.fromDatetime(current_time)
		radec = sidereal.AltAz(alt_pos, az_pos).raDec(gst.lst(lat_lon.lon), lat_lon)
		strategy.append(CoordSpec(sec, alt_pos, az_pos, radec.ra, radec.dec))

		if sec_us % sec_us_per_day == 0 and sec_us != 0:
			alt_pos += alt_delta

		if alt_pos > alt_max:
			alt_pos = np.radians(ctx.params.altitude_start_pos)

	return strategy