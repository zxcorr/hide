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

author: jakeret
'''
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import healpy as hp
from hope import jit
import matplotlib.pyplot as plt

def load_beam_profile(beam_spec, frequencies, params):
	"""
	Creates a 2d gaussian beam profile for the given params definition 
	
	:param params: The params instance with the paramterization
	
	:returns profile: A list of callable beam profiles
	"""

	beam_profiles = []
	beam_norms = []
	for k, frequency in enumerate(frequencies):
		fwhm = params.fwhm_0 * params.beam_frequency_min / frequency
		sigma = fwhm2sigma(fwhm) #fwhm / (2. * np.sqrt(2. * np.log(2)))
		beam_norms.append(normalization(sigma, params.beam_nside))
		beam_profiles.append(gauss_wrapper(sigma, params.beam_response, params))
	return beam_profiles, beam_norms

# def _2d_gauss(sigma, beam_response):
#	 i2sigma2 = 1. / (2*sigma**2)
#	 def wrapped(i,j):
# #		 Z =  np.exp(-i**2/(2*sigma**2) - j**2/(2*sigma**2))
#		 Z =  np.exp(i2sigma2 * (-i**2 - j**2))
#		 return Z / Z.sum() * beam_response
#	 return wrapped


def fwhm2sigma(FWHM):
        return FWHM / (2.*np.sqrt(2*np.log(2)))

def sigma2fwhm(sigma):
        return sigma * 2.*np.sqrt(2*np.log(2))


# cumbersome call to avoid mem leak in hope
def gauss_wrapper(sigma, beam_response, params):
	i2sigma2 = 1. / (2*sigma**2)
	def wrapped(i,j):
		Z = np.empty_like(i)
		hope_gauss(i2sigma2, beam_response, i, j, Z)
		
		plot_beam = True
		fits_writting = True
		if plot_beam:
			base_name = "beam_nside{nside}_gauss".format(nside=params.beam_nside)
			if fits_writting:
				import astropy.io.fits as pyfits
				fitsfile = params.output_path + base_name + "_healpix.fits"
			
			reso = hp.nside2resol(params.beam_nside, arcmin=True)/10
			xysize = 30*5*int(params.beam_nside/128)
			beam = Z
			plt.figure(5)
			theta, phi = np.pi/2+i, j
			idxs = hp.ang2pix(params.beam_nside, theta,phi)
			mapa = np.zeros(12*params.beam_nside**2)
			mapa[idxs]=beam#*normalization(sigma, params.beam_nside)
			total = np.sum(mapa)#*normalization(sigma, params.beam_nside)
			if fits_writting:
				pyfits.writeto(fitsfile, mapa,overwrite=True)
			hp.gnomview(20*np.log(abs(mapa)), fig=5, cbar=True, unit="Amplitude (dB)", #20*np.log10(abs(mapa))
									title="HIDE Beam - Gaussian FWHM={} Sum={} N={}".format(sigma2fwhm(sigma),total,len(idxs)),
									rot=(0,0), xsize=xysize, ysize=xysize,
									reso=reso)#, min=vmin, max=vmax)
			hp.graticule(local=True) # default interval = 1 deg
			plt.savefig(params.output_path + base_name.format(interp="linear") + ".png")
			exit()
		
		return Z
	return wrapped

def normalization(sigma, nside):
	n = (4 * np.pi) * sigma * sigma
	pixarea = hp.nside2pixarea(nside, degrees=False)
	return pixarea/n

@jit
def hope_gauss(i2sigma2, beam_response, i, j, Z):
	Z[:] =  np.exp(i2sigma2 * (-i**2 - j**2))
	#Z /= np.sum(Z) * beam_response
