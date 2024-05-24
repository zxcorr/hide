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
Created on 2021

authors: Joao Alberto, Bruno Ponquio

Update: May, 2024
authors: Thiago Pena
'''
from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import healpy as hp
from hope import jit
import zmod.zernike_fit as zfit
import astropy.io.fits as pyfits
from scipy import interpolate
import os
from pkg_resources import resource_filename
import hide

BEAM_FITS_PATH = "/data/1.123deg_fit/beta22/"

def load_beam_profile(beam_spec, frequencies, params):
        """
        Creates a 2d gaussian beam profile for the given params definition 
        
        :param params: The params instance with the paramterization
        
        :returns profile: A list of callable beam profiles
        """
        
        frequencies = frequencies/1000. #conv. to GHz
        beam_profiles = []
        beam_norms = []
        Npoints = [len(beam_spec.ra),len(beam_spec.dec)]

        fits_file = os.path.join(BEAM_FITS_PATH, params.zernike_coefficients_file_name)
        fits_path = resource_filename(hide.__name__, fits_file)#params.zernike_coefficients_file_name)
        with pyfits.open(fits_path) as zernike_fits:
                Coeffs_unsorted = zernike_fits[0].data
                zernike_freqs = np.array([f[0] for f in zernike_fits[1].data])
                Coeffs = [x for _,x in sorted(zip(zernike_freqs, Coeffs_unsorted))]
                zernike_freqs.sort()
                R = float(zernike_fits[0].header["RADIUS"])

        assert (R<=np.amax(beam_spec.dec) and 
                        R<=np.amax(beam_spec.ra)), \
                        ("Beam support must be greater than "
                        "the beam radius ({} deg).".format(np.degrees(R)))
        assert (np.amin(zernike_freqs)<=np.amin(frequencies) and
                        np.amax(zernike_freqs)>=np.amax(frequencies)), \
                        ("Frequency range chosen is out of the frequency range given by the Zernike"
                        "coefficients file ({}-{} GHz).".format(np.amin(zernike_freqs),
                                                                                                        np.amax(zernike_freqs)))
        if set(frequencies).issubset(set(zernike_freqs)):
                matching_freqs = np.array([f for f in frequencies if (zernike_freqs==f).any()])
                Coeffs_interp = np.array(Coeffs)[zernike_freqs==matching_freqs]
        else:
                Coeffs_interp = coeffs_interpolation(zernike_freqs, Coeffs, frequencies)


        beam_profiles = [zernike_wrapper(R, Coeffs_interp[f_idx],
                                                                         beam_spec, params.beam_response,
                                                                         params)
                                         for f_idx in range(len(frequencies))] 
        beam_norms = [normalization() for _ in frequencies]
        return beam_profiles, beam_norms



def coeffs_interpolation(data_freqs, data_coeffs, new_freqs):

        Alphabeta_idxs = data_coeffs[0,:,1:3]
        # Saving the numerical values separetely and transposed, in which each row is a list of
        # fixed (alpha, beta) for all frequencies.
        coeffs_values = data_coeffs[:,:,0].T
        
        interp_wrappers = np.array([interpolate.interp1d(data_freqs, C_ab, kind="cubic")
                                                            for C_ab in coeffs_values])
        interp_coeffs_values = [interp_wrapper(new_freqs) 
                                                        for interp_wrapper in interp_wrappers]
        interp_coeffs_values = np.array(interp_coeffs_values).T.reshape(len(new_freqs),
                                                                                                                                        len(Alphabeta_idxs), 1)
                                                                                                                                        
        # Reconstructing the coefficients list back into the original shape.
        # shape = (n_freqs, n_coeffs, 3)
        interp_coeffs = np.array([np.hstack((coeff_f,Alphabeta_idxs)) 
                                                          for coeff_f in interp_coeffs_values])
        return interp_coeffs


# cumbersome call to avoid mem leak in hope
def zernike_wrapper(R, Coeffs, beam_spec, beam_response, params):
        ra, dec = map(np.ndarray.flatten, np.meshgrid(beam_spec.ra, beam_spec.dec))
        beam_0 = zfit.radec_beam(Coeffs,ra,dec,R).filled(fill_value=0)
        
        def wrapped(i,j): # (i,j) -> (ra,dec)
#               from matplotlib.path import Path
#               rect_verts = [(beam_spec.ra[0],beam_spec.dec[0]),
#                                         (beam_spec.ra[0],beam_spec.dec[-1]),
#                                         (beam_spec.ra[-1],beam_spec.dec[-1]),
#                                         (beam_spec.ra[-1],beam_spec.dec[0]),]
#               p = Path(rect_verts)
#               hp_points = np.vstack((i,j)).T
#               grid_mask = p.contains_points(hp_points)
#               valid_idxs = np.where(grid_mask==True)
#               valid_i, valid_j = i,j#i[valid_idxs], j[valid_idxs]

                i0, j0 = i,j
                #i, j = i+np.pi/2, j+np.pi/2
                beam = np.nan_to_num(interpolate.griddata((ra,dec), 
                                                                                                  beam_0,
                                                                                                  (i,j),
                                                                                                  method=params.interpolation_scheme))
                norm = np.sum(beam)#*hp.nside2pixarea(params.beam_nside, degrees=False)
                beam /= norm #2*norm*hp.nside2pixarea(params.beam_nside, degrees=False)
                
                plot=False
                fits_writing = False
                field_directions_file = False
                
                if field_directions_file:
                
                        output_path = "/home/joao/Documentos/cosmologia/grasp/hide_field_directions/teste_inicial/"
                
                        #uv_center = [-0.004,0.0475] #(30,-304)
                        uv_center = [0.1571,-0.0475] #(-990,304)
                        
                        # 1.(ra,dec) -> (theta,phi)
                        pointing_thetaphi = zfit.radec2thetaphi(np.vstack((i0,j0)).T)
                        
                        # 2.(u,v) centro -> (theta,phi) centro
                        thetaphi_center = zfit.uv2thetaphi(uv_center)
                        print("Beam center (theta,phi) = ", np.degrees(thetaphi_center))
                        
                        rotation_angles = [-thetaphi_center[0], thetaphi_center[1]]
                        # 3.(theta,phi) -> rot para (theta,phi) centro
                        rot_thetaphi = zfit.rotate_coords(pointing_thetaphi, 
                                                                                          rotation_angles, 
                                                                                          verbose=True)
                        
                        # 4.(theta,phi) rot -> (u,v)
                        pointing_u = np.sin(rot_thetaphi.T[0]) * np.cos(rot_thetaphi.T[1])
                        pointing_v = np.sin(rot_thetaphi.T[0]) * np.sin(rot_thetaphi.T[1])
                        pointing_uv = np.vstack((pointing_u,pointing_v)).T
                        
                        if False:
                                import matplotlib.pyplot as plt
                                
                                plt.figure(1)
                                plt.title("(i,j) original (ra-dec)")
                                plt.scatter(i0,j0)
                                
                                fig2 = plt.figure(2)
                                ax2 = fig2.add_subplot(111, projection="polar")
                                ax2.set_title("(i,j) converted to thetaphi")
                                ax2.scatter(pointing_thetaphi.T[1],pointing_thetaphi.T[0])
                                
                                fig3 = plt.figure(3)
                                ax3 = fig3.add_subplot(111, projection="polar")
                                ax3.set_title("Rotated thetaphi")
                                ax3.scatter(rot_thetaphi.T[1],rot_thetaphi.T[0])
                                plt.scatter(thetaphi_center[1], thetaphi_center[0], c="red")
                                
                                plt.figure(4)
                                plt.title("Rotated uv (hide uv grid around beam center)")
                                plt.scatter(pointing_uv.T[0],pointing_uv.T[1])
                                plt.scatter(uv_center[0],uv_center[1],c="red")
                                
                                plt.show()
                                exit()
                        
                        # 5. write .sta
                        import os
                        file_path = os.path.join(output_path,"field_storage.sta")
                        with open(file_path, "w") as directions_f:
                        
                                directions_f.write(str(len(pointing_uv))+"\n")
                                for uv_coords in pointing_uv:
                                        directions_f.write(" " + str(uv_coords[0]) +
                                                                           " " + str(uv_coords[1]) +
                                                                           " 0 0 1 0 0\n")
                                        # cx?:
                                        #directions_f.write(" " + str(uv_coords[0]) +
                                        #                                  " " + str(uv_coords[1]) +
                                        #                                  " 0 0 2 0 0\n")
                        from stat import S_IRWXO
                        os.chmod(file_path, S_IRWXO)
                        exit()
                                
                
                
                if plot:
                
                        base_name = "beam_nside{nside}_{{interp}}_pixarea_norm".format(nside=params.beam_nside)
                        if fits_writing:
                                import astropy.io.fits as pyfits
                                fitsfile = params.output_path + base_name + "_healpix.fits" #"/hide_beams/MINUS990_304/nside{nside}/healpix_MINUS990_304_beam_nside{nside}_{interp}.fits"
                                # 30_MINUS304
                        pixarea = hp.nside2pixarea(params.beam_nside, degrees=False)
                        import matplotlib.pyplot as plt

                        # Grid Plotting
                        plt.figure(0)
                        plt.scatter(ra,dec,c="blue",label="Pre Interp Grid")
                        plt.scatter(i,j,c="red",label="Post Interp Grid")
                        plt.legend()
                        plt.title("HIDE Beam Interpolation Grid")
                        plt.xlabel("Dec (rad)")
                        plt.ylabel("RA (rad)")
                        plt.savefig(params.output_path + base_name.format(interp="") + "grid.png")
                        
                        # Pre-Interp Beam
                        new_shape = (int(np.sqrt(ra.shape[0])), int(np.sqrt(ra.shape[0])))
                        plt.figure(1)
                        plt.pcolormesh(ra.reshape(new_shape),
                                                   dec.reshape(new_shape),
                                                   20*np.log10(abs(beam_0)).reshape(new_shape), 
                                                   shading="auto") #, vmin=0,vmax=50)
                        plt.title("HIDE Beam - Pre Interp Beam")
                        plt.xlabel("Dec (rad)")
                        plt.ylabel("RA (rad)")
                        plt.colorbar(label="Amplitude (dB)")
                        plt.savefig(params.output_path + base_name.format(interp="preinterp") + ".png")
                        
                        # Pre-Interp Normalised Beam
                        norm_beam_0 = beam_0/np.sum(beam_0)/pixarea
                        vmin = 20*np.log10(min(abs(norm_beam_0[np.nonzero(abs(norm_beam_0))])))
                        vmax = 20*np.log10(max(abs(norm_beam_0)))
                        plt.figure(2)
                        plt.pcolormesh(ra.reshape(new_shape),
                                       dec.reshape(new_shape),
                                       20*np.log10(abs(norm_beam_0)).reshape(new_shape), 
                                       #norm_beam_0.reshape(new_shape),
                                       shading="auto", vmin=vmin, vmax=vmax)
                        plt.title("HIDE Beam - Pre Interp Normalised Beam")
                        plt.xlabel("Dec (rad)")
                        plt.ylabel("RA (rad)")
                        plt.colorbar(label="Intensity")#"Amplitude (dB)")
                        plt.savefig(params.output_path + base_name.format(interp="preinterp_norm") + ".png")
                        
                        
                        reso = hp.nside2resol(params.beam_nside, arcmin=True)/10
                        xysize = 30*5*int(params.beam_nside/128)

#                       # Post-Interp Beam
#                       plt.figure(3)
#                       theta, phi = np.pi/2+i, j
#                       idxs = hp.ang2pix(params.beam_nside, theta,phi)
#                       mapa = np.zeros(12*params.beam_nside**2)
#                       mapa[idxs]=beam
#                       reso = hp.nside2resol(params.beam_nside, arcmin=True)/10
#                       hp.gnomview(20*np.log10(abs(mapa)), fig=23 cbar=True, unit="Amplitude (dB)",
#                                               title=("HIDE Beam - {}"
#                                               " Interpolation".format(params.interpolation_scheme)),
#                                               rot=(0,0), xsize=xysize, ysize=xysize,
#                                               reso=reso, min=vmin, max=vmax)
#                       hp.graticule(local=True) # default interval = 1 deg
                        
                        # Nearest Beam
                        beam = np.nan_to_num(interpolate.griddata((ra,dec), 
                                                                                                          beam_0,
                                                                                                          (i,j),
                                                                                                          method="nearest"))
                        norm = np.sum(beam)*pixarea
                        beam /= norm
                        plt.figure(4)
                        theta, phi = np.pi/2+i, j
                        idxs = hp.ang2pix(params.beam_nside, theta,phi)
                        mapa = np.zeros(12*params.beam_nside**2)
                        mapa[idxs]=beam
                        
                        if fits_writing:
                                pyfits.writeto(fitsfile.format(interp="nearest"), 
                                                           mapa,overwrite=True)
                        hp.gnomview(20*np.log10(abs(mapa)), fig=4, cbar=True, unit="Amplitude (dB)",
                                                title="HIDE Beam - Nearest Interpolation",
                                                rot=(0,0), xsize=xysize, ysize=xysize,
                                                reso=reso, min=vmin, max=vmax)
                        hp.graticule(local=True) # default interval = 1 deg
                        plt.savefig(params.output_path + base_name.format(interp="nearest") + ".png")
                        
                        # Linear Beam
                        beam = np.nan_to_num(interpolate.griddata((ra,dec), 
                                                                                                          beam_0,
                                                                                                          (i,j),
                                                                                                          method="linear"))
                        norm = np.sum(beam)*pixarea
                        beam /= norm
                        plt.figure(5)
                        theta, phi = np.pi/2+i, j
                        idxs = hp.ang2pix(params.beam_nside, theta,phi)
                        mapa = np.zeros(12*params.beam_nside**2)
                        mapa[idxs]=beam
                        
                        if fits_writing:
                                pyfits.writeto(fitsfile.format(interp="linear"),
                                                           mapa,overwrite=True)
                        hp.gnomview(20*np.log10(abs(mapa)), fig=5, cbar=True, unit="Amplitude (dB)",
                                                title="HIDE Beam - Linear Interpolation",
                                                rot=(0,0), xsize=xysize, ysize=xysize,
                                                reso=reso, min=vmin, max=vmax)
                        hp.graticule(local=True) # default interval = 1 deg
                        plt.savefig(params.output_path + base_name.format(interp="linear") + ".png")
                        
                        # Cubic Beam
                        beam = np.nan_to_num(interpolate.griddata((ra,dec), 
                                                                                                          beam_0,
                                                                                                          (i,j),
                                                                                                          method="cubic"))
                        norm = np.sum(beam)*pixarea
                        beam /= norm
                        plt.figure(6)
                        theta, phi = np.pi/2+i, j
                        idxs = hp.ang2pix(params.beam_nside, theta,phi)
                        mapa = np.zeros(12*params.beam_nside**2)
                        mapa[idxs]=beam
                        
                        if fits_writing:
                                pyfits.writeto(fitsfile.format(nside=params.beam_nside,interp="cubic"),
                                                           mapa,overwrite=True)
                        hp.gnomview(20*np.log10(abs(mapa)), fig=6, cbar=True, unit="Amplitude (dB)",
                                                title="HIDE Beam - Cubic Interpolation",
                                                rot=(0,0), xsize=xysize, ysize=xysize,
                                                reso=reso, min=vmin, max=vmax)
                        hp.graticule(local=True) # default interval = 1 deg
                        plt.savefig(params.output_path + base_name.format(interp="cubic") + ".png")
                        
                        exit()

                return beam
        return wrapped

def normalization():
        return 1

