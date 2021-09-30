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
Created on Nov 10, 2015
Modified on August, 2018 by Lucas Olivari

authors: jakeret, Lucas Olivari
'''

from __future__ import division

import numpy as np

def noise_amplitude(delta_nu, t_sys):

    return t_sys / np.sqrt(delta_nu)

def thermal_noise_tod(sfreq, size):

    samples = size[1]

    # frequencies
    f = np.fft.fftfreq(samples, d = 1. / sfreq)
    
    # scaling factor for all frequencies 
    s_scale = abs(np.concatenate([f[f<0], [f[-1]]]))
    
    s_scale = 1. + np.zeros(s_scale.size)
    
    # scale random power + phase
    sr = s_scale * np.random.normal(size=len(s_scale))
    si = s_scale * np.random.normal(size=len(s_scale))
    if not (samples % 2): si[0] = si[0].real

    s = sr + 1J * si
    
    # for odd sample numbers, there is one less positive 
    ## freq than for even sample numbers
    s = np.concatenate([s[1-(samples % 2):][::-1], s[:-1].conj()])

    # time series
    y = np.fft.ifft(s).real

    return y

def color_noise_tod(alpha, fknee, beta, delta_nu, sfreq, size):

    n_nu = size[0]
    samples = size[1]

    w_freq_0 = 1. / (n_nu * delta_nu)	

    # frequencies
    f_time = np.fft.fftfreq(samples, d = 1. / sfreq)
    w_freq = np.fft.fftfreq(n_nu, d = delta_nu)

    # scaling factor for all frequencies 
    s_scale_time = abs(np.concatenate([f_time[f_time<0], [f_time[-1]]]))
    s_scale_freq = abs(np.concatenate([w_freq[w_freq<0], [w_freq[-1]]]))

    ps_scale_time = (fknee / s_scale_time)**(alpha/2.)
    ps_scale_freq = (w_freq_0 / s_scale_freq)**((1. - beta)/(2. * beta))		

    # Stuart's normalization	

    ps_scale_freq_norm = (w_freq_0 / s_scale_freq)**((1. - beta)/(beta))
    sum_norm = np.sum(w_freq_0 * ps_scale_freq_norm[:-1])
    int_norm = (n_nu - 1.)/(2. * n_nu * delta_nu) 		
    c_norm = n_nu / (1. + (sum_norm / int_norm) * (n_nu - 1.))	

    test = np.zeros((ps_scale_freq.size, samples), dtype=np.complex64)	

    for nu in range(0, ps_scale_freq.size):
        psr = ps_scale_time * np.random.normal(size=len(ps_scale_time))
        psi = ps_scale_time * np.random.normal(size=len(ps_scale_time))
        if not (samples % 2): psi[0] = psi[0].real
        psc = psr + 1J * psi
        psc = np.concatenate([psc[1-(samples % 2):][::-1], psc[:-1].conj()])
  
        test[nu, :] = ps_scale_freq[nu] * psc

    test_out = np.zeros((n_nu, samples), dtype=np.complex64)	

    for s in range(0, samples):
	test_out[:, s] = np.concatenate([(test[:, s])[1-(n_nu % 2):][::-1], (test[:, s])[:-1]]) #np.concatenate([(test[:, s])[1-(n_nu % 2):][::-1], (test[:, s])[:-1].conj()]) 
		
    # time series -- note the normalization
    y = np.fft.ifftn(test_out).real * np.sqrt(n_nu)

    return np.sqrt(c_norm) * y

def noisegen(beta=0, N=2**13):
    """
    Noise will be generated that has spectral densities that vary as powers of inverse frequency,
    more precisely, the power spectra P(f) is proportional to 1 / fbeta for beta >= 0. 
    When beta is 0 the noise is referred to white noise, when it is 2 it is referred to 
    as Brownian noise, and when it is 1 it normally referred to simply as 1/f noise 
    which occurs very often in processes found in nature. 
    
    The basic method involves creating frequency components which have a magnitude that is 
    generated from a Gaussian white process and scaled by the appropriate power of f. 
    The phase is uniformly distributed on 0, 2pi.
    
    from http://paulbourke.net/fractals/noise/
    
    :param beta:
    :param N: number of samples (can also be shape of array)
    
    :returns out: the sampled noise
    """
    assert np.all(beta>=0) and np.all(beta <= 3), "Beta must be between 0 and 3"
    
    if type(N) is tuple:
        nc, nn_ = N
        if nn_%2 == 1:
            nn = nn_ + 1
        else:
            nn = nn_
        real = np.zeros((nc, nn/2+1))
        imag = np.zeros((nc, nn/2+1))
        shape = (nc, nn/2)
        beta = np.atleast_1d(beta).reshape(-1,1)
    else:
        nn_ = N
        if nn_%2 == 1:
            nn = nn_ + 1
        else:
            nn = nn_
        real = np.zeros(nn/2+1)
        imag = np.zeros(nn/2+1)
        shape = nn/2
    
    freq = np.fft.helper.rfftfreq(nn)[1:]
    mag = np.power(freq, -beta / 2.0) * np.random.normal(0.0, 1.0, shape) # Note to self a number of years later, why "i+1"
    pha = np.random.uniform(0, 2 * np.pi, size = shape)
    real[...,1:] = mag * np.cos(pha)
    imag[...,1:] = mag * np.sin(pha)
    b = real + imag*1j

    out = np.fft.irfft(b, norm = "ortho")
    return out.real[...,:nn_]

# translated, not vectorized version
# def noisegen(N=2**13, beta=0, seed=42):
#     real = np.empty(N)
#     imag = np.empty(N)
# 
#     np.random.seed = seed
# 
#     real[0] = 0;
#     imag[0] = 0;
# 
#     for i in range(1,N/2):
#         mag = pow(i+1.0,-beta/2) * np.random.normal(0.0,1.0); # Note to self a number of years later, why "i+1"
#         pha = 2*np.pi * np.random.uniform()
#         real[i] = mag * np.cos(pha);
#         imag[i] = mag * np.sin(pha);
# 
#         real[N-i] =  real[i];
#         imag[N-i] = -imag[i];
# 
# 
#     imag[N/2] = 0;
#     b = real + imag*1j
# 
#     out = np.fft.ifft(b)
#     return out
