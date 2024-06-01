# modified after seisfwd of seispy [Xu and He, 2022]
# Under the GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
#  Everyone is permitted to copy and distribute verbatim copies
#  of this license document, but changing it is not allowed.

#                             Preamble

#   The GNU General Public License is a free, copyleft license for
# software and other kinds of works.

#   The licenses for most software and other practical works are designed
# to take away your freedom to share and change the works.  By contrast,
# the GNU General Public License is intended to guarantee your freedom to
# share and change all versions of a program--to make sure it remains free
# software for all its users.  We, the Free Software Foundation, use the
# GNU General Public License for most of our software; it applies also to
# any other work released this way by its authors.  You can apply it to
# your programs, too.

#   When we speak of free software, we are referring to freedom, not
# price.  Our General Public Licenses are designed to make sure that you
# have the freedom to distribute copies of free software (and charge for
# them if you wish), that you receive source code or can get it if you
# want it, that you can change the software or use pieces of it in new
# free programs, and that you know you can do these things.

#   To protect your rights, we need to prevent others from denying you
# these rights or asking you to surrender the rights.  Therefore, you have
# certain responsibilities if you distribute copies of the software, or if
# you modify it: responsibilities to respect the freedom of others.

#   For example, if you distribute copies of such a program, whether
# gratis or for a fee, you must pass on to the recipients the same
# freedoms that you received.  You must make sure that they, too, receive
# or can get the source code.  And you must show them these terms so they
# know their rights.

#   Developers that use the GNU GPL protect your rights with two steps:
# (1) assert copyright on the software, and (2) offer you this License
# giving you legal permission to copy, distribute and/or modify it.

#   For the developers' and authors' protection, the GPL clearly explains
# that there is no warranty for this free software.  For both users' and
# authors' sake, the GPL requires that modified versions be marked as
# changed, so that their problems will not be attributed erroneously to
# authors of previous versions.

#   Some devices are designed to deny users access to install or run
# modified versions of the software inside them, although the manufacturer
# can do so.  This is fundamentally incompatible with the aim of
# protecting users' freedom to change the software.  The systematic
# pattern of such abuse occurs in the area of products for individuals to
# use, which is precisely where it is most unacceptable.  Therefore, we
# have designed this version of the GPL to prohibit the practice for those
# products.  If such problems arise substantially in other domains, we
# stand ready to extend this provision to those domains in future versions
# of the GPL, as needed to protect the freedom of users.

#   Finally, every program is threatened constantly by software patents.
# States should not allow patents to restrict development and use of
# software on general-purpose computers, but in those that do, we wish to
# avoid the special danger that patents applied to a free program could
# make it effectively proprietary.  To prevent this, the GPL assures that
# patents cannot be used to render the program non-free.


import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.util import next_pow_2
from math import pi
from numpy.fft import fft, ifft

from obspy import Trace, Stream
class ValueError(Exception): ...

def gaussFilter(dt, nft, gaussian):
    df = 1.0 / (nft * dt)
    nft21 = 0.5 * nft + 1
    f = df * np.arange(0, nft21)
    w = 2 * pi * f

    gauss = np.zeros([nft, 1])
    gauss1 = np.exp(-0.25 * (w / gaussian) ** 2) / dt
    gauss1.shape = (len(gauss1), 1)
    gauss[0:int(nft21)] = gauss1
    gauss[int(nft21):] = np.flipud(gauss[1:int(nft21) - 1])
    gauss = gauss[:, 0]

    return gauss

def gfilter(x, nfft, gauss, dt):
    Xf = fft(x, nfft)
    Xf = Xf * gauss * dt
    xnew = ifft(Xf, nfft).real
    return xnew

def correl(R, W, nfft):
    x = ifft(fft(R, nfft) * np.conj(fft(W, nfft)), nfft)
    x = x.real
    return x

def phaseshift(x, nfft, dt, preonset):
    Xf = fft(x, nfft)
    shift_i = int(preonset / dt)
    p = 2 * pi * np.arange(1, nfft + 1) * shift_i / nfft
    Xf = Xf * np.vectorize(complex)(np.cos(p), -np.sin(p))
    x = ifft(Xf, nfft) / np.cos(2 * pi * shift_i / nfft)
    x = x.real
    return x

def deconv_iterative(uin, win, dt, nt=None, preonset=10, gaussian=2.0, 
                     itmax=400, minderr=0.001, phase='P'):
    """
    Iterative Deconvolution of Receiver Function (Ligorria & Ammon, 1999)

    :param uin: numerator (radial for PdS) R or Q component for the 
                response function
    :type uin: np.ndarray
    :param win: denominator (vertical component for PdS) Z or L component for 
                the source function
    :type uin: np.ndarray
    :param dt: sample interval (s)
    :type dt: float
    :param nt: number of samples
    :type nt: int
    :param preonset: Time until beginning of receiver function (s)
    :type preonset: float
    :param gaussian: width of gaussian filter
    :type gaussian: float
    :param itmax: max # iterations
    :type itmax: int
    :param minderr: Min change in error required for stopping iterations
    :type minderr: float

    Returns:
        RFI: receiver function
        :type RFI: np.ndarray
        rms: Root mean square error for predicting numerator after each iteration
        :type rms: np.ndarray
        it: number of iterations
        :type it: int
    """

    # print('Iterative Decon (Ligorria & Ammon):\n')
    if len(uin) != len(win):
        raise ValueError('The two input trace must be in same length')
    elif nt is None:
        nt = len(uin)
    else:
        pass

    rms = np.zeros(itmax)
    nfft = next_pow_2(nt)
    p0 = np.zeros(nfft)

    u0 = np.zeros(nfft)
    w0 = np.zeros(nfft)

    u0[0:nt] = uin
    w0[0:nt] = win

    gaussF = gaussFilter(dt, nfft, gaussian)
    # gaussF = _gauss_filter(dt, nfft, gaussian)

    u_flt = gfilter(u0, nfft, gaussF, dt)
    w_flt = gfilter(w0, nfft, gaussF, dt)

    wf = fft(w0, nfft)
    r_flt = u_flt

    powerU = np.sum(u_flt ** 2)

    it = 0
    sumsq_i = 1
    d_error = 100 * powerU + minderr
    maxlag = 0.5 * nfft
    # print('\tMax Spike Display is ' + str((maxlag) * dt))

    while np.abs(d_error) > minderr and it < itmax:
        rw = correl(r_flt, w_flt, nfft)
        rw = rw / np.sum(w_flt ** 2)

        if phase == 'P':
            i1 = np.argmax(np.abs(rw[0:int(maxlag) - 1]))
        else:
            i1 = np.argmax(np.abs(rw))
        amp = rw[i1] / dt

        p0[i1] = p0[i1] + amp
        p_flt = gfilter(p0, nfft, gaussF, dt)
        p_flt = gfilter(p_flt, nfft, wf, dt)

        r_flt = u_flt - p_flt
        sumsq = np.sum(r_flt ** 2) / powerU
        rms[it] = sumsq
        d_error = 100 * (sumsq_i - sumsq)

        sumsq_i = sumsq

        it = it + 1

    p_flt = gfilter(p0, nfft, gaussF, dt)
    p_flt = phaseshift(p_flt, nfft, dt, preonset)
    RFI = p_flt[0:nt]
    rms = rms[0:it - 1]

    return RFI, rms, it 
    

def deconv_waterlevel(uin, win, dt, preonset=10., wlevel=0.05, 
                      gaussian=2.0, normalize=False, phase='P'):
    """
    Frequency-domain deconvolution using waterlevel method.

    :param uin: numerator (radial for PdS) R or Q component for 
                the response function
    :type uin: np.ndarray
    :param win: denominator (vertical component for PdS) Z or L component 
                for the source function
    :type uin: np.ndarray
    :param dt: sample interval (s)
    :type dt: float
    :param preonset: Time shift before P arrival, defaults to 10.
    :type preonset: float, optional
    :param wlevel: Waterlevel to stabilize the deconvolution, defaults to 0.05
    :type wlevel: float, optional
    :param gaussian: Gauss factor, defaults to 2.0
    :type gaussian: float, optional
    :param normalize: If normalize the amplitude of the RF, defaults to False
    :type normalize: bool, optional

    Returns: 
        (rft, rms) receiver function and root mean square error for 
        predicting numerator after each iteration
        :type: (np.ndarray, float)
    """
    if uin.size != win.size:
        raise ValueError('The length of the \'uin\' must be same as the \'win\'')
    nt = uin.size
    nft = next_pow_2(nt)
    nfpts = nft / 2 + 1     # number of freq samples
    fny = 1. / (2.* dt);     # nyquist
    delf = fny / (0.5 * nft)
    freq = delf * np.arange(nfpts)
    w = 2 * pi * freq

    # containers
    # rff = np.zeros(nft); # Rfn in freq domain
    upf = np.zeros(nft); # predicted numer in freq domain

    # Convert seismograms to freq domain
    uf = fft(uin, nft)
    wf = fft(win, nft)

    # denominator
    df = wf * wf.conjugate()
    dmax = max(df.real);   

    # add water level correction
    phi1 = wlevel * dmax # water level
    # nwl = length( find(df<phi1) ) # number corrected
    df[np.where(df.real < phi1)[0]] = phi1
    gaussF = gaussFilter(dt, nft, gaussian)
    nf = gaussF * uf * wf.conjugate()

    # compute RF
    rff = nf / df
    
    # compute predicted numerator
    upf = rff * wf

    # add phase shift to RF
    w = np.append(w, -np.flipud(w[1:-1]))
    rff = rff * np.exp(-1j * w * preonset)

    # back to time domain
    rft = ifft(rff , nft)
    rft = rft[0:nt].real

    # compute the fit
    uf = gaussF * uf # compare to filtered numerator
    ut = ifft(uf, nft).real
    upt = ifft(upf, nft).real

    powerU = np.sum(ut[0:nt] ** 2)
    rms = np.sum((upt[0:nt] - ut[0:nt]) ** 2 )/powerU

    if normalize:
        gnorm = np.sum(gaussF) * delf * dt
        rft = rft.real / gnorm

    return rft, rms

def deconvolve(uin, win, dt, method='iter', **kwargs):
    """
    Deconvolution of earthwuake data to create receiver function.

    :param uin: numerator (radial for PdS) R or Q component for the 
                response function
    :type uin: np.ndarray
    :param win: denominator (vertical component for PdS) Z or L component 
                for the source function
    :type uin: np.ndarray
    :param dt: sample interval (s)
    :type dt: float
    :param method: method for deconvolution('iter' or 'water'), 
                    defaults to 'iter'
    :type method: str, optional

    Returns:
        RFI: receiver function
        :type RFI: np.ndarray
        rms: Root mean square error for predicting numerator after each iteration
        :type rms: np.ndarray

    """	
    if method.lower() == 'iter':
        return deconv_iterative(uin, win, dt, **kwargs)
    elif method.lower() == 'water':
        return deconv_waterlevel(uin, win, dt, **kwargs)
    else:
        raise ValueError('method must be \'iter\' or \'water\'')

class RFTrace(obspy.Trace):
    def __init__(self, data=..., header=None):
        super().__init__(data=data, header=header)

    @classmethod
    def deconvolve(cls, utr, wtr, method='iter', **kwargs):
        header = utr.stats.__getstate__()
        for key, value in kwargs.items():
            header[key] = value
        if method.lower() == 'iter':
            rf, rms, it = deconv_iterative(utr.data, wtr.data, utr.stats.delta, 
                                           **kwargs)
            header['rms'] = rms
            header['iter'] = it
        elif method.lower() == 'water':
            rf, rms = deconv_waterlevel(utr.data, wtr.data, utr.stats.delta, 
                                        **kwargs)
            header['rms'] = rms
            header['iter'] = np.nan
        else:
            raise ValueError('method must be \'iter\' or \'water\'')
        return cls(rf, header)
    
    
ei = 0+1j

def e_inverse(omega, rho, alpha, beta, p):
    """ 
    E_inverse (Aki & Richards, pp. 161, Eq. (5.71))

    """

    e_inv = np.zeros([4,4], dtype=complex)
    eta = np.sqrt(1.0/(beta*beta) - p*p)
    xi  = np.sqrt(1.0/(alpha*alpha) - p*p)
    bp = 1.0 - 2.0*beta*beta*p*p

    e_inv[0,0] = beta*beta*p/alpha
    e_inv[0,1] = bp/(2.0*alpha*xi)
    e_inv[0,2] = -p/(2.0*omega*rho*alpha*xi) * ei
    e_inv[0,3] = -1.0/(2.0*omega*rho*alpha) * ei
    e_inv[1,0] = bp / (2.0*beta*eta)
    e_inv[1,1] = -beta*p
    e_inv[1,2] = -1.0/(2.0*omega*rho*beta) * ei
    e_inv[1,3] = p/(2.0*omega*rho*beta*eta) * ei
    e_inv[2,0] = e_inv[0,0]
    e_inv[2,1] = - e_inv[0,1]
    e_inv[2,2] = - e_inv[0,2]
    e_inv[2,3] = e_inv[0,3]
    e_inv[3,0] = e_inv[1,0]
    e_inv[3,1] = - e_inv[1,1]
    e_inv[3,2] = - e_inv[1,2]
    e_inv[3,3] = e_inv[1,3]
    return e_inv


def propagator_sol(omega, rho, alpha, beta, p, z):
    """ 
    Propagator (Aki & Richards, pp. 398, Eq. (3) in Box 9.1)
    
    """
    
    p_mat = np.zeros([4,4], dtype=complex)
    beta2 = beta*beta
    p2 = p*p
    bp = 1.0 -2.0*beta2*p2
    eta = np.sqrt(1.0/(beta2) - p2)
    xi  = np.sqrt(1.0/(alpha*alpha) - p2)
    cos_xi = np.cos(omega*xi*z)
    cos_eta = np.cos(omega*eta*z)
    sin_xi = np.sin(omega*xi*z)
    sin_eta = np.sin(omega*eta*z)

    p_mat[0,0] = 2.0*beta2*p2*cos_xi + bp*cos_eta
    p_mat[0,1] = p*( bp/xi*sin_xi - 2.0*beta2*eta*sin_eta ) * ei
    p_mat[0,2] = (p2/xi*sin_xi + eta*sin_eta)/(omega*rho)
    p_mat[0,3] = p*(-cos_xi + cos_eta)/(omega*rho) * ei  
    p_mat[1,0] = p*( 2.0*beta2*xi*sin_xi - bp/eta*sin_eta ) * ei
    p_mat[1,1] = bp*cos_xi + 2.0*beta2*p2*cos_eta
    p_mat[1,2] = p_mat[0,3]
    p_mat[1,3] = (xi*sin_xi + p2/eta*sin_eta)/(omega*rho)
    p_mat[2,0] = omega*rho*( -4.0*beta2*beta2*p2*xi*sin_xi - bp*bp/eta*sin_eta )
    p_mat[2,1] = 2.0*omega*beta2*rho*p*bp*( cos_xi - cos_eta ) * ei
    p_mat[2,2] = p_mat[0,0]
    p_mat[2,3] = p_mat[1,0]
    p_mat[3,0] = p_mat[2,1]
    p_mat[3,1] = -omega*rho*( bp*bp/xi*sin_xi + 4.0*beta2*beta2*p2*eta*sin_eta  )
    p_mat[3,2] = p_mat[0,1]  
    p_mat[3,3] = p_mat[1,1]

    return p_mat

def haskell(omega, p, nl, ipha, alpha, beta, rho, h):

    i0 = 0
    e_inv = e_inverse(omega, rho[-1], alpha[-1], beta[-1], p)
    p_mat = propagator_sol(omega, rho[i0], alpha[i0], beta[i0], p, h[i0] )
    for i in range(i0+1, nl):
        p_mat2 = propagator_sol(omega, rho[i], alpha[i], beta[i], p, h[i])
        p_mat = np.matmul(p_mat2, p_mat)
    if nl > i0+1:
        sl = np.matmul(e_inv, p_mat)
    else:
        sl = e_inv
    denom = sl[2,0] * sl[3,1] - sl[2,1] * sl[3,0]
    if ipha >= 0:
        ur = sl[3,1] / denom
        uz = - sl[3,0] / denom
    else:
        ur = - sl[2,1] / denom
        uz = sl[2,0] / denom
    return ur, uz

def fwd_seis(rayp, dt, npts, ipha, alpha, beta, rho, h):
    nlay = h.size
    npts_max = next_pow_2(npts)
    ur_freq = np.zeros(npts_max, dtype=complex)
    uz_freq = np.zeros(npts_max, dtype=complex)
    nhalf = int(npts_max / 2 + 1)
    for i in range(1, nhalf):
        omg = 2*np.pi * i / (npts_max * dt)
        ur_freq[i], uz_freq[i] = haskell(omg, rayp, nlay, ipha, 
                                         alpha, beta, rho, h)
    ur = ifft(ur_freq).real[::-1]/npts_max
    uz = -ifft(uz_freq).real[::-1]/npts_max

    return ur[0:npts], uz[0:npts]

class synrf():
    """
    Class for calculating synthetic receiver function using the forward 
    modelling method.

    
    Example
    -------
    Initialize the synrf module:
    >>> from rfsed.synrf import synrf
    Define the Earth model parameters
    >>> import numpy as np
    >>> depth = np.array([2, 35, 77.5])
    >>> vp = np.array([2.0, 6.5, 8.045])
    >>> vs = np.array([1.36, 3.75, 4.485])
    >>> rho=np.array([2.72, 2.92, 3.3455])
    >>> preonset=5
    >>> n=2100
    >>> rayp=0.04
    >>> gaussian=1.25
    >>> delta=0.025
    Call the synrf class
    >>> Synth=synrf(depth, vp, vs, rho, rayp, dt=delta, npts=n, ipha=1)
    >>> Synth.run_fwd()
    >>> Synth.filter(freqmin=0.05, freqmax=1.25, order=2, zerophase=True)
    >>> rf_synth=Synth.run_deconvolution(pre_filt=[0.05, 1.25], 
                                        preonset=preonset, gaussian=gaussian)
    >>> trdata=(rf_synth[0]).data
    """

    def __init__(self, depth, vp, vs, rho, rayp, dt, npts=2500, ipha=1, 
                 filter=None) -> None:
        """
        Initialize the synthetic receiver function
        """
        self.depth = depth
        self.vp = vp
        self.vs = vs
        self.rho = rho
        self.dt = dt
        self.npts = npts
        self.thickness=[]
        # convert depths in the model to layer thicknesses
        for i in range(len(self.depth)):
            if i == 0:
                thick=self.depth[i]
                self.thickness.append(thick)
            else:
                thick=self.depth[i]-self.depth[i-1]
                self.thickness.append(thick)
        self.thickness=np.array(self.thickness)
        if not isinstance(rayp, (float, list, np.ndarray)):
            raise TypeError('The rayp should be in float, list and np.ndarray')
        if isinstance(rayp, float):
            self.rayp = [rayp]
        else:
            self.rayp = rayp
        self.ipha = ipha

    def run_fwd(self):
        """
        Forward modelling synthetic seismograms.

        Returns:
            Radial and vertical Seismograms in ``Obspy.Stream`` type.
        """

        self.rstream = Stream()
        self.zstream = Stream()
        for _, rayp in enumerate(self.rayp):
            ur, uz = fwd_seis(rayp, self.dt, self.npts, self.ipha,
                            self.vp, self.vs, self.rho,
                            self.thickness)
            tr = Trace(data=ur)
            tr.stats.delta = self.dt
            self.rstream.append(tr)
            tr = Trace(data=uz)
            tr.stats.delta = self.dt
            self.zstream.append(tr)

    def filter(self, freqmin, freqmax, order=2, zerophase=True):
        """
        Apply a bandpass filter on synthetic waveforms

        :param freqmin: minimum cut-off frequency
        :type freqmin: float
        :param freqmax: maximum cut-off frequency
        :type fremax: float
        :param order: order of the filter, defaults to 2
        type order: int 
        :param zerophase: whether use a zero-phase filter, by default True
        :type zerophase: bool
        

        Returns:
            Filtered radial and vertical Seismograms in ``Obspy.Stream`` type.
    
        """
        for st in [self.rstream, self.zstream]:
            st.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                      corners=order, zerophase=zerophase)

    def run_deconvolution(self, pre_filt=[0.05, 2], preonset=10, gaussian=2.0, 
                          **kwargs):
        """
        Deconvolve receiver function from the synthetic seismograms

        :param pre_filt: pre-filter for deconvolution, defaults to [0.05, 2]
        :type pre_filt: list
        :param preonset: time shift before P arrival (in seconds), defaults to 10
        :type preonset: int
        :param gaussian: width of gaussian filter, defaults to 2.0
        :type gaussian: float

        Returns:
        Deconvolved receiver function in ``Obspy.Stream`` type.

        """	
        if pre_filt is not None:
            self.filter(*pre_filt)
        rfstream = Stream()
        for i, _ in enumerate(self.rayp):
            rftr = RFTrace.deconvolve(self.rstream[i], self.zstream[i], 
                                      preonset=preonset, gaussian=gaussian, 
                                      **kwargs)
            rfstream.append(rftr)
        return rfstream