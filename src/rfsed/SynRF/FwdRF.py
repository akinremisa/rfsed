import numpy as np
import obspy
from obspy.io.sac import SACTrace
from obspy.signal.util import next_pow_2
from math import pi
from numpy.fft import fft, ifft, ifftshift
import matplotlib.pyplot as plt
from obspy import Trace, Stream
from scipy.interpolate import interp1d
from os.path import join, dirname, exists
from scipy.io import loadmat
from matplotlib.colors import ListedColormap
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

def phaseshift(x, nfft, dt, tshift):
    Xf = fft(x, nfft)
    shift_i = int(tshift / dt)
    p = 2 * pi * np.arange(1, nfft + 1) * shift_i / nfft
    Xf = Xf * np.vectorize(complex)(np.cos(p), -np.sin(p))
    x = ifft(Xf, nfft) / np.cos(2 * pi * shift_i / nfft)
    x = x.real
    return x

def deconit(uin, win, dt, nt=None, tshift=10, gaussian=2.0, itmax=400, minderr=0.001, phase='P'):
    """
    Created on Wed Sep 10 14:21:38 2014

    In:
    uin = numerator (radial for PdS)
    win = denominator (vertical component for PdS)
    dt = sample interval (s)
    nt = number of samples
    tshift = Time until beginning of receiver function (s)
    gaussian = width of gaussian filter
    itmax = max # iterations
    minderr = Min change in error required for stopping iterations

    Out:
    RFI = receiver function
    rms = Root mean square error for predicting numerator after each iteration

    @author: Mijian Xu @ NJU
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
    p_flt = phaseshift(p_flt, nfft, dt, tshift)
    RFI = p_flt[0:nt]
    rms = rms[0:it - 1]

    return RFI, rms, it 
    

def deconwater(uin, win, dt, tshift=10., wlevel=0.05, gaussian=2.0, normalize=False, phase='P'):
    """
    Frequency-domain deconvolution using waterlevel method.

    :param uin: R or Q component for the response function
    :type uin: np.ndarray
    :param win: Z or L component for the source function
    :type win: np.ndarray
    :param dt: sample interval in second 
    :type dt: float
    :param tshift: Time shift before P arrival, defaults to 10.
    :type tshift: float, optional
    :param wlevel: Waterlevel to stabilize the deconvolution, defaults to 0.05
    :type wlevel: float, optional
    :param gaussian: Gauss factor, defaults to 2.0
    :type gaussian: float, optional
    :param normalize: If normalize the amplitude of the RF, defaults to False
    :type normalize: bool, optional

    :return: (rf, rms) RF and final rms.
    :rtype: (np.ndarray, float)
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
    rff = rff * np.exp(-1j * w * tshift)

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

def deconvolute(uin, win, dt, method='iter', **kwargs):
    if method.lower() == 'iter':
        return deconit(uin, win, dt, **kwargs)
    elif method.lower() == 'water':
        return deconwater(uin, win, dt, **kwargs)
    else:
        raise ValueError('method must be \'iter\' or \'water\'')

class RFTrace(obspy.Trace):
    def __init__(self, data=..., header=None):
        super().__init__(data=data, header=header)

    @classmethod
    def deconvolute(cls, utr, wtr, method='iter', **kwargs):
        header = utr.stats.__getstate__()
        for key, value in kwargs.items():
            header[key] = value
        if method.lower() == 'iter':
            rf, rms, it = deconit(utr.data, wtr.data, utr.stats.delta, **kwargs)
            header['rms'] = rms
            header['iter'] = it
        elif method.lower() == 'water':
            rf, rms = deconwater(utr.data, wtr.data, utr.stats.delta, **kwargs)
            header['rms'] = rms
            header['iter'] = np.nan
        else:
            raise ValueError('method must be \'iter\' or \'water\'')
        return cls(rf, header)
    
    
ei = 0+1j

def e_inverse(omega, rho, alpha, beta, p):
    """ E_inverse (Aki & Richards, pp. 161, Eq. (5.71))

    Parameters
    ----------
    omega : _type_
        _description_
    rho : _type_
        _description_
    alpha : _type_
        _description_
    beta : _type_
        _description_
    p : _type_
        _description_
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
    """ propagator (Aki & Richards, pp. 398, Eq. (3) in Box 9.1)

    Parameters
    ----------
    omega : _type_
        _description_
    rho : _type_
        _description_
    alpha : _type_
        _description_
    beta : _type_
        _description_
    p : _type_
        _description_
    z : _type_
        _description_

    Returns
    -------
    _type_
        _description_
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


class SynRF():
    def __init__(self, depth, vp, vs, rho, rayp, dt, npts=2500, ipha=1, filter=None) -> None:
        """_summary_

        Parameters
        ----------
        depth : _type_
            np.array of depth
        vp : _type_
            np.array of vp
        vs : _type_
            np.array of vs
        rho : _type_
            np.array of rho
        rayp : _type_
            Ray-parameter in s/km
        dt : _type_
            Time interval
        npts : _type_
            samples of synthetic waveform
        ipha : _type_
            Specify incident wave 1 for P and -1 for S
        """
        self.depth = depth
        self.vp = vp
        self.vs = vs
        self.rho = rho
        self.dt = dt
        self.npts = npts
        self.thickness=[]
        for i in range(len(self.depth)):
            if i == 0:
                thick=self.depth[i]
                self.thickness.append(thick)
            else:
                if i==1:
                    thick=self.depth[i]-self.depth[i-1]
                    self.thickness.append(thick)
                else:
                    if i==2:
                        thick=self.depth[i]-self.depth[i-1]
                        self.thickness.append(thick)
                    else:
                        if i==3:
                            thick=self.depth[i]-self.depth[i-1]
                            self.thickness.append(thick)
                        else:
                            if i==4:
                                thick=self.depth[i]-self.depth[i-1]
                                self.thickness.append(thick)
                            else:
                                if i==5:
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
        """Forward modelling synthetic seismograms.

        ``SynSeis.rstream`` and ``SynSeis.zstream`` are generated as 
        radial and vertical Seismograms in ``Obspy.Stream`` type.
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
        """Apply a bandpass filter on synthetic waveforms

        Parameters
        ----------
        freqmin : float
            Minimum cut-off frequency
        freqmax : float
            maximum cut-off frequency
        order : int, optional
            Order of filter, by default 2
        zerophase : bool, optional
            whether use a zero-phase filter, by default True
        """
        for st in [self.rstream, self.zstream]:
            st.filter('bandpass', freqmin=freqmin, freqmax=freqmax,
                      corners=order, zerophase=zerophase)

    def run_deconvolution(self, pre_filt=[0.05, 2], shift=10, gaussian=2.0, **kwargs):
        if pre_filt is not None:
            self.filter(*pre_filt)
        rfstream = Stream()
        for i, _ in enumerate(self.rayp):
            rftr = RFTrace.deconvolute(self.rstream[i], self.zstream[i], tshift=shift,
                                       gaussian=gaussian, **kwargs)
            rfstream.append(rftr)
        return rfstream