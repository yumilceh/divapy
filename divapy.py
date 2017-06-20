"""
Created on Feb 3, 2017

@author: Juan Manuel Acevedo Valle
@department: Automatic Control Department (Knowledge Engineering Research Group)
@institution: Universitat Polit\`{e}cnica de Catalunya
@email: jmavbpl@gmail.com

This code is based on the matlab code:
http://sites.bu.edu/guentherlab/files/2016/10/DIVAsimulink1.zip
DIVA Source Code
Speechlab - Boston University
"""

import os
import random
import numpy as np
from numpy import tanh, array, fft
from scipy.io import loadmat as loadmat
from numpy_groupies.aggregate_weave import aggregate
from voiceboxpy import glotlf
from scipy.io.wavfile import write
import time

import matplotlib.pyplot as plt

eps = np.finfo(np.float32).eps
global_noise = False

diva_output_scale = [100.0, 500.0, 1500.0, 3000.0]

class Object(object):
    pass


class Diva(object):
    """
    Python implementation of the Diva Synthesizer
    """

    def __init__(self):
        self.vt = None
        self.fmfit = None

        abs_path = os.path.dirname(os.path.abspath(__file__))

        self.diva_synth_vt_file = abs_path + '/data/vt_py.mat'
        self.diva_synth_fmfit_file = abs_path + '/data/fmfit_py.mat'

        self.diva_hanning_file = abs_path + '/data/hanning.mat'
        self.hanning = loadmat(self.diva_hanning_file)['h']

        vt = loadmat(self.diva_synth_vt_file)
        fmfit = loadmat(self.diva_synth_fmfit_file)

        keys = ['vt_scale', 'vt_base', 'vt_average', 'vt_box']
        for key in keys:
            vt[key] = array(vt[key])
        keys = ['fmfit_beta_fmt', 'fmfit_p', 'fmfit_beta_som', 'fmfit_mu', 'fmfit_isigma']
        for key in keys:
            fmfit[key] = array(fmfit[key])
        self.vt = vt
        self.fmfit = fmfit

    def get_audsom(self, art, scale = False):
        """
            Input:
                art n_samples x n_articulators(13) [ndarray]
            Output:
                aud n_samples x formant_frequencies(4)
                som n_samples x somato out (DIVA) (6)
                outline n_samples x max_outline_shape_dim
                af  n_saples x max_af_shape_dim
        """
        art = tanh(art)

        if len(art.shape) > 1:
            n_samples = art.shape[0]
            aud = np.zeros((n_samples,4))
            som = np.zeros((n_samples,8))
            outline = []
            af = []
            d = []
            for i in range(n_samples):
                aud_, som_, outline_, af_, d_ = self.get_sample(art[i,:])
                if scale:
                    aud[i,:] = np.divide(aud_, diva_output_scale)
                else:
                    aud[i,:] = aud_
                som[i,:] = som_
                outline += [outline_]
                af += [af_]
                # d += [d_]
        else:
            aud, som, outline, af, d = self.get_sample(art)
            if scale:
                aud = np.divide(aud, diva_output_scale)

        return aud, som, outline, af

    def get_sound(self, art):
        """
            Input:
                art n_samples x n_articulators(13)
            Output:
                s sound signal from series of articulations art
        """
        art = tanh(art)
        synth = Object()
        synth.fs = 11025.
        synth.update_fs = 200.  # Sample frequency
        synth.f0 = 120.
        synth.samplesperperiod = np.ceil(synth.fs / synth.f0)

        glt_in = array(np.arange(0, 1, 1 / synth.samplesperperiod))
        synth.glottalsource = glotlf(0, glt_in)

        synth.f = array([0, 1])
        synth.filt = array([0, 0])
        synth.pressure = 0.

        synth.voicing = 1.
        synth.pressurebuildup = 0.
        synth.pressure0 = 0.
        synth.sample = np.zeros((synth.samplesperperiod,))#VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
        synth.k1 = 1
        synth.numberofperiods = 1
        synth.samplesoutput = 0

        self.vt['idx'] = range(10)
        self.vt['pressure'] = 0.
        self.vt['f0'] = 120
        self.vt['closed'] = False
        self.vt['closure_time'] = 0.
        self.vt['closure_position'] = 0.
        self.vt['opening_time'] = 0.

        voices = [{'F0': 120., 'size': 1.}, {'F0': 340., 'size': 0.7}]

        opt = Object()
        opt.voices = 0

        ndata = art.shape[0]
        dt = 0.005
        time = 0.
        s = np.zeros((np.ceil((ndata + 1) * dt * synth.fs),))#VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future

        while time < (ndata) * dt:
            t0 = np.floor(time / dt)
            t1 = (time - t0 * dt) / dt

            dummy1, dummy2, dummy3, af1, d = self.get_sample(art[np.minimum(ndata - 1, 0 + t0), :])#VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
            dummy1, dummy2, dummy3, af2, d = self.get_sample(art[np.minimum(ndata - 1, 1 + t0), :])#_VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future

            naf1 = len(af1)
            naf2 = len(af2)
            if naf2 < naf1:
                af2 = np.concatenate((af2, np.repeat(af2[-1], naf1 - naf2)))
            if naf1 < naf2:
                af1 = np.concatenate((af1, np.repeat(af1[-1], naf2 - naf1)))
            af = af1 * (1 - t1) + af2 * t1

            FPV = art[np.minimum(ndata - 1, t0), -3:] * (1 - t1) + art[np.minimum(ndata - 1, 1 + t0), -3:] * t1#VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
            FPV = np.minimum(np.ones((1, 3)), FPV)
            FPV = np.maximum(-1 * np.ones((1, 3)), FPV)
            FPV = FPV.flatten()

            self.vt['voicing'] = (1 + np.tanh(3 * FPV[2])) / 2.
            self.vt['pressure'] = FPV[1]
            self.vt['pressure0'] = self.vt['pressure'] > 0.01
            self.vt['f0'] = 100 + 20 * FPV[0]

            af0 = np.maximum(0 * af, af)
            k = 0.025

            idx_af0 = np.where(np.logical_and(np.greater(af0, 0 * af0), np.less(af0, k * np.ones(af0.shape))))
            af0[idx_af0] = k * np.ones((len(af0),))
            minaf = np.min(af)
            minaf0 = np.min(af0)
            self.vt['af'] = af

            # ======================================================================
            #    Tracks place of articulation
            # ======================================================================

            if minaf0 == 0:
                release = 0
                self.vt['opening_time'] = 0
                self.vt['closure_time'] = self.vt['closure_time'] + 1  # Here
                self.vt['closure_position'] = np.where(af0 == 0)[0][-1]
                if not self.vt['closed']:
                    # closure = self.vt['closure_position']
                    pass
                else:
                    pass
                    # closure = False
                self.vt['closed'] = True
            else:
                if self.vt['closed']:
                    release = self.vt['closure_position']
                    # release_closure_time = self.vt['closure_time']
                else:
                    release = 0
                if self.vt['pressure0'] and not synth.pressure0:
                    self.vt['opening_time'] = 0
                self.vt['opening_time'] = self.vt['opening_time'] + 1
                self.vt['closure_time'] = 0
                self.vt['closure_position'] = np.argmin(af)
                # closure = 0
                self.vt['closed'] = False

            if release > 0:
                af = np.maximum(k * np.ones(af.shape), af)
                minaf = np.max((k, minaf))
                minaf0 = np.max((k, minaf0))

            if release > 0:
                self.vt['f0'] = (0.95 + 0. * random.random()) * voices[opt.voices][
                    'F0']
                synth.pressure = 0
            elif (self.vt['pressure0'] and not synth.pressure0):
                self.vt['f0'] = (0.95 + 0. * random.random()) * voices[opt.voices]['F0']
                synth.pressure = self.vt['pressure']
                synth.f0 = 1.25 * self.vt['f0']
                synth.pressure = 1.
            elif (not self.vt['pressure0'] and synth.pressure0 and not self.vt['closed']):
                synth.pressure = synth.pressure / 10.

            # ======================================================================
            #  Computes glottal source
            # ======================================================================

            synth.samplesperperiod = int(np.ceil(synth.fs / synth.f0))
            pp = [0.6, 0.2 - 0.1 * synth.voicing, 0.25]

            glt_in = array(
                np.arange(0., 1.0 - 1.0e-10, 1. / synth.samplesperperiod))  # arange sometime includes [,stop]
            synth.glottalsource = 10 * 0.25 * glotlf(0, glt_in, pp) + \
                                  10 * 0.025 * synth.k1 * glotlf(1, glt_in, pp)
            numberofperiods = int(synth.numberofperiods)

            # ======================================================================
            # Computes vocal tract filter
            # ======================================================================

            synth.filt, synth.f, synth.filt_closure = \
                self.a2h(af0, d, synth.samplesperperiod, synth.fs, self.vt['closure_position'], minaf0)

            # synth.filt = synth.filt[:, 0]  # Change dimensions from (n,1) to (n,)
            # synth.filt_closure = synth.filt_closure[:, 0]  # Change dimensions from (n,1) to (n,)

            synth.filt = 2. * synth.filt / np.maximum(eps, synth.filt[0])
            synth.filt[0] = synth.filt[0] * 0.
            synth.filt_closure = 2. * synth.filt_closure / np.maximum(eps, synth.filt_closure[0])
            synth.filt_closure[0] = synth.filt_closure[0] * 0

            # ======================================================================
            #  Computes sound signal
            # ======================================================================

            w = np.linspace(0., 1., synth.samplesperperiod)
            if release > 0:
                u = synth.voicing * 1. * 0.010 * (synth.pressure + 20 * synth.pressurebuildup) \
                    * synth.glottalsource + (1 - synth.voicing) * 1. * 0.010 \
                                            * (synth.pressure + 20. * synth.pressurebuildup) \
                                            * array([0. * random.random() for i in range(synth.samplesperperiod)])

                v0 = np.real(fft.ifft(np.multiply(fft.fft(u), synth.filt_closure)))

                numberofperiods = numberofperiods - 1
                synth.pressure = synth.pressure / 10.
                vnew = v0[range(synth.samplesperperiod)]
                v0 = np.multiply(array([1 - x for x in w]), \
                                 synth.sample[(np.ceil(len(synth.sample) \
                                                       * array(range(synth.samplesperperiod)) \
                                                       / synth.samplesperperiod)).astype(int)]) + \
                     np.multiply(w, vnew)
                synth.sample = vnew

            else:
                v0 = array([])

            if numberofperiods > 0:
                u = 0.25 * synth.pressure * np.multiply(synth.glottalsource, array(
                    [0.0 * random.random() + 1 for i in range(synth.samplesperperiod)]))

                u = synth.voicing * u + (1 - synth.voicing) * synth.pressure * array(
                    [0.0 * random.random() for i in range(synth.samplesperperiod)])

                if minaf0 > 0 and minaf0 <= k:
                    u = minaf / k * u + (1 - minaf / k) * 0.02 * synth.pressure * array(
                        [0.0 * random.random() for i in range(synth.samplesperperiod)])

                v = np.real(fft.ifft(np.multiply(fft.fft(u), synth.filt)))

                vnew = v[0:synth.samplesperperiod]
                v1_ = [1 - x for x in w]
                v2_idx = [int(np.ceil(float(len(synth.sample) * float(x) / float(synth.samplesperperiod)))) - 1 for x in
                          range(synth.samplesperperiod)]

                while v2_idx[0] < 0:
                    last = v2_idx[0] + len(synth.sample)
                    v2_idx[:-1] = v2_idx[1:]  # might be slow
                    v2_idx[-1] = last
                v = np.multiply(v1_, synth.sample[v2_idx]) + np.multiply(w, vnew)

                synth.sample = vnew

                if numberofperiods > 1:
                    v = np.concatenate((v, np.repeat(vnew, numberofperiods - 1)))

            else:
                v = array([])

            v = np.concatenate((v0, v))
            # print(len(v))
            v = v + [.0000 * random.random() for i in range(len(v))]
            v = np.divide([1 - np.exp(-x) for x in v], [1 + np.exp(-x) for x in v])

            try:
                s[[synth.samplesoutput + x for x in range(len(v))]] = v
            except IndexError:
                s = np.concatenate((s, np.zeros(synth.samplesoutput + len(v) - len(s), )))
                s[[synth.samplesoutput + x for x in range(len(v))]] = v

            time = time + len(v) / synth.fs
            synth.samplesoutput = synth.samplesoutput + len(v)

            # ======================================================================
            #  Computes f0/amp/voicing/pressurebuildup modulation
            # ======================================================================

            synth.pressure0 = self.vt['pressure0']
            alpha = np.min((1, (0.1) * synth.numberofperiods))
            beta = 100. / synth.numberofperiods
            synth.pressure = synth.pressure + alpha * (
                self.vt['pressure'] * (np.max((1, 1.5 - self.vt['opening_time'] / beta))) - synth.pressure)
            alpha = np.min((1, .5 * synth.numberofperiods))
            beta = 100. / synth.numberofperiods
            synth.f0 = synth.f0 + 2 * np.square(alpha) * 0. * random.random() + alpha * (
                self.vt['f0'] * np.max((1., 1.25 - self.vt['opening_time'] / beta)) - synth.f0)  # %147;%120;
            synth.voicing = np.max((0, np.min((1, synth.voicing + 0.5 * (self.vt['voicing'] - synth.voicing)))))
            alpha = np.min((1, 0.1 * synth.numberofperiods))
            synth.pressurebuildup = np.max((0, np.min(
                (1, synth.pressurebuildup + alpha * (2 * float(self.vt['pressure'] > 0 and minaf < 0) - 1)))))
            synth.numberofperiods = np.max((1, numberofperiods))

        s = s[0:int(np.ceil(synth.fs * ndata * dt))]
        return s#, af

    def get_sample(self, art):
        """
            Input:
                art 1 x n_articulators(13)
                    art(1:10) vocaltract shape params
                    art(11:13) F0/P/V params
            Output:
                aud n_samples x [pitch(1) formant_frequencies(3)]
                som n_samples x somato out (DIVA) (6)
                    som(1:6) place of articulation (~ from pharyngeal to labial closure)
                    som(7:8) P/V params (pressure,voicing)
                outline n_samples x max_outline_shape_dim
                af  n_saples x max_af_shape_dim
                d
        """
        # ======================================================================
        # Computes vocal tract configuration
        # ======================================================================
        idx = range(10)

        x = np.multiply(self.vt['vt_scale'][idx].flatten(), art[idx])
        outline = self.vt['vt_average'] + np.dot(self.vt['vt_base'][:, idx].reshape((396, 10), order='F'),
                                                 x.reshape((10, 1)))
        outline = outline.flatten()

        # =======================================================================
        #  Computes somatosensory output (explicitly from vocal tract configuration)
        # =======================================================================

        som = np.zeros((8,))
        a, b, sc, af, d = self.xy2ab(outline)
        som[0:6] = np.maximum(-1 * np.ones((len(sc),)), np.minimum(np.ones((len(sc),)), -tanh(1 * sc)))
        som[6:] = art[-2:]

        aud = np.zeros((4,))
        aud[0] = 100 + 50 * art[-3]
        dx = art[idx] - self.fmfit['fmfit_mu']
        p = -1 * np.sum(np.multiply(np.dot(dx, self.fmfit['fmfit_isigma']), dx), axis=1) / 2
        p = np.multiply(self.fmfit['fmfit_p'].flatten(), np.exp(p - (np.max(p) * np.ones(p.shape))))
        p = p / np.sum(p)
        px = np.dot(p.reshape((len(p), 1)), np.append(art[idx], 1).reshape((1, len(idx) + 1)))
        aud[1:4] = np.dot(self.fmfit['fmfit_beta_fmt'], px.flatten(1))

        ####################### Not implemente yet (nargout = 2 or 3)
        # ---------------------------- if ~isempty(fmfit)&&nargout>1&&nargout<=3,
        # ------------------------------------ Som(1:6)=fmfit.beta_som*px(:)
        # ------------------------------------------ Som(7:8)=art(end-1:end)
        # ------------------------------------------------------------------- end

        return aud, som, outline, af, d

    def xy2ab(self, x, y=False):
        # x -> Outline (column matrix)
        if not hasattr(self, 'ab_alpha'):
            amax = 220
            alpha = array([1, 1, 1, 1, 1, 1, 1])
            beta = array([.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25])
            idx = [range(60), range(60, 70), range(70, 80), range(80, 120), range(120, 150), range(150, 190),
                   range(190, amax)]
            ab_alpha = np.zeros((amax,))
            ab_beta = np.zeros((amax,))
            for n1 in range(len(idx)):
                ab_alpha[idx[n1]] = alpha[n1]
                ab_beta[idx[n1]] = beta[n1]

            h = self.hanning.flatten()
            # h = hanning(51)/np.sum(hanning(51)) #Not same result as in hanning (matlab hann vs hanning)
            idx_2 = np.zeros((25,))
            idx_2 = np.concatenate((idx_2, np.array(range(amax))))
            idx_2 = np.concatenate((idx_2, (amax - 1) * np.ones((25,))))
            idx_2 = idx_2.tolist()

            ab_alpha = np.convolve(ab_alpha[idx_2], h, 'valid')#VisibleDeprecationWarning: non integer (and non boolean) array-likes will not be accepted as indices in the future
            ab_beta = np.convolve(ab_beta[idx_2], h, 'valid')#VisibleDeprecationWarning: non integer (and non boolean) array-likes will not be accepted as indices in the future
            self.ab_alpha = ab_alpha
            self.ab_beta = ab_beta

        if not y:
            i = np.array([0 + 1j])[0]
            x = np.exp(-i * np.pi / 12) * x
            y = np.imag(x)
            x = np.real(x)

        ab_alpha = self.ab_alpha
        ab_beta = self.ab_beta
        # =======================================================================
        # Grid
        # =======================================================================
        x0 = 45  # %90
        y0 = -100  # %-60
        r = 60  # %30
        k = np.pi * (r / 2)
        d = 0.75 / 10  # %unitstocm

        a = np.zeros(x.shape)
        b = np.zeros(x.shape)
        i1 = np.where(np.less(y, y0))
        i2 = np.where(np.less(x, x0))
        i3 = np.where(np.logical_and(np.greater_equal(y, y0), np.greater_equal(x, x0)))

        # =======================================================================
        # % a,b: "linearized" coordinates along vocal tract
        # =======================================================================

        a[i1] = y[i1] - y0
        b[i1] = x[i1] - x0
        a[i2] = k + x0 - x[i2]
        b[i2] = y[i2] - y0
        z = x[i3] - x0 + i * (y[i3] - y0)
        a[i3] = r * np.angle(z)
        b[i3] = np.abs(z)
        # =======================================================================
        # Tube area
        # =======================================================================

        olips = range(29, 45)
        ilips = range(256, 302)
        # ------------------------------------------------- owall = range(44,164)  #Not used in Matlab
        iwall = range(163 + 10, 257)
        oall = range(29, 164)
        iall = range(163, 302)
        xmin = -20
        ymin = -160
        amin = ymin - y0
        amax = np.ceil((x0 - xmin + k - amin))  # here is the variability of the af vector

        fact = 3

        idx_wallab1 = np.int_(np.maximum(np.zeros((len(oall),)), np.minimum(((fact * 9) - 1) * np.ones((len(oall),)),
                                                                            np.ceil(fact * 9 * (
                                                                                a[oall] - amin) / amax) - 1)))
        idx_wallab2 = np.int_(np.maximum(np.zeros((len(iwall),)), np.minimum(((fact * 9) - 1) * np.ones((len(iwall),)),
                                                                             np.ceil(fact * 9 * (
                                                                                 a[iwall] - amin) / amax) - 1)))
        wallab1 = aggregate(idx_wallab1, b[oall], size=fact * 9, func='min', fill_value=None)
        wallab2 = aggregate(idx_wallab2, b[iwall], size=fact * 9, func='max', fill_value=None)

        lipsab1 = np.nanmin(b[olips])
        lipsab2 = np.nanmax(b[ilips])

        mind_precursor = wallab1[range(fact * 2, fact * 8)] - wallab2[range(fact * 2, fact * 8)]
        mind = np.nanmin(mind_precursor.reshape((fact, 6), order='F'), axis=0)
        sc = mind[range(4)]
        sc = np.append(sc, np.nanmin(mind[range(4, 6)]))
        sc = np.append(sc, lipsab1 - lipsab2)
        sc = d * sc  # In Matlab this is a column vector

        w = 2

        idx_ab1 = np.int_(np.maximum(np.zeros((len(oall),)),
                                     np.minimum((amax - 1) * np.ones((len(oall),)), np.round(a[oall] - amin - 1))))
        idx_ab2 = np.int_(np.maximum(np.zeros((len(iall),)),
                                     np.minimum((amax - 1) * np.ones((len(iall),)), np.round(a[iall] - amin - 1))))

        ab1 = aggregate(idx_ab1, b[oall], size=amax, func='min', fill_value=None)
        ab2 = aggregate(idx_ab2, b[iall], size=amax, func='max', fill_value=None)

        ab1[np.isnan(ab1)] = np.inf
        ab2[np.isnan(ab2)] = -np.inf
        for n1 in range(w):
            ab1[1:-1] = np.minimum(np.minimum(ab1[1:-1], ab1[0:-2]), ab1[2:])
        for n1 in range(w):
            ab2[1:-1] = np.maximum(np.maximum(ab2[1:-1], ab2[0:-2]), ab2[2:])

        ab1[np.isinf(ab1)] = np.NINF
        idx_af = np.where(np.logical_and(np.greater(ab1, 0), np.greater(ab2, 0)))[0]
        af = d * (ab1[idx_af] - ab2[idx_af])
        idx = None
        for ii in range(len(af)):
            if af[ii] > 0:
                idx = ii
                break

        # =======================================================================
        # af: area function
        # =======================================================================
        af_tmp = np.minimum(np.zeros((len(af),)), af)
        af_power = np.power(np.maximum(np.zeros((len(af),)), af), ab_beta[idx_af])
        af = af_tmp + np.multiply(ab_alpha[idx_af],
                                  af_power)

        af = af[idx:]
        for ii in range(len(af) - 1, -1, -1):
            if np.isinf(af[ii]):
                af = np.delete(af, -1)
            else:
                break
        return a, b, sc, af, d

    # def a2h_(self, a, l, n, fs=None, closure=None, mina=None):  # To be deleted (Original function)
    #     if not isinstance(mina, float):
    #         mina = np.amin(a, axis=0)
    #     if not isinstance(closure, float) and not isinstance(closure, int):
    #         closure = 0
    #     if not isinstance(fs, float) and not isinstance(fs, int):
    #         fs = 11025.
    #     # if np.sum(len(np.where(a.shape>1))) <= 1:
    #     if not isinstance(a, np.ndarray):
    #         # a = a.flatten('F')
    #         N, M = (1, 1)
    #     else:
    #         if len(a.shape) == 1:
    #             N = a.shape[0]
    #             M = 1
    #         else:
    #             N, M = a.shape
    #
    #     # if not np.sum(len(np.where(l.shape>1))) <= 1:
    #     if not isinstance(l, np.ndarray):
    #         # l = l.flatten('F')
    #         NL, ML = (1, 1)
    #     else:
    #         NL, ML = l.shape
    #
    #     c = 34326.  # % speed of sound (cm/s)
    #
    #     m = int(np.ceil(n / 2.) + 1)
    #     f = fs * array(range(int(np.ceil(n / 2.) + 1))) / float(n)
    #     t = l / c
    #     Rrad = 0.9 * np.exp(-(np.abs(f) / 4.0e3) ** 2)  # % reflection at lips (low pass)
    #     H = np.zeros((m, M)) + 0j * np.zeros((m, M))  # %H=zeros(n,M);
    #     Hc = np.zeros((m, M)) + 0j * np.zeros((m, M))
    #     if mina == 0:
    #         a = array([np.max((0.05, x)) for x in a])
    #     # %k=0.995;%.999;
    #     for nM in range(M):
    #         # if 1: #%mina(nM)>0,
    #         coswt = np.cos(2 * np.pi * f * sub_array(t, idx_y=np.min((ML, nM)) - 1))
    #         sinwt = 1j * np.sin(2 * np.pi * f * sub_array(t, idx_y=np.min((ML, nM)) - 1))
    #         # R=[.9;(a(2:N,nM)-a(1:N-1,nM))./max(eps,a(2:N,nM)+a(1:N-1,nM))];
    #         R = np.concatenate(([0.9], np.divide(sub_array(a, idx_x=range(1, N, 1), idx_y=nM - 1) \
    #                                              - sub_array(a, idx_x=range(N - 1), idx_y=nM - 1), \
    #                                              np.maximum(eps * np.ones((N - 1,)),
    #                                                         sub_array(a, idx_x=range(1, N, 1), idx_y=nM - 1) \
    #                                                         + sub_array(a, idx_x=range(N - 1), idx_y=nM - 1)))))
    #
    #         U = coswt + sinwt
    #         V = coswt - sinwt
    #
    #         h1 = np.ones((m,))  # % signal at glottis   #In this part of the code there are diferences w.r.t. MATLAB,
    #         h2 = np.zeros((m,))
    #         for nN in range(N - 1):
    #             RnN = -R[nN]
    #             u = h1 + RnN * h2
    #             v = h2 + RnN * h1  # % reflection
    #             if closure == nN:
    #                 Hc[:, nM] = u - v
    #             if NL == 1:
    #                 h1 = np.multiply(U, u)
    #                 h2 = np.multiply(V, v)  # % delay
    #             else:  # This part of the code has not been tested
    #                 h1 = np.multiply(U[:, nN], u)
    #                 h2 = np.multiply(V[:, nN], v)
    #
    #         u = h1 - np.multiply(Rrad, h2)
    #         h = u
    #         if closure >= N - 1:  # Indexing here must be debugged
    #             Hc[:, nM] = u - (h2 - np.multiply(Rrad, h1))
    #
    #         H[:, nM] = np.divide(np.multiply((np.ones((Rrad.shape)) + Rrad), np.prod(np.ones((R.shape)) + R)),
    #                              sub_array(u, idx_y=0))
    #         if closure > 0 - 1:  # Hot point for indexing
    #             Hc[:, nM] = np.divide(np.multiply((np.ones((Rrad.shape)) + Rrad), \
    #                                               np.prod(np.ones((N - closure - 1,)) + R[closure + 1:N]) * Hc[:, nM]), \
    #                                   sub_array(h, idx_y=0))
    #
    #     idxH_x = [x + 1 for x in range(n - m - 1, -1, -1)]
    #
    #     H = np.concatenate((H, np.conjugate(H[idxH_x, :])), axis=0)
    #
    #     Hc = np.concatenate((Hc, np.conjugate(Hc[idxH_x, :])), axis=0)
    #
    #     f = np.concatenate((f, -f[idxH_x]))
    #     if mina == 0:
    #         H = 0 * H
    #
    #     return H, f, Hc%

    def a2h(self, a, l, n, fs=None, closure=None, mina=None):  # New function
        if not isinstance(mina, float):
            mina = np.amin(a, axis=0)
        if not isinstance(closure, float) and not isinstance(closure, int):
            closure = 0
        if not isinstance(fs, float) and not isinstance(fs, int):
            fs = 11025.
        # if np.sum(len(np.where(a.shape>1))) <= 1:
        if not isinstance(a, np.ndarray):
            # a = a.flatten('F')
            N, M = (1, 1)
        else:
            if len(a.shape) == 1:
                N = a.shape[0]
                M = 1
            else:
                N, M = a.shape

        # if not np.sum(len(np.where(l.shape>1))) <= 1:
        if not isinstance(l, np.ndarray):
            # l = l.flatten('F')
            NL, ML = (1, 1)
        else:
            NL, ML = l.shape

        c = 34326.  # % speed of sound (cm/s)

        m = int(np.ceil(n / 2.) + 1)
        f = fs * array(range(int(np.ceil(n / 2.) + 1))) / float(n)
        t = l / c
        Rrad = 0.9 * np.exp(-(np.abs(f) / 4.0e3) ** 2)  # % reflection at lips (low pass)
        H = np.zeros((m, M)) + 0j * np.zeros((m, M))  # %H=zeros(n,M);
        Hc = np.zeros((m, M)) + 0j * np.zeros((m, M))
        if mina == 0:
            a = array([np.max((0.05, x)) for x in a])
        # %k=0.995;%.999;

        R = np.zeros((1 + N - 1,))

        for nM in range(M):
            # if 1: #%mina(nM)>0,
            coswt = np.cos(2 * np.pi * f * sub_array(t, idx_y=np.min((ML, nM)) - 1))
            sinwt = 1j * np.sin(2 * np.pi * f * sub_array(t, idx_y=np.min((ML, nM)) - 1))
            # R=[.9;(a(2:N,nM)-a(1:N-1,nM))./max(eps,a(2:N,nM)+a(1:N-1,nM))];
            R[0] = 0.9
            R[1:] = np.divide(sub_array(a, idx_x=range(1, N, 1), idx_y=nM - 1) \
                              - sub_array(a, idx_x=range(N - 1), idx_y=nM - 1), \
                              np.maximum(eps * np.ones((N - 1,)),
                                         sub_array(a, idx_x=range(1, N, 1), idx_y=nM - 1) \
                                         + sub_array(a, idx_x=range(N - 1), idx_y=nM - 1)))

            U = coswt + sinwt
            V = coswt - sinwt

            h1 = np.ones((m,))  # % signal at glottis   #In this part of the code there are diferences w.r.t. MATLAB,
            h2 = np.zeros((m,))
            for nN in range(N - 1):
                RnN = -R[nN]
                u = h1 + RnN * h2
                v = h2 + RnN * h1  # % reflection
                if closure == nN:
                    Hc[:, nM] = u - v
                if NL == 1:
                    h1 = np.multiply(U, u)
                    h2 = np.multiply(V, v)  # % delay
                else:  # This part of the code has not been tested
                    h1 = np.multiply(U[:, nN], u)
                    h2 = np.multiply(V[:, nN], v)

            u = h1 - np.multiply(Rrad, h2)
            h = u
            if closure >= N - 1:  # Indexing here must be debugged
                Hc[:, nM] = u - (h2 - np.multiply(Rrad, h1))

            H[:, nM] = np.divide(np.multiply((np.ones((Rrad.shape)) + Rrad), np.prod(np.ones((R.shape)) + R)),
                                 sub_array(u, idx_y=0))
            if closure > 0 - 1:  # Hot point for indexing
                Hc[:, nM] = np.divide(np.multiply((np.ones((Rrad.shape)) + Rrad), \
                                                  np.prod(np.ones((N - closure - 1,)) + R[closure + 1:N]) * Hc[:, nM]), \
                                      sub_array(h, idx_y=0))

        idxH_x = [x + 1 for x in range(n - m - 1, -1, -1)]

        H_tmp = np.zeros((len(H) + len(idxH_x),), dtype=np.complex_)
        H_tmp[0:len(H)] = H.flatten()
        H_tmp[len(H):] = np.conjugate(H[idxH_x, :]).flatten()

        Hc_tmp = np.zeros((len(Hc) + len(idxH_x),), dtype=np.complex_)
        Hc_tmp[0:len(Hc)] = Hc.flatten()
        Hc_tmp[len(Hc):] = np.conjugate(Hc[idxH_x, :]).flatten()

        f_tmp = np.zeros((len(f) + len(idxH_x),))
        f_tmp[0:len(f)] = f
        f_tmp[len(f):] = -f[idxH_x]
        if mina == 0:
            H_tmp = 0 * H_tmp

        return H_tmp, f_tmp, Hc_tmp

    def plot_outline(self, art, axes=None):
        return_fig = False
        if axes is None:
            fig, axes = plt.subplots(1,1)
            return_fig = True
        a,b,outline,d = self.get_audsom(art)

        axes.plot(np.real(outline),np.imag(outline))
        plt.axis('off')
        # plt.plot(-np.flipud(np.real(outline)), np.flipud(np.imag(outline)))
        if return_fig:
            return fig, axes
        else:
            return axes

    def get_static_sound(self, art, play=False, ts = 0.005, time_=0.4):
        n_samples = time_/ts
        arts  = np.array([list(art)] * int(n_samples))
        return self.get_sound(arts)

    def play_sound(self, sound):  # keep in mind that DivaMatlab works with ts=0.005
        import pyaudio
        self.pa = pyaudio.PyAudio()  # If pa and stream are not elements of the self object then sound does not play
        self.stream = self.pa.open(format=pyaudio.paFloat32,
                                   channels=1,
                                   rate=11025,
                                   output=True)
        self.stream.start_stream()
        self.stream.write(sound.astype(np.float32).tostring())
        # time.sleep(0.4)

    def release_audio_device(self):  # any sound in the buffer will be removed
        try:
            self.stream.close()
            self.pa.terminate()
        except:
            pass


def arr_minus_noise(arr, noise_magnitude=0.):
    if global_noise:
        return array([x - random.random() for x in arr])

def arr_plus_cte(arr, cte):
    return array([x + cte for x in arr])

def sub_array(x, idx_x=np.inf, idx_y=np.inf):
    # This function supports floats, array(n,1) or array(n > 1, m > 1)
    if isinstance(x, float):
        return x

    n_dims = len(x.shape)
    if n_dims == 1:
        if not idx_x == np.inf:
            x = x[idx_x]
        return x

    if not idx_x == np.inf:
        x = x[idx_x, :]
    if not idx_y == np.inf:
        x = x[:, idx_y]

    return x

    # print('time {}: sum_v0 {}'.format(time, sum(v0)))
