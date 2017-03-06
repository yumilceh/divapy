'''
Created on Feb 3, 2017

@author: Juan Manuel Acevedo Valle
'''
import os

import random
import numpy.linalg  as linalg
import numpy as np      
from numpy import tanh, matrix, array
from scipy.io import loadmat as loadmat
from numpy_groupies.aggregate_weave import aggregate 
    
eps = np.finfo(np.float32).eps
global_noise = False

        
class Object(object):
    pass

class Diva(object):
    '''
    Python implementation of the Diva Synthesizer
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.vt = None
        self.fmfit = None
        
        abs_path = os.path.dirname(os.path.abspath(__file__))

        self.diva_synth_vt_file = abs_path + '/DIVA/vt_py.mat'
        self.diva_synth_fmfit_file = abs_path + '/DIVA/fmfit_py.mat'
        
        self.diva_hanning_file = abs_path + '/DIVA/hanning.mat'
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
        
    
    def get_audsom(self, Art):
        '''
             Art n_samples x n_articulators(13)
        '''
        Art = tanh(Art)
        
        if len(Art.shape)>1:
            n_samples = Art.shape[0]    
        else:
            Aud, Som, Outline, af, d = self.get_sample(Art)
        return Aud, Som, Outline, af
    
    def get_sound(self, art):
        art = tanh(art)
        synth = Object()
        synth.fs = 11025.
        synth.update_fs = 200. # Modify sample frequency
        synth.f0 = 120.
        synth.samplesperperiod = np.ceil(synth.fs/synth.f0)
        
        glt_in = array(np.arange(0, 1 ,1/synth.samplesperperiod))
        synth.glottalsource = glotlf(0, glt_in)
        
        synth.f = array([0,1])
        synth.filt = array([0,0])
        synth.pressure = 0.
        #-------------------------------------------------- %synth.modulation=1;
        synth.voicing = 1.
        synth.pressurebuildup = 0.
        synth.pressure0 = 0.
        synth.sample = np.zeros((synth.samplesperperiod,))
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
        
        voices = [{'F0': 120.,'size':1.}, {'F0': 340.,'size':0.7}]

        opt = Object()
        opt.voices = 0
        
        
        ndata = art.shape[0]
        dt = 0.005
        time = 0.
        s = np.zeros((np.ceil((ndata+1)*dt*synth.fs),))
        
        while time<(ndata)*dt:
            #---------------------------------- % sample articulatory parameters
            t0 = np.floor(time/dt)
            print(t0)
            t1 = (time-t0*dt)/dt
            dummy1,dummy2,dummy3,af1,d = self.get_sample(art[np.minimum(ndata-1,0+t0),:])
            dummy1,dummy2,dummy3,af2,d = self.get_sample(art[np.minimum(ndata-1,1+t0),:])
            #print(np.minimum(ndata-1,1+t0))
            naf1 = len(af1)
            naf2 = len(af2)
            if naf2<naf1:
                af2 = np.concatenate((af2, np.repeat(af2[-1], naf1-naf2)))
            if naf1<naf2:
                af1 = np.concatenate((af1, np.repeat(af1[-1], naf2-naf1)))
            af = af1*(1-t1) + af2*t1
            
            
            FPV = art[np.minimum(ndata-1,t0),-3:]*(1-t1) + art[np.minimum(ndata-1,1+t0),-3:]*t1
            FPV = np.minimum(np.ones((1,3)), FPV)
            FPV = np.maximum(-1*np.ones((1,3)), FPV)
            FPV = FPV.flatten()
                                           
            self.vt['voicing'] = (1+np.tanh(3 * FPV[2]))/2
            self.vt['pressure'] = FPV[1]
            self.vt['pressure0'] = self.vt['pressure'] > 0.01
            self.vt['f0'] = 100+20*FPV[0]
            
            af0 = np.maximum(0*af,af)
            k = 0.025
                  
            idx_af0 = np.where(np.logical_and(np.greater(af0,0*af0), np.less(af0,k*np.ones(af0.shape)))) 
            af0[idx_af0] = k*np.ones((len(af0),))
            minaf = np.min(af)
            minaf0 = np.min(af0)
            self.vt['af'] = af
            
        #---------------------------------------------------- %      display(af)
        #---------------------------- %      DIVA_x.af_sample=DIVA_x.af_sample+1
        #------------------------------- %      DIVA_x.af(:,DIVA_x.af_sample)=af
        
        #======================================================================
        # %    tracks place of articulation
        #======================================================================
        
            if minaf0 == 0: 
                release = 0
                self.vt['opening_time'] = 0
                self.vt['closure_time'] = self.vt['closure_time'] + 1
                self.vt['closure_position'] = np.where(af0 == 0)[0][-1]
                if not self.vt['closed']:
                    closure = self.vt['closure_position']
                else:
                    closure = False
                self.vt['closed'] = True
            else:
                if self.vt['closed']:
                    release = self.vt['closure_position']
                    release_closure_time = self.vt['closure_time']
                else: 
                    release = 0
                if self.vt['pressure0'] and not synth.pressure0:
                    self.vt['opening_time'] = 0
                self.vt['opening_time'] = self.vt['opening_time'] + 1
                self.vt['closure_time'] = 0
                self.vt['closure_position'] = np.argmin(af)
                closure = 0
                self.vt['closed'] = False
        
        #%     display(vt.closed)
            if release > 0:
                af = np.maximum(k*np.ones(af.shape),af)
                minaf = np.max((k, minaf))
                minaf0 = np.max((k, minaf0))
            
            if release > 0: 
                self.vt['f0'] = (0.95 + 0. * random.random()) * voices[opt.voices]['F0']  #Maybe we must have more control over randomness
                synth.pressure = 0   #%modulation=0 
            elif (self.vt['pressure0']  and not synth.pressure0): 
                self.vt['f0'] = (0.95 + 0. * random.random()) * voices[opt.voices]['F0']
                synth.pressure = self.vt['pressure']
                synth.f0 = 1.25 * self.vt['f0'] 
                synth.pressure = 1.     #%synth.modulation=1 
            elif (not self.vt['pressure0'] and synth.pressure0 and not self.vt['closed']):
                synth.pressure = synth.pressure/10.
            
            
            #% computes glottal source
            synth.samplesperperiod = int(np.ceil(synth.fs/synth.f0))
            pp = [0.6 , 0.2 - 0.1 * synth.voicing, 0.25] #%10+.15*max(0,min(1,1-vt.opening_time/100))]
            
            glt_in = array(np.arange(0.,1.,1./synth.samplesperperiod))
            synth.glottalsource = 10 * 0.25 * glotlf(0,glt_in,pp) + \
                                  10 * 0.025 * synth.k1 * glotlf(1,glt_in,pp)
            numberofperiods = int(synth.numberofperiods)
                
            #% computes vocal tract filter
            
            synth.filt, synth.f, synth.filt_closure = \
                                 self.a2h(af0, d, synth.samplesperperiod, synth.fs,self.vt['closure_position'], minaf0)
                                 
            synth.filt = synth.filt[:,0] #Changes dimensions from (n,1) to (n,)
            synth.filt_closure = synth.filt_closure[:,0] #Changes dimensions from (n,1) to (n,)
            
            
            synth.filt = 2. * synth.filt / np.maximum(eps,synth.filt[0])
            synth.filt[0] = synth.filt[0] * 0. 
            synth.filt_closure = 2. * synth.filt_closure/np.maximum(eps,synth.filt_closure[0])
            synth.filt_closure[0] = synth.filt_closure[0] * 0 
              
            #% computes sound signal
            w = np.linspace(0., 1., synth.samplesperperiod)
            if release > 0: #,%&&synth.pressure>.01,
                u = synth.voicing * 1.* 0.010 * (synth.pressure + 20 * synth.pressurebuildup)\
                                  * synth.glottalsource + (1-synth.voicing) * 1. * 0.010\
                                  * (synth.pressure + 20. * synth.pressurebuildup)\
                                  * array([0. * random.random() + 1 for i in range(synth.samplesperperiod)])
        #%         if release_closure_time<40
        #%             u=1*.010*synth.pressure*synth.glottalsource;%.*(0.25+.025*randn(synth.samplesperperiod,1)); % vocal tract filter
        #%         else
        #%             u=1*.010*(synth.pressure+synth.pressurebuildup)*randn(synth.samplesperperiod,1)
        #%         end
                v0 = np.real(   np.fft.ifft(np.multiply(np.fft.fft(u), synth.filt_closure))  )
                numberofperiods = numberofperiods - 1
                synth.pressure = synth.pressure/10.
                vnew = v0[range(synth.samplesperperiod)]
                v0 = np.multiply( array([1- x for x in w]) ,\
                                 synth.sample[(np.ceil(len(synth.sample)\
                                                      * array(range(synth.samplesperperiod))\
                                                      / synth.samplesperperiod)).astype(int)]) + \
                                                      np.multiply(w,vnew)
                synth.sample=vnew   
                 
            else:
                v0=array([])
          
            
            if numberofperiods>0:
                # %u=0.25*synth.modulation*synth.pressure*synth.glottalsource.*(1+.1*randn(synth.samplesperperiod,1)) % vocal tract filter
                u = 0.25 * synth.pressure * np.multiply(synth.glottalsource, array([ 0.0 * random.random() + 1 for i in range(synth.samplesperperiod)])) #% vocal tract filter #NOISE
                u = synth.voicing * u + (1 - synth.voicing) * 0.025 * synth.pressure * array([ 0.0 * random.random() + 1 for i in range(synth.samplesperperiod)])
                if minaf0>0 and minaf0<=k:
                    u = minaf / k * u + (1-minaf/k) * 0.02 * synth.pressure * array([ 0.0 * random.random() + 1 for i in range(synth.samplesperperiod)])
                    
                v = np.real(np.fft.ifft(np.multiply(np.fft.fft(u), synth.filt)))
                
                vnew = v[0:synth.samplesperperiod]
                v = np.multiply([1 - x for x in w],\
                                  synth.sample[[int(np.ceil(len(synth.sample) * (x)/synth.samplesperperiod)) for x in range(synth.samplesperperiod)]])\
                                  + np.multiply(w, vnew)
                synth.sample = vnew
                
                if numberofperiods>1:
                    v= np.concatenate((v,np.repeat(vnew,numberofperiods-1)))
                
            else:
                v = array([])
                
            print(time, ': ', len(v0), ', ,', len(v))
                
            v = np.concatenate((v0,v))
            v = v + [.0000 * random.random() for i in range(len(v))]    #Does not support matices
            v = np.divide([1-np.exp(-x) for x in v]   , [1+np.exp(-x) for x in v] )
            try:
                s[[synth.samplesoutput + x for x in range(len(v))]] = v
            except IndexError:
                s = np.concatenate((s, np.zeros(synth.samplesoutput + len(v)-len(s),)))
                s[[synth.samplesoutput + x for x in range(len(v))]] = v
                
            time = time + len(v)/synth.fs
            synth.samplesoutput = synth.samplesoutput + len(v)
            
            # % computes f0/amp/voicing/pressurebuildup modulation
            synth.pressure0 = self.vt['pressure0']
            alpha = np.min((1,(0.1)*synth.numberofperiods))
            beta= 100. / synth.numberofperiods
            synth.pressure = synth.pressure + alpha*(self.vt['pressure']*(np.max((1, 1.5 - self.vt['opening_time']/beta)))-synth.pressure)
            alpha = np.min((1,.5*synth.numberofperiods))
            beta = 100./synth.numberofperiods
            synth.f0 = synth.f0 + 2 * np.square(alpha) * 0. * random.random() + alpha * (self.vt['f0'] * np.max((1. ,1.25-self.vt['opening_time']/beta))-synth.f0) #%147;%120;
            synth.voicing = np.max((0, np.min((1, synth.voicing + 0.5 * (self.vt['voicing']-synth.voicing) ))))
            #%synth.modulation=max(0,min(1, synth.modulation+.1*(2*(vt.pressure>0&&minaf>-k)-1) ))
            alpha = np.min((1 , 0.1 * synth.numberofperiods))
            synth.pressurebuildup = np.max((0, np.min((1, synth.pressurebuildup + alpha * (2 * float(self.vt['pressure'] > 0 and minaf < 0 ) -1) ))))
            synth.numberofperiods = np.max((1,numberofperiods))

        s = s[0:int(np.ceil(synth.fs*ndata*dt))]
        return s, af
    
        
    def get_sample(self, Art):
        '''
            Art numpy array 1 x n_articulators(13) 
            % computes auditory/somatosensory representations
            % Art(1:10) vocaltract shape params
            % Art(11:13) F0/P/V params
            % Aud(1:4) F0-F3 pitch&formants
            % Som(1:6) place of articulation (~ from pharyngeal to labial closure)
            % Som(7:8) P/V params (pressure,voicing)
        '''
        
        #======================================================================
        # computes vocal tract configuration
        #======================================================================
        idx = range(10)
                    
        x = np.multiply(self.vt['vt_scale'][idx].flatten(), Art[idx])
        Outline = self.vt['vt_average'] + np.dot(self.vt['vt_base'][:,idx].reshape((396,10), order = 'F'), x.reshape((10,1))) #Keep in mind reshaping order
        Outline = Outline.flatten()
        
        #=======================================================================
        # % computes somatosensory output (explicitly from vocal tract configuration)        
        #=======================================================================
        
        Som = np.zeros((8,))
        a, b, sc, af, d = self.xy2ab(Outline)
        Som[0:6] = np.maximum(-1*np.ones((len(sc),)),np.minimum(np.ones((len(sc),)), -tanh(1*sc) ))
        Som[6:] = Art[-2:]
        
        Aud = np.zeros((4,))
        Aud[0] = 100 + 50 * Art[-3]
        dx = Art[idx] - self.fmfit['fmfit_mu']
        p = -1 * np.sum(np.multiply(np.dot(dx, self.fmfit['fmfit_isigma']), dx),axis=1)/2
        p = np.multiply(self.fmfit['fmfit_p'].flatten(),np.exp(p-(np.max(p)*np.ones(p.shape))))
        p = p / np.sum(p)
        px = np.dot(p.reshape((len(p),1)) , np.append(Art[idx],1).reshape((1,len(idx)+1)))
        Aud[1:4] = np.dot(self.fmfit['fmfit_beta_fmt'], px.flatten(1))
        
        ####################### Not implemente yet (nargout = 2 or 3)
        #---------------------------- if ~isempty(fmfit)&&nargout>1&&nargout<=3,
            #------------------------------------ Som(1:6)=fmfit.beta_som*px(:)
            #------------------------------------------ Som(7:8)=Art(end-1:end)
        #------------------------------------------------------------------- end
        
        return Aud, Som, Outline, af, d 

    def xy2ab(self, x, y = False):  #x -> Outline (column matrix)
        if not hasattr(self, 'ab_alpha'):
            amax = 220
            alpha = array([1, 1, 1, 1, 1, 1, 1])
            beta = array([.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25])
            idx = [range(60), range(60,70), range(70,80), range(80,120), range(120,150), range(150,190), range(190,amax)]
            ab_alpha = np.zeros((amax,))
            ab_beta = np.zeros((amax,))
            for n1 in range(len(idx)):
                ab_alpha[idx[n1]] = alpha[n1]
                ab_beta[idx[n1]] = beta[n1]
            
            h = self.hanning.flatten()
            #h = hanning(51)/np.sum(hanning(51)) #Not same result as in hanning (matlab hann vs hanning) Source of numnerial differences
            idx_2 = np.zeros((25,))
            idx_2 = np.concatenate((idx_2, np.array(range(amax))))
            idx_2 = np.concatenate((idx_2, (amax-1)*np.ones((25,))))
            idx_2 = idx_2.tolist()
            
            ab_alpha = np.convolve(ab_alpha[idx_2],h,'valid')
            ab_beta = np.convolve(ab_beta[idx_2],h,'valid')
            self.ab_alpha = ab_alpha
            self.ab_beta = ab_beta
            
        if not y:
            i = np.array([0 + 1j])[0]
            x=np.exp(-i*np.pi/12) * x
            y=np.imag(x)
            x=np.real(x)
        
        ab_alpha = self.ab_alpha
        ab_beta = self.ab_beta
        #=======================================================================
        # % Grid
        #=======================================================================
        x0 = 45       #%90
        y0 = -100     #%-60
        r = 60        #%30
        k = np.pi * (r/2)
        d = 0.75/10    #%unitstocm
        
        a = np.zeros(x.shape)
        b = np.zeros(x.shape)
        i1 = np.where(np.less(y,y0))
        i2 = np.where(np.less(x,x0))
        i3 = np.where(np.logical_and(np.greater_equal(y, y0), np.greater_equal(x, x0)))
        
        #=======================================================================
        # % a,b: "linearized" coordinates along vocal tract
        #======================================================================= 
        
        a[i1] = y[i1]-y0
        b[i1] = x[i1]-x0
        a[i2] = k+x0-x[i2]
        b[i2] = y[i2]-y0
        z = x[i3]-x0 + i*(y[i3]-y0)
        a[i3] = r*np.angle(z)
        b[i3] = np.abs(z)
        #=======================================================================
        # % tube area
        #=======================================================================
         
        olips = range(29,45) 
        ilips = range(256,302)
        #------------------------------------------------- owall = range(44,164)  #Not used in Matlab
        iwall = range(163+10,257)
        oall = range(29,164)
        iall = range(163,302)
        xmin = -20 
        ymin = -160
        amin = ymin-y0
        amax=np.ceil((x0-xmin+k-amin))  #here is the variability of the af vector
        
        fact = 3
        
        # wallab1 = accumarray(max(1,min(fact*9, ceil(fact*9*(a(oall)-amin)/amax))),b(oall),[fact*9,1],@min,nan)
        # wallab2 = accumarray(max(1,min(fact*9, ceil(fact*9*(a(iwall)-amin)/amax))),b(iwall),[fact*9,1],@max,nan)
        
        idx_wallab1 = np.int_(np.maximum(np.zeros((len(oall),)),np.minimum(((fact*9)-1)*np.ones((len(oall),)), np.ceil(fact*9*(a[oall]-amin)/amax)-1)))
        idx_wallab2 = np.int_(np.maximum(np.zeros((len(iwall),)),np.minimum(((fact*9)-1)*np.ones((len(iwall),)), np.ceil(fact*9*(a[iwall]-amin)/amax)-1)))
        wallab1 = aggregate(idx_wallab1,b[oall], size = fact*9, func = 'min', fill_value=None)
        wallab2 = aggregate(idx_wallab2,b[iwall], size = fact*9, func = 'max', fill_value=None)
        
        lipsab1 = np.nanmin(b[olips])
        lipsab2 = np.nanmax(b[ilips])
        
        mind_precursor = wallab1[range(fact*2,fact*8)]-wallab2[range(fact*2,fact*8)]
        mind = np.nanmin(mind_precursor.reshape((fact,6), order = 'F'), axis = 0)
        sc = mind[range(4)]    
        sc = np.append(sc,np.nanmin(mind[range(4,6)]))
        sc = np.append(sc,lipsab1-lipsab2) 
        sc = d * sc #In Matlab this is a column vector
        
        w = 2
        
        #ab1 = aggregate(max(1,min(amax, round((a(oall)-amin)))),b(oall),[amax,1],min,nan)
        #ab2 = aggregate(max(1,min(amax, round((a(iall)-amin)))),b(iall),[amax,1],max,nan)
        
        idx_ab1 = np.int_(np.maximum(np.zeros((len(oall),)),np.minimum((amax-1)*np.ones((len(oall),)), np.round(a[oall]-amin-1))))
        idx_ab2 = np.int_(np.maximum(np.zeros((len(iall),)),np.minimum((amax-1)*np.ones((len(iall),)), np.round(a[iall]-amin-1))))
        
        ab1 = aggregate(idx_ab1,b[oall],size = amax, func= 'min', fill_value=None)
        ab2 = aggregate(idx_ab2,b[iall],size = amax, func = 'max', fill_value=None)
        
        ab1[np.isnan(ab1)] = np.inf
        ab2[np.isnan(ab2)] = -np.inf
        for n1 in range(w):
            ab1[1:-1] = np.minimum(np.minimum(ab1[1:-1], ab1[0:-2]), ab1[2:])
        for n1 in range(w):
            ab2[1:-1] = np.maximum(np.maximum(ab2[1:-1], ab2[0:-2]), ab2[2:])
        
        ab1[np.isinf(ab1)] = np.NINF
        idx_af = np.where(np.logical_and(  np.greater(ab1,0)   ,np.greater(ab2,0)))[0]
        af = d*(ab1[idx_af]-ab2[idx_af])
        idx = None
        for ii in range(len(af)):
            if af[ii] > 0: 
                idx =  ii
                break
       
        #=======================================================================
        # % af: area function
        #=======================================================================
        af_tmp = np.minimum(np.zeros((len(af),)),af)
        af_power = np.power(np.maximum(np.zeros((len(af),)),af), ab_beta[idx_af])
        af = af_tmp + np.multiply(ab_alpha[idx_af], af_power ) # np.power(np.maximum(np.zeros((len(af),)),af), ab_beta[idx_af]) introduces an small error
        
        
        af = af[idx:]
        for ii in range(len(af)-1,-1,-1):
            if np.isinf(af[ii]):
                af = np.delete(af, -1) 
            else:
                break
        return a, b, sc, af, d
    
    def a2h(self, a,l,n,fs = None, closure = None, mina = None):   
        if not isinstance(mina, float):
            mina = np.amin(a,axis=0)
        if not isinstance(closure, float) and not isinstance(closure, int):
            closure = 0
        if not isinstance(fs, float) and not isinstance(fs, int):
            fs = 11025
        #if np.sum(len(np.where(a.shape>1))) <= 1:
        if not isinstance(a, np.ndarray):
            #a = a.flatten('F')
            N, M = (1,1)
        else:
            if len(a.shape)==1:
                N = a.shape[0]
                M = 1
            else:
                N, M = a.shape
            
            
        #if not np.sum(len(np.where(l.shape>1))) <= 1:
        if not isinstance(l, np.ndarray):
            #l = l.flatten('F')
            NL, ML =(1, 1)
        else:
            NL, ML = l.shape
            
            
        c=34326 #% speed of sound (cm/s)

        m = int(np.ceil(n/2)+1)
        f = fs*array(range(int(np.ceil(n/2)+1)))/n
        t = l / c
        Rrad = 0.9 * np.exp(-(np.abs(f)/4e3)**2) #% reflection at lips (low pass)
        H = np.zeros((m,M)) + 0j*np.zeros((m,M)) #%H=zeros(n,M);
        Hc = np.zeros((m,M)) + 0j*np.zeros((m,M))
        if mina==0:
            a = array([np.max((0.05, x)) for x in a])
        #%k=0.995;%.999;
        for nM in range(M):
            #if 1: #%mina(nM)>0,
            coswt = np.cos(2 * np.pi * f * sub_array(t, idx_y = np.min((ML,nM))-1))
            sinwt = 1j * np.sin(2 * np.pi * f * sub_array(t, idx_y = np.min((ML,nM))-1))
            #R=[.9;(a(2:N,nM)-a(1:N-1,nM))./max(eps,a(2:N,nM)+a(1:N-1,nM))];
            R = np.concatenate(([0.9], np.divide(sub_array(a,idx_x = range(1,N,1), idx_y = nM-1) \
                                                                         - sub_array(a,idx_x=range(N-1), idx_y = nM-1), \
                                                    np.maximum(eps*np.ones((N-1,)),sub_array(a,idx_x = range(1,N,1), idx_y = nM-1) \
                                                                         + sub_array(a,idx_x=range(N-1), idx_y = nM-1)))))
            
            U = coswt+sinwt
            V = coswt-sinwt
            #--------------------------------------------------------- if 1:
            h1 = np.ones((m,))                  # % signal at glottis   #In this part of the code there are diferences w.r.t. MATLAB,
            h2 = np.zeros((m,))
            for nN in range(N-1):
                RnN = -R[nN]
                u = h1 + RnN * h2
                v = h2 + RnN * h1   # % reflection
                if closure==nN:
                    Hc[:,nM] = u - v
                if NL == 1:
                    h1 = np.multiply(U, u)
                    h2 = np.multiply(V, v)      #  % delay
                else:   # This case might not be working properly, It hasnot been tested
                    h1 = np.multiply(U[:,nN],u)
                    h2 = np.multiply(V[:,nN],v)
                #print(sum(h1))
                #print(sum(h2))
                #print(RnN)   
                #%h(:,1)=h(:,1)/k;h(:,2)=h(:,2)*k;
            
            u = h1 - np.multiply(Rrad,h2) # %v=h2-Rrad.*h1; #   % reflection
            h = u;              #%h(:,2)=v;
            if closure >= N:  #Indexing here must be debugged
                Hc[:,nM] = u - np.multiply(h2 - np.multiply(Rrad,h1))
            
            ###### CODE IN MATLAB THAT IS NEVER EXECUTED #########3
            
            H[:,nM] = np.divide(np.multiply((np.ones((Rrad.shape)) + Rrad), np.prod(np.ones((R.shape))+R)), sub_array(u, idx_y = 0))
            #Small difference in decimal w.r.t MATLAB
            if closure>0:
                Hc[:,nM] = np.divide(np.multiply((np.ones((Rrad.shape))+Rrad),\
                                                 np.prod(np.ones((N-closure-1,))+R[closure+1:N])*Hc[:,nM]),\
                                                            sub_array(h,idx_y=0))
        
        idxH_x =  [x + 1 for x in range(n-m-1,-1,-1)]
              
        H = np.concatenate( (H , np.conjugate(H[idxH_x,:]) ) ,axis = 0)
        
        Hc = np.concatenate( ( Hc, np.conjugate(Hc[idxH_x,:])),axis = 0)
        
        f = np.concatenate((f,-f[idxH_x]))
        if mina == 0:
            H = 0 * H
        
        return  H, f, Hc
        
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
        x = x[idx_x,:]
    if not idx_y == np.inf:
        x = x[:,idx_y]
    
    return x
                
                
def glotlf(d, t = None, p = None):
    if t is None:
        tt = array(range(99))/100
    else:
        tt = t-np.floor(t)
    
    u = np.zeros((len(tt),))
    de = array([0.6, 0.1, 0.2])
    
    if p is None: 
        p = de 
    else:
        p = np.concatenate((p, de[len(p):1]))
    
    te = p[0]
    mtc = te-1
    e0 = 1
    wa = np.pi / (te * ( 1 - p[2] ))
    a = -np.log(-p[1]* np.sin(wa * te)) / te
    inta = e0 * ((wa / np.tan(wa*te)-a)/p[1]+wa) / (a**2.0 + wa**2.0)
    
    #----------------------------------------- % if inta<0 we should reduce p(2)
    #------------------------- % if inta>0.5*p(2)*(1-te) we should increase p(2)
    
    rb0 = p[1] * inta
    rb = rb0
    
    #--------------------------- % Use Newton to determine closure time constant
    #----------------------------------- % so that flow starts and ends at zero.
    
    for i in range(4):
        kk = 1 - np.exp(mtc/rb)
        err = rb + mtc * (1 /kk - 1 ) - rb0
        derr = 1 - (1 - kk) * (mtc/ rb / kk)**2.
        rb = rb - err/derr
    e1 = 1 /(p[1]*(1 - np.exp( mtc / rb ) ))
    
    pre_ta_tb = np.less(tt, te)
    ta = np.where(pre_ta_tb)
    tb = np.where(np.logical_not(pre_ta_tb))
    
    if d == 0:
        u[ta] = e0 * (np.multiply(np.exp(a*tt[ta]),(a*np.sin(wa*tt[ta])-wa*np.cos(wa*tt[ta])))+wa)/(a**2.+wa**2.)
        u[tb] = e1 * (np.exp(mtc/rb)*(tt[tb]-1-rb)+np.exp((te-tt[tb])/rb)*rb)
    elif d==1:
        u[ta] = e0 * np.multiply(np.exp(a*tt[ta]),np.sin(wa*tt[ta]))
        u[tb] = e1 * (np.exp(mtc/rb)-np.exp((te-tt[tb])/rb))
    elif d==2:
        u[ta] = e0 * np.multiply(np.exp(a*tt[ta]),(a*np.sin(wa*tt[ta])+wa*np.cos(wa*tt[ta])))
        u[tb] = e1 * np.exp((te-tt[tb])/rb)/rb
    else:
        print('Derivative must be 0,1 or 2')
        raise ValueError
    
    return u
     
    doc = '''
        %GLOTLF   Liljencrants-Fant glottal model U=(D,T,P)
        % d is derivative of flow waveform: must be 0, 1 or 2
        % t is in fractions of a cycle
        % p has one row per output point
        %    p(:,1)=open phase [0.6]
        %    p(:,2)=+ve/-ve slope ratio [0.1]
        %    p(:,3)=closure time constant/closed phase [0.2]
        % Note: this signal has not been low-pass filtered
        % and will therefore be aliased
        %
        % Usage example:    ncyc=5;
        %            period=80;
        %            t=0:1/period:ncyc;
        %            ug=glotlf(0,t);
        %            plot(t,ug)
        
        
        %      Copyright (C) Mike Brookes 1998
        %
        %      Last modified Thu Apr 30 17:22:00 1998
        %
        %   VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   This program is free software; you can redistribute it and/or modify
        %   it under the terms of the GNU General Public License as published by
        %   the Free Software Foundation; either version 2 of the License, or
        %   (at your option) any later version.
        %
        %   This program is distributed in the hope that it will be useful,
        %   but WITHOUT ANY WARRANTY; without even the implied warranty of
        %   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        %   GNU General Public License for more details.
        %
        %   You can obtain a copy of the GNU General Public License from
        %   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
        %   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        '''
            
        