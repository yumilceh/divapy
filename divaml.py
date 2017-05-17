'''
Created on Feb 5, 2016
This sensorimor system defines the DIVA agent used for the CCIA 2015's paper
@author: Juan Manuel Acevedo Valle
'''

# import sys
# import wave
import subprocess as sp
import math
import numpy as np
import pymatlab as ml
import matplotlib.pyplot as plt
from matplotlib import animation
from scipy.integrate import odeint
from scipy import linspace
from scipy.io.wavfile import write
import os

# from matplotlib.pyplot import autoscale
# from matplotlib.animation import Animation
# from scipy.interpolate.interpolate_wrapper import block

diva_output_scale = [100.0, 500.0, 1500.0, 3000.0]


class Diva(object):
    def __init__(self):
        n_motor = 13
        n_sensor = 4
        output_scale = diva_output_scale
        min_motor_values = np.array([-3, -3, -3, -3, -3, -3, -3, -3, -3, -3, 0., 0., 0.])
        max_motor_values = np.array([3.3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1., 1., 1.])

        name = 'DivaWrapper'
        self.name = name

        self.n_motor = n_motor
        self.n_sensor = n_sensor

        self.min_motor_values = min_motor_values
        self.max_motor_values = max_motor_values

        self.motor_command = [0.0] * n_motor
        self.sensor_out = [0.0] * n_sensor

        self.matlabSession = ml.session_factory()
        abs_path = os.path.dirname(os.path.abspath(__file__))

        command_ = 'cd ' + abs_path + '/DivaMatlab/'

        self.matlabSession.run(command_)  # Path to DivaMatlab functions
        self.matlabSession.putvalue('outputScale', output_scale)

    def get_audsom(self, art, scale=True):
        self.matlabSession.putvalue('art', art)
        self.matlabSession.run('[aud, Som, outline, af] = diva_synth(art\',\'audsom\')')
        aud =  np.transpose(self.matlabSession.getvalue('aud'))
        som = np.transpose(self.matlabSession.getvalue('Som'))
        outline = np.transpose(self.matlabSession.getvalue('outline'))
        af_ = self.matlabSession.getvalue('af')
        af = []
        for i in range(art.shape[0]):
            af += [af_[:,i]]

        return aud, som, outline, af

    def get_sound(self, arts, play=0, save=0, file_name='vt'):  # based on explauto
        ts = 0.005
        # arts ->   n_samples * 13

        self.matlabSession.putvalue('artStates', arts)

        self.matlabSession.run('sound_wave = diva_synth(artStates\', \'sound\')')
        if (play):
            self.play_sound()
        if (save):
            scaled = np.int16(self.sound_wave / np.max(np.abs(self.sound_wave)) * 32767)
            write(file_name + '.wav', 11025, scaled)

        sound = self.matlabSession.getvalue('sound_wave')
        return sound


    def get_static_sound(self, art, ts = 0.005, time=0.4):
        n_samples = time/ts
        arts  = np.repeat([art], n_samples, axis=1)
        #print(arts)
        return self.get_sound(arts)

    def plot_sound(self, ax=None):
        if type(ax) is type(None):
            fig, ax = plt.plot(np.float128(xrange(0, len(self.sound_wave))) * self.ts, self.sound_wave)
            return fig, ax
        else:
            ax.plot(np.float128(xrange(0, len(self.sound_wave))) * self.ts, self.sound_wave)

    def play_sound(self):  # keep in mind that DivaMatlab works with ts=0.005
        import pyaudio
        self.pa = pyaudio.PyAudio()  # If pa and stream are not elements of the self object then sound does not play
        self.stream = self.pa.open(format=pyaudio.paFloat32,
                                   channels=1,
                                   rate=11025,
                                   output=True)
        self.stream.start_stream()
        self.stream.write(self.sound_wave.astype(np.float32).tostring())

    def releaseAudioDevice(self):  # any sound in the buffer will be removed
        try:
            self.stream.close()
            self.pa.terminate()
        except:
            pass

    def stop(self):
        del self.matlabSession
