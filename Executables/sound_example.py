'''
Created on Feb 4, 2017

@author: Juan Manuel Acevedo Valle

Exampleon how to use Diva's implementation in python
'''

import h5py, time, pyaudio
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from divapy import Diva


def h5_to_ndarray(key_, file_name):
    with h5py.File(file_name, 'r') as hf:
        data = hf[key_][:]
    return data

if __name__ == '__main__':
    diva_synth = Diva()

    # arts = [-2.82586455, 0.6511439, 0.41705491, -0.16636287, 1.52612695, -0.65089551,
    #         2.71660497, -2.8833028, -1.46663703, -2.10624563, 0.21783021, 0.74918289,
    #         0.48187699]
    # arts = np.tile(arts, (60,1))

    arts = h5_to_ndarray('arts','../data/hard_arts.h5')

    aud, af = diva_synth.get_sound(arts)
    sio.savemat('../DivaMatlab/soundtest.mat', {'Arts2': arts, 'Aud2':aud})
     
    pa = pyaudio.PyAudio()
    stream = pa.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=11025,
                     output=True)
    stream.start_stream()
    stream.write(aud.astype(np.float32).tostring())

    plt.plot(aud)
    plt.show()

    time.sleep(1)
    stream.close()
    pa.terminate()