'''
Created on Feb 4, 2017

@author: Juan Manuel Acevedo Valle

Exampleon how to use Diva's implementation in python
'''

import h5py, time
import numpy as np
import matplotlib.pyplot as plt
from divapy import Diva as Divapy
from divaml import Diva as Divaml

def h5_to_ndarray(key_, file_name):
    with h5py.File(file_name, 'r') as hf:
        data = hf[key_][:]
    return data

if __name__ == '__main__':
    diva_synth = Divapy()
    diva_synth_ml = Divaml()

    art = [-2.82586455, 0.6511439, 0.41705491, -0.16636287, 1.52612695, -0.65089551,
            2.71660497, -2.8833028, -1.46663703, -2.10624563, 0.21783021, 0.74918289,
            0.48187699]

    art = h5_to_ndarray('arts','../data/hard_arts.h5')

    # art = [0.1]*10 + [1]*3

    # Using get_sample
    aud, som, outline, af, d = diva_synth.get_sample(np.array(art[0,:]))
    print(aud)
    print(som)
    print(af)

    # Using get_audsom
    arts = np.tile(art, (100,1))
    aud, som, outline, af = diva_synth.get_audsom(arts.transpose())
    fig1, ax1 = plt.subplots(2,2)
    plt.sca(ax1[0, 0])
    plt.plot(aud)
    plt.sca(ax1[0, 1])
    plt.plot(som)

    # Using plot_outline
    diva_synth.plot_outline(art, axes=ax1[1, 0])

    # Using get_sound  and play_sound
    sound = diva_synth.get_sound(art)
    sound2 = diva_synth_ml.get_sound(art)
    plt.sca(ax1[1, 1])
    plt.plot(sound)
    plt.hold(True)
    plt.plot(sound2)
    diva_synth.play_sound(sound)

    plt.show()
    time.sleep(1)
