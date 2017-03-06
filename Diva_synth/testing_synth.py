'''
Created on Feb 4, 2017

@author: Juan Manuel Acevedo Valle
'''

import numpy as np
import matplotlib.pyplot as plt
import time
import os, sys
import scipy.io as sio


if __name__ == '__main__':
    sys.path.append("../../")
    
    from SensorimotorSystems.Diva_Synth import Diva
    diva_synth = Diva()
    #------------------------------------------------------------ art = [0.1]*13
#------------------------------------------------------------------------------ 
    #------------------------------------------------------- art = np.array(art)
#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 
    #------------------------ Aud, Som, Outline, af = diva_synth.get_audsom(art)

    
    arts = [[3.]*13, [0.7]*13]
    arts[0][11:] = ([1]*3)
    arts[1][11:] = ([1]*3)

    arts = np.array(arts).transpose()
    arts = np.tile(arts, 40).transpose()
    Aud, af = diva_synth.get_sound(arts)
    sio.savemat('soundtest.mat', {'Aud2':Aud})
     
    import pyaudio 
    pa = pyaudio.PyAudio() #If pa and stream are not elements of the self object then sound does not play
    stream = pa.open(format=pyaudio.paFloat32,
                     channels=1,
                     rate=11025,
                     output=True)
    stream.start_stream()
    
    stream.write(Aud.astype(np.float32).tostring())
    
    #stream.close()
    #pa.terminate()
    #--------------------------------------------------------------- x = input()
    
    
    plt.plot(af)
    plt.show()
    #------------------------------------------------ arts = [[.2]*13, [0.2]*13]
    #---------------------------------------------------- arts[0][11:] = ([1]*3)
    #---------------------------------------------------- arts[1][11:] = ([1]*3)
#------------------------------------------------------------------------------ 
    #----------------------------------------------------- arts = np.array(arts)
#------------------------------------------------------------------------------ 
    #-------------------------------------- Aud, af = diva_synth.get_sound(arts)