'''
Created on April 2017, 2017

@author: Juan Manuel Acevedo Valle

Script used to compare the performances beetween Diva's python implementation
and original implementation Matlab
'''
import numpy as np
import sys
from matplotlib import pyplot as plt
from divapy import Diva as PyDiva
from divaml import Diva as MlDiva
import time

def get_random_motor_set(system, n_samples):
    n_motor = system.n_motor
    raw_rnd_data = np.random.random((n_samples, n_motor))
    min_values = system.min_motor_values
    max_values = system.max_motor_values

    min_values = np.array(n_samples * [np.array(min_values)])
    max_values = np.array(n_samples * [np.array(max_values)])

    motor_commands = min_values + raw_rnd_data * (max_values - min_values)

    return motor_commands

pydiva_synth1 = PyDiva()
pydiva_synth2 = PyDiva()
mldiva_synth = MlDiva()

n_static_sound_test = 100

static_sound_test = True

if __name__ == '__main__':

    static_commands = get_random_motor_set(mldiva_synth,
                                           n_static_sound_test)
    print("Evaluating performance with with static sound samples...")

    t = time.time()
    for i in range(n_static_sound_test):
        m_tmp = static_commands[i, :]
        arts = np.tile(m_tmp, (80, 1))
        sound_wave2= pydiva_synth2.get_sound(arts)

    elapsed_ml = time.time() - t

    t = time.time()
    for i in range(n_static_sound_test):
        try:
            m_tmp = static_commands[i, :]
            arts = np.tile(m_tmp, (80, 1))
            sound_wave1= pydiva_synth1.get_sound_opt(arts)
        except ValueError:
            print(m_tmp)
            sys.exit()
    elapsed_py = time.time()-t

    print("{} static sound samples. Matlab running time: {}. Python running time {}".
          format(n_static_sound_test, elapsed_ml, elapsed_py))

