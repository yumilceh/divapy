'''
Created on Feb 20, 2017

@author: Juan Manuel Acevedo Valle

Script used to compare the Diva's python implementation against the original implementation
in Matlab
'''
import numpy as np
import sys
from matplotlib import pyplot as plt
from divapy import Diva as PyDiva
from divaml import Diva as MlDiva

def get_random_motor_set(system, n_samples):
    n_motor = system.n_motor
    raw_rnd_data = np.random.random((n_samples, n_motor))
    min_values = system.min_motor_values
    max_values = system.max_motor_values

    min_values = np.array(n_samples * [np.array(min_values)])
    max_values = np.array(n_samples * [np.array(max_values)])

    motor_commands = min_values + raw_rnd_data * (max_values - min_values)

    return motor_commands

pydiva_synth = PyDiva()
mldiva_synth = PyDiva()

foo = MlDiva()

n_static_test = 1
n_static_sound_test = 1
n_random_sound_test = 200

static_test = False
static_sound_test = True
random_sound_test = True



if __name__ == '__main__':
    static_commands = get_random_motor_set(foo,
                                           n_static_sound_test)

    print("Validating with static sound samples...")
    for i in range(n_static_sound_test):
        m_tmp = static_commands[i, :]
        arts = np.tile(m_tmp, (80, 1))

        sound_wave1 = mldiva_synth.get_sound(arts)
        sound_wave2 = pydiva_synth.get_sound_opt(arts)
        # try:
        #
        # except ValueError:
        #     print(m_tmp)
        #     sys.exit()
        err_sw = np.linalg.norm(sound_wave1 - sound_wave2)

        if err_sw > 1e-3:
            print("Validation failed in static sound test.")
            print(err_sw)
            plt.figure()
            plt.plot(sound_wave1, 'b')
            plt.hold(True)
            plt.plot(sound_wave2, 'r')
            plt.legend(['Matlab', 'Python'])
            plt.show()
            static_sound_test = False

    if static_sound_test:
        print("Validation in static sound test OK.")

    print("Validating with random sound samples...")
    for i in range(n_random_sound_test):
        arts = get_random_motor_set(foo, 80)

        try:
            sound_wave1 = pydiva_synth.get_sound(arts)
            sound_wave2 = mldiva_synth.get_sound_opt(arts)
        except ValueError:
            print(arts)
            sys.exit()
        err_sw = np.linalg.norm(sound_wave2 - sound_wave1)

        if err_sw > 1e-3:
            print("Validation failed in random sound test.")
            print(err_sw)
            print(arts)
            plt.figure()
            plt.plot(sound_wave2, 'b')
            plt.hold(True)
            plt.plot(sound_wave1, 'r')
            plt.legend(['Matlab', 'Python'])
            plt.show(block=True)

            static_sound_test = False
            x  = raw_input('Press [key]+[Enter] to cancel or [Enter] to continue...')

            if not len(x)==0:
                break

    if static_sound_test:
        print("Validation in random sound test OK.")
