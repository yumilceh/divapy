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
mldiva_synth = MlDiva()

n_static_test = 100
n_static_sound_test = 100
n_random_sound_test = 100

static_test = True
static_sound_test = True
random_sound_test = True

if __name__ == '__main__':
    static_commands = get_random_motor_set(mldiva_synth,
                                           n_static_test)

    print("Validating with static samples...")
    for i in range(n_static_test):
        try:
            Aud_ml, Som_ml, Outline_ml, af_ml = mldiva_synth.get_audsom(static_commands[i, :])
            Aud, Som, Outline, af = pydiva_synth.get_audsom(static_commands[i, :])
        except ValueError:
            print(static_commands[i, :])
            sys.exit()
        err_aud = np.linalg.norm(Aud_ml - Aud)
        err_som = np.linalg.norm(Som_ml - Som)
        err_ol = np.linalg.norm(np.nan_to_num(Outline_ml - Outline))
        err_af = np.linalg.norm(af_ml - af)

        if err_aud > 1e-3 or err_som > 1e-3 or err_ol > 1e-3 or err_af > 1e-3:
            print("Validation failed in static test.")
            print(err_aud)
            print(err_som)
            print(err_ol)
            print(err_af)
            print(static_commands[i, :])
            static_test = False
    if static_test:
        print("Validation in static test OK.")

    static_commands = get_random_motor_set(mldiva_synth,
                                           n_static_sound_test)
    print("Validating with static sound samples...")
    for i in range(n_static_sound_test):
        m_tmp = static_commands[i, :]
        arts = np.tile(m_tmp, (80, 1))

        sound_wave_ml = mldiva_synth.get_sound(arts)

        sound_wave = pydiva_synth.get_sound(arts)
        # try:
        #
        # except ValueError:
        #     print(m_tmp)
        #     sys.exit()
        err_sw = np.linalg.norm(sound_wave_ml - sound_wave)

        if err_sw > 1e-3:
            print("Validation failed in static sound test.")
            print(err_sw)
            plt.figure()
            plt.plot(sound_wave_ml, 'b')
            plt.hold(True)
            plt.plot(sound_wave, 'r')
            plt.legend(['Matlab', 'Python'])
            plt.show()
            static_sound_test = False

    if static_sound_test:
        print("Validation in static sound test OK.")

    print("Validating with random sound samples...")
    for i in range(n_random_sound_test):
        arts = get_random_motor_set(mldiva_synth, 80)

        sound_wave_ml = mldiva_synth.get_sound(arts)
        try:
            sound_wave = pydiva_synth.get_sound(arts)
        except ValueError:
            print(arts)
            sys.exit()
        err_sw = np.linalg.norm(sound_wave_ml - sound_wave)

        if err_sw > 1e-3:
            print("Validation failed in random sound test.")
            print(err_sw)
            print(arts)
            plt.figure()
            plt.plot(sound_wave_ml, 'b')
            plt.hold(True)
            plt.plot(sound_wave, 'r')
            plt.legend(['Matlab', 'Python'])
            plt.show(block=True)

            static_sound_test = False
            x  = raw_input('Press [key]+[Enter] to cancel or [Enter] to continue...')

            if not len(x)==0:
                break

    if static_sound_test:
        print("Validation in random sound test OK.")
