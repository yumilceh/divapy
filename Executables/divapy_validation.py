'''
Created on Feb 20, 2017

@author: Juan Manuel Acevedo Valle

Script used to compare the Diva's python implementation agains the original implementation
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

n_static_test = 1
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
        mldiva_synth.set_action(static_commands[i, :])
        mldiva_synth.vocalize()

        try:
            Aud, Som, Outline, af = pydiva_synth.get_audsom(static_commands[i, :])
        except ValueError:
            print(static_commands[i, :])
            sys.exit()
        err_aud = np.linalg.norm(mldiva_synth.aud - Aud)
        err_som = np.linalg.norm(mldiva_synth.som - Som)
        err_ol = np.linalg.norm(np.nan_to_num(mldiva_synth.outline - Outline))
        err_af = np.linalg.norm(mldiva_synth.af - af)

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
        mldiva_synth.get_sound(arts)

        sound_wave, af = pydiva_synth.get_sound(arts)
        # try:
        #
        # except ValueError:
        #     print(m_tmp)
        #     sys.exit()
        err_sw = np.linalg.norm(mldiva_synth.sound_wave - sound_wave)

        if err_sw > 1e-3:
            print("Validation failed in static sound test.")
            print(err_sw)
            plt.figure()
            plt.plot(mldiva_synth.sound_wave, 'b')
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

        mldiva_synth.get_sound(arts)
        try:
            sound_wave, af = pydiva_synth.get_sound(arts)
        except ValueError:
            print(arts)
            sys.exit()
        err_sw = np.linalg.norm(mldiva_synth.sound_wave - sound_wave)

        if err_sw > 1e-3:
            print("Validation failed in random sound test.")
            print(err_sw)
            print(arts)
            plt.figure()
            plt.plot(mldiva_synth.sound_wave, 'b')
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
