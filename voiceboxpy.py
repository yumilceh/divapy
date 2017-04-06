"""
Created on Feb 3, 2017

@author: Juan Manuel Acevedo Valle
@department: Automatic Control Department (Knowledge Engineering Research Group)
@institution: Universitat Polit\`{e}cnica de Catalunya
@email: jmavbpl@gmail.com

Python code for function(s) in:
VOICEBOX home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html

"""

import numpy as np
from numpy import array


def glotlf(d, t=None, p=None):
    if t is None:
        tt = array(range(99)) / 100.
    else:
        tt = t - np.floor(t)

    u = np.zeros((len(tt),))
    de = array([0.6, 0.1, 0.2])

    if p is None:
        p = de
    else:
        p = np.concatenate((p, de[len(p):1]))

    te = p[0]
    mtc = te - 1
    e0 = 1
    wa = np.pi / (te * (1 - p[2]))
    a = -np.log(-p[1] * np.sin(wa * te)) / te
    inta = e0 * ((wa / np.tan(wa * te) - a) / p[1] + wa) / (a ** 2.0 + wa ** 2.0)

    # ----------------------------------------- % if inta<0 we should reduce p(2)
    # ------------------------- % if inta>0.5*p(2)*(1-te) we should increase p(2)

    rb0 = p[1] * inta
    rb = rb0

    # --------------------------- % Use Newton to determine closure time constant
    # ----------------------------------- % so that flow starts and ends at zero.

    for i in range(4):
        kk = 1 - np.exp(mtc / rb)
        err = rb + mtc * (1 / kk - 1) - rb0
        derr = 1 - (1 - kk) * (mtc / rb / kk) ** 2.
        rb = rb - err / derr
    e1 = 1 / (p[1] * (1 - np.exp(mtc / rb)))

    pre_ta_tb = np.less(tt, te)
    ta = np.where(pre_ta_tb)
    tb = np.where(np.logical_not(pre_ta_tb))

    if d == 0:
        u[ta] = e0 * (np.multiply(np.exp(a * tt[ta]), (a * np.sin(wa * tt[ta]) - wa * np.cos(wa * tt[ta]))) + wa) / (
            a ** 2. + wa ** 2.)
        u[tb] = e1 * (np.exp(mtc / rb) * (tt[tb] - 1 - rb) + np.exp((te - tt[tb]) / rb) * rb)
    elif d == 1:
        u[ta] = e0 * np.multiply(np.exp(a * tt[ta]), np.sin(wa * tt[ta]))
        u[tb] = e1 * (np.exp(mtc / rb) - np.exp((te - tt[tb]) / rb))
    elif d == 2:
        u[ta] = e0 * np.multiply(np.exp(a * tt[ta]), (a * np.sin(wa * tt[ta]) + wa * np.cos(wa * tt[ta])))
        u[tb] = e1 * np.exp((te - tt[tb]) / rb) / rb
    else:
        # print('Derivative must be 0,1 or 2')
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
