# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/03/2018 12:45 PM
"""

from __future__ import division
from core import env
import math
import numpy as np

# I am pulling in infromation from the hydraulics VBA editors in stream depetionv3.xls


def Q_4(t, K, epsilon, l):
    """
    This calculates the total flow depletion lost from a stream when a well abstracts a flow
    Q from a delayed-yield (semi-confined) aquifer. All input and output variables are
    dimensionless with the following definitions:
    Q_4=flow depletion/Q  t'=t*T/(S*L^2)  Lambda'=Lambda*L/T  Epsilon=S/Sigma  K=(K'/B')*L^2/T
    NOTE: Setting K'=0 for any value of epsilon gives the solution obtained by B.Hunt(1999)
    GROUND WATER,37(1),98-102 for either a completely confined or completely unconfined aquifer.

    :param t:
    :param K:
    :param epsilon:
    :param l:
    :return:
    """
    if t <= 0:
        Q_4 = 0
    else:
        n = 8
        Q_4 = 0
        for i in range(1, n + 1):
            Q_4 = Q_4 + Stehcoef(i, n) * g_4(i * math.log(2) / t, K, epsilon, l)  # todo check what log acutally is
        Q_4 = Q_4 * math.log(2) / t
    return Q_4

def g_4(p, K, epsilon, l):
    m0 = math.sqrt(p * (p + K * (1 + epsilon)) / (p + epsilon * K))
    g_4 = l * math.e**(-m0) / (p * (l + 2 * m0))
    return g_4

def Stehcoef(i, n):
    """
    'This calculates the coefficient c(i) in the Stehfest algorithm for Laplace transform
    'inversion. The integer n must be even.

    :param i:
    :param n:
    :return:
    """
    M = int(round(n / 2, 0))
    upperlimit = min([i, M])
    lowerlimit = int(math.floor((i + 1) / 2))
    Stehcoef = 0
    for K in range(lowerlimit, upperlimit+1):
        num = math.factorial(2 * K) * K**M
        denom = math.factorial(M - K) * math.factorial(K) * math.factorial(K - 1) * math.factorial(i - K) * math.factorial(2 * K - i)
        Stehcoef = Stehcoef + num / denom
    Stehcoef = Stehcoef * (-1)**(i + M)
    return Stehcoef

def hunt2003(discharge, time, trans, s, kaqt, baqt, sy, lam, l, return_rate=True):
    """
    compute the hunt 2003 stream depletion analysis note this does not account for application efficiency.
    the approach is vectorised, where a value is consistent across all entries it need only be specified once
    note setting kaqt to 0 gives hunt 1999

    :param discharge: pumping rate L/s
    :param time: time in days
    :param trans: transmissivity of the pumped aquifer m2/day
    :param s: storage coefficient of the pumped aquifer unitless
    :param kaqt: hydraulic conductivity of the aquitard (K') m/d
    :param baqt: thickness of the aquitard (B') m
    :param sy: specific yield of the aquifer unitless
    :param lam: stream bed conductance (bed conductance * bed width * bed thickness) m/day
    :param l: separation distance m
    :param return_rate: boolean if True return the rates false return the fraction
    :return:
    """
    discharge = np.atleast_1d(discharge)
    time=np.atleast_1d(time)
    idx = np.isclose(time,0)
    time[idx] = 1 # handle a weird vectorizing error
    trans=np.atleast_1d(trans)
    s=np.atleast_1d(s)
    kaqt=np.atleast_1d(kaqt)
    baqt=np.atleast_1d(baqt)
    sy=np.atleast_1d(sy)
    lam=np.atleast_1d(lam)
    l=np.atleast_1d(l)

    kb = kaqt/baqt
    t = time*trans/(s*l**2)
    k = kb*l**2/trans
    epsilon = s/sy
    lamd = lam*l/trans
    if return_rate:
        output = discharge * np.vectorize(Q_4)(t, k, epsilon, lamd)
    else:
        output = np.vectorize(Q_4)(t, k, epsilon, lamd)
    output[idx] = 0  # handle a weird vectorizing error
    return output


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    distances = [10**e for e in np.arange(0,5,0.02)]
    fig, ax = plt.subplots()
    for t in [7, 150, 365]:
        data = hunt2003(discharge=100,
                 time=t,
                 trans=1000,
                 s=0.1,
                 kaqt=1,
                 baqt=1,
                 sy= 0.1,
                 lam=1000,
                 l=distances,
                return_rate=False)
        ax.scatter(distances, data, label=t)
        ax.legend()
    plt.show()
