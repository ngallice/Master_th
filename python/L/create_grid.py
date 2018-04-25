import numpy as np

import pandas as pd
from pandas.plotting import scatter_matrix

from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns

from random import randrange

import re

import random
import pickle,pprint

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


U_t=15    # mV
D=2       # mV
C=1   # 10Hz-> 0.01ms^{-1}
dt=1    #ms
tau=20    #ms

def mu(t):
    return 16


t_window = 100
dt = 0.1
N = int(t_window / dt)

ur = np.zeros((N + 1, 1))
rho = np.zeros((N + 1, 1))

mu_ = 16

for r in range(N):
    ur[r + 1] = ur[r] * (1 - dt / tau) + mu_ * dt / tau

for r in range(N + 1):
    rho[r] = C * np.exp((ur[r] - U_t) / D)

int_rho = np.zeros((N + 1, 1))

for r in range(N):
    int_rho[r + 1] = int_rho[r] + (rho[r] + rho[r + 1]) / 2 * dt

exp_rho = np.exp(-int_rho[1:])


def int_tau(l,exp_rho):
    sum_=0
    for r in range(N):
        sum_+=np.exp(-l*dt*(r+0.5))*exp_rho[r]*(rho[r]+rho[r+1])/2*dt
    return abs(sum_)


Grid_size=500

real_array=np.linspace(-1.5,0.,Grid_size)
im_array=np.linspace(0,3,Grid_size)


Grid=np.zeros((Grid_size,Grid_size))
for re in range(Grid_size):
    print(re)
    for im in range(Grid_size):
        Grid[re,im]=int_tau(real_array[re]+im_array[im]*1j,exp_rho)


Grid_L = {'Grid': Grid}

output = open('Grid_laplace.pkl', 'wb')

# Pickle dictionary using protocol 0.
pickle.dump(Grid_L, output)


output.close()