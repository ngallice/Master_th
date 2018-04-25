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



v=2.
a=0.1


t_window = 100
dt = 0.05
N = int(t_window / dt)


rho = np.zeros((N + 1, 1))

mu_ = 16


for r in range(N + 1):
    rho[r] = v*(1-np.exp(-a*r*dt))

int_rho = np.zeros((N + 1, 1))

for r in range(N):
    int_rho[r + 1] = int_rho[r] + (rho[r] + rho[r + 1]) / 2 * dt

exp_rho = np.exp(-int_rho[1:])


def int_tau(l,exp_rho):
    sum_=0
    for r in range(N):
        sum_+=np.exp(-l*dt*(r+0.5))*exp_rho[r]*(rho[r]+rho[r+1])/2*dt
    return abs(sum_)


L = np.zeros((N + 1, N + 1))
for i in range(N - 1):
    L[0, i + 1] = rho[i] * dt
    L[i + 1, i + 1] = -1 - rho[i] * dt
    L[i + 2, i + 1] = 1

L[0, 0] = -1
L[1, 0] = 1
L[0, N] = 1
L[N, N] = -1

eigv, v = np.linalg.eig(L)
eigv = eigv / dt


Grid_size=100

real_array=np.linspace(-2.5,0.,Grid_size)
im_array=np.linspace(0,2,Grid_size)

Grid=np.zeros((Grid_size,Grid_size))
for re in range(Grid_size):
    print(re)
    for im in range(Grid_size):
        Grid[re,im]=int_tau(real_array[re]+im_array[im]*1j,exp_rho)


plt.figure(figsize=(8,6))
plt.scatter(eigv.real,eigv.imag,s=20,marker='*',c='k')

max_ = 1.01
min_ = 0.99

cm = plt.cm.get_cmap('RdYlBu')

plt.figure(figsize=(8, 6))
for re in range(Grid_size):
    for im in range(Grid_size):
        if (Grid[re, im] >= min_) and (Grid[re, im] <= max_):
            sc = plt.scatter(real_array[re], im_array[im], s=30, c=abs(1 - Grid[re, im]), vmin=0, vmax=1 - min_,
                             cmap='spring')

plt.scatter(eigv.real, eigv.imag, s=20, marker='*', c='k')
plt.xlim(-1.5, 0.1)
plt.ylim(-0.1, 3)
plt.colorbar(sc)