import matplotlib.pyplot as plt
from drawnow import drawnow
import numpy as np
from pylab import *

def makeFig():
    #plt.figure(figsize=(8, 8))
    plt.plot(t_, v_0, label='$v_0$')
    #plt.plot(t_, v_0t, label='$v^{*}_0$')
    plt.plot(t_, v_1, label='$v_1$')
    plt.plot(t_, v_2, label='$v_2$')

    #plt.plot(t_q, q[:, time], c='k', ls=':', label='$q_{steady}$')

    plt.legend()
    plt.xlim(0, t_window)

    plt.xlabel('r [ms]')
    plt.ylabel('$q_0$')
    plt.title('$dt=$' + str(t))


plt.ion() # enable interactivity
fig=plt.figure() # make a figure


t_window=50
dt=0.1

t_final=200
N=int(t_window/dt)
N_sim=int(200/dt)

tau_m=20
gamma=0.01




t_ = [j * dt for j in range(N + 1)]



L=np.load('Lt.npy')

e0=[]
e1=[]
e2=[]
et=[]



for it in range(N_sim):

    t = it * dt

    # Phi


    #drawnow(makeFig)

    if it%50==0:
        eigv, vec = np.linalg.eig(L[:, :, it])
        eigv = eigv / dt

        real_eigv = eigv.real
        x = np.sort(real_eigv)

        eig_0 = np.where(eigv.real == x[-1])
        v_0 = vec[:, eig_0[0][0]].real
        e0.append(x[-1])

        if v_0[0] < 0:
            v_0 = -1 * v_0

        norm = np.sum(v_0) * dt
        v_0 = v_0 / norm

        eig_1 = np.where(eigv.real == x[-2])
        v_1 = vec[:, eig_1[0][0]].real
        e1.append(x[-2])

        if np.max(v_1) != np.max(abs(v_1)):
            v_1 = -1 * v_1

        eig_2 = np.where(eigv.real == x[-4])
        v_2 = vec[:, eig_2[0][0]].real
        e2.append(x[-4])

        if np.max(v_2) != np.max(abs(v_2)):
            v_2 = -1 * v_2

        drawnow(makeFig)

        et.append(t)




plt.figure(figsize=(8, 8))
plt.plot(et, e0, label='$\lambda_0$')
plt.plot(et, e1, label='$\lambda_1$')
plt.plot(et, e2, label='$\lambda_2$')


plt.legend()
plt.xlim(0, t_window)

plt.xlabel('t [ms]')
plt.ylabel('$\lambda$')
plt.show()







