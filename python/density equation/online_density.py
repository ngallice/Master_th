import matplotlib.pyplot as plt
from drawnow import drawnow
import numpy as np
from pylab import *

def makeFig():
    num_fig=2

    plt.subplot(num_fig, 1, 1)
    plt.plot(qt, 2000*q_new)
    plt.xlabel('$τ$ [ms]')
    plt.ylabel('$q(τ)$ [Hz]')
    plt.ylim(0, 60)
    plt.title('$t=$'+str(t)+'ms, $\mu(t)=$'+str(mu(t))+' mV')


    plt.subplot(num_fig, 1, 2)
    plt.plot(xt, A)
    plt.xlabel('t [ms]')
    plt.ylabel('A(t) [Hz]')
    plt.xlim(0,160)
    plt.ylim(0, 60)

    #plt.subplot(num_fig, 1, 3)
    #plt.plot(qt, u_old)
    #plt.xlabel('r [ms]')
    #plt.ylabel('u')

    #plt.subplot(num_fig, 1, 4)
    #plt.plot(xt, u)
    #plt.xlabel('t [ms]')
    #plt.ylabel('U')


plt.ion() # enable interactivity
fig=plt.figure() # make a figure


U_t=15    # mV
D=2 # 0.5       # mV
C=1   # 10Hz-> 0.01ms^{-1}
dt=0.5    #ms
tau=20    #ms

t_window=50


def mu(t):
    if t<70:  #ms
        return 8 #mV
    else :
        return 16 #mV



N = int(t_window / dt)

qt=[j*dt for j in range(N+1)]



#period=80
t_max=160
N_sim = int(t_max/dt)

u_old = np.zeros((N + 1, 1))
u_new = np.zeros((N + 1, 1))





q_old = np.zeros((N + 1, 1))
q_old[0] = 1



A =list()
xt=list()
u=list()


q_new = np.zeros((N + 1, 1))



for it in range(N_sim):

    t = it * dt


    if it%5==0:
        A.append( 2000*q_new[0])
        xt.append(t)
        u.append(u_new[N])
        drawnow(makeFig)





    q_new = np.zeros((N + 1, 1))
    u_new = np.zeros((N + 1, 1))

    for ir in range(N):
        u_new[ir+1 ] = u_old[ir] * (1 - dt / tau) + mu(t) * dt / tau

        rho = C * np.exp((u_old[ir] - U_t) / D)

        q_new[ir + 1] = q_old[ir] * (1 - dt * rho)
        q_new[0] += dt * rho * q_old[ir]

    u_new[N] = u_old[N] * (1 - dt / tau) + mu(t) * dt / tau
    rho = C * np.exp((u_old[N] - U_t) / D)
    q_new[N] += q_old[N] * (1 - dt * rho)
    q_new[0] += dt * rho * q_old[N]


    q_old = q_new
    u_old = u_new






