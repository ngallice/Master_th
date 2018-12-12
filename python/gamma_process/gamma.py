import numpy as np
import math as math

from scipy import stats
import scipy as sc

import random
import pickle,pprint

from scipy.stats import gamma
from scipy.optimize import minimize


def S(t, gamma, beta):
    return 1 - sc.special.gammainc(gamma, beta * t)


def S_sum(t, gamma, beta):
    summ = 0
    for i in range(gamma):
        summ += (beta * t) ** i / math.factorial(i)

    return summ * np.exp(-beta * t)


def P(t, gamma, beta):
    return (beta ** gamma) * (t ** (gamma - 1)) * np.exp(-beta * t) / sc.special.gamma(gamma)




def RHO(t, gamma, beta):
    summ = 0
    for i in range(gamma):
        summ += (beta * t) ** i / math.factorial(i)

    return beta ** gamma * t ** (gamma - 1) / (summ * math.factorial(gamma - 1))


def sum_exp(gamma, n, m):
    summ = 0
    for i in range(gamma):
        summ += i * np.exp(2 * np.pi * 1j * (n - m) * i / gamma)
    return summ


def phih(h, beta0, kappa, h_0):
    return beta0 * np.exp(kappa * (h - h_0))


def eigenvalue(n, gamma, beta):
    return beta * (np.exp(2 * np.pi * 1j * n / gamma) - 1)


def phi_0(n, gamma, beta):
    return beta * np.exp(2 * np.pi * 1j * n / gamma) / gamma


def Cnm(n, m, kappa):
    return kappa * (1 - np.exp(-2 * np.pi * 1j * n / gamma)) * (1 + sum_exp(gamma, n, m) / gamma)


def KLM(h, X, Y, beta0, kappa, h_0, tau_m, t):
    beta = phih(h, beta0, kappa, h_0)

    hpoint = 1 / tau_m * (-h + mu(t))

    alpha_r = beta * (np.cos(2 * np.pi / gamma) - 1)
    alpha_i = beta * np.sin(2 * np.pi / gamma)
    beta_r = (gamma + 1) * np.sin(np.pi / gamma) ** 2
    beta_i = (gamma + 1) * np.sin(np.pi / gamma) * np.cos(np.pi / gamma)
    gamma_r = 0.5
    gamma_i = 0.5 * np.tan(np.pi / gamma)
    eta_r = 1
    # eta_i=0

    K = (alpha_r + (beta_r + gamma_r) * hpoint) * X - (alpha_i + (beta_i - gamma_i) * hpoint) * Y + eta_r * hpoint
    L = (alpha_r + (beta_r - gamma_r) * hpoint) * Y + (alpha_i + (beta_i + gamma_i) * hpoint) * X
    M = 1 / tau_m * (-h + mu(t))

    return K, L, M


def KLMs(h, X, Y, beta0, kappa, h_0, tau_m, t):
    beta = phih(h, beta0, kappa, h_0)
    l1 = eigenvalue(1, gamma, beta)

    hpoint = 1 / tau_m * (-h + mu(t))

    K = l1.real * X - l1.imag * Y
    L = l1.imag * X + l1.real * Y
    M = 1 / tau_m * (-h + mu(t))

    return K, L, M


#PARAMETERS:

#SIMULATION PARAMETERS

dt=0.1
t_max=400
N_sim=int(t_max/dt)
tsim=[i*dt for i in range(N_sim)]

#----------------------------------------------------------

#NEURON PARAMETER:

gamma=25.7
tau_m=15 #ms
#----------------------------------------------------------

#RECOVERY FUNCTION PARAMETERS Phi(h)=nu_max*exp(beta(h-h0)):

kappa=1 #mv^(-1)
beta0=0.02*gamma #kHZ : firing rate for h=h_0
h_0=0 #mV
#----------------------------------------------------------



#EXTERNAL INPUT:  mu(t)=mu_1+epsilon*sin(omega*t)
#different frequency trials
fr_vec=np.logspace(-3,-1,5) #KHz
omega_vec=2*np.pi*fr_vec #kHz

#different epsilon trials
epsilon_vec=np.copy([1])# 5  #mV

#baseline of the external input
mu_1=h_0 #mV
#----------------------------------------------------------


tau_max=400
N=int(tau_max/dt)
tau_vec=[i*dt for i in range(N+1)]
S_init=[S(tau,gamma,beta0) for tau in tau_vec]

#TO STORE THE ACTIVITY and the external INPUT:

mu_Matrix  =np.zeros((len(omega_vec),len(epsilon_vec),N_sim))
A_Matrix   =np.zeros((len(omega_vec),len(epsilon_vec),N_sim))
A_int_Matrix   =np.zeros((len(omega_vec),len(epsilon_vec),N_sim))
As_Matrix   =np.zeros((len(omega_vec),len(epsilon_vec),N_sim))
#----------------------------------------------------------


for i_o, omega in enumerate(omega_vec):
    for i_e, epsilon in enumerate(epsilon_vec):
        print(i_o, i_e)


        def mu(t):
            return mu_1 - epsilon * np.sin(omega * t)


        # H PARMETER
        H = np.zeros((N_sim + 1, 1))
        H[0] = h_0

        A = np.zeros((N_sim, 1))
        A[0] = phih(H[0], beta0, kappa, h_0) / gamma

        X = np.zeros((N_sim + 1, 1))  # X in the thesis
        Y = np.zeros((N_sim + 1, 1))  # Y in the thesis

        Hs = np.zeros((N_sim + 1, 1))
        Hs[0] = h_0

        As = np.zeros((N_sim, 1))
        As[0] = phih(H[0], beta0, kappa, h_0) / gamma

        Xs = np.zeros((N_sim + 1, 1))  # X in the thesis
        Ys = np.zeros((N_sim + 1, 1))  # Y in the thesis

        # INTEGRAL:
        q_old = beta0 / gamma * np.copy(S_init)
        A_int = np.zeros((N_sim, 1))

        q_new = beta0 / gamma * np.copy(S_init)

        h_int = np.zeros((N_sim + 1, 1))
        h_int[0] = h_0

        for i in range(N_sim - 1):
            t = i * dt

            mu_Matrix[i_o, i_e, i] = mu(t)

            # INTEGRAL
            A_int[i] = q_new[0]
            q_new = np.zeros((N + 1, 1))
            hpoint = 1 / tau_m * (-h_int[i] + mu(t))
            h_int[i + 1] = h_int[i] + dt * hpoint
            beta_int = phih(h_int[i + 1], beta0, kappa, h_0)[0]

            for ir in range(N):
                rho_ = RHO(ir * dt, gamma, beta_int)
                q_new[ir + 1] = q_old[ir] * (1 - dt * rho_)
                q_new[0] += dt * rho_ * q_old[ir]

                #if q_old[ir] == 0:
                #   break

            q_new /= (np.sum(q_new * dt))
            q_old = q_new
            # ___END INTEGRAL______


            K1, L1, M1 = KLM(H[i], X[i], Y[i], beta0, kappa, h_0, tau_m, t)
            K2, L2, M2 = KLM(H[i] + 0.5 * dt * M1, X[i] + 0.5 * dt * K1, Y[i] + 0.5 * dt * L1, beta0, kappa, h_0, tau_m,
                             t)
            K3, L3, M3 = KLM(H[i] + 0.5 * dt * M2, X[i] + 0.5 * dt * K2, Y[i] + 0.5 * dt * L2, beta0, kappa, h_0, tau_m,
                             t)
            K4, L4, M4 = KLM(H[i] + dt * M3, X[i] + dt * K3, Y[i] + dt * L3, beta0, kappa, h_0, tau_m, t)

            X[i + 1] = X[i] + dt / 6 * (K1 + 2 * K2 + 2 * K3 + K4)
            Y[i + 1] = Y[i] + dt / 6 * (L1 + 2 * L2 + 2 * L3 + L4)
            H[i + 1] = H[i] + dt / 6 * (M1 + 2 * M2 + 2 * M3 + M4)

            K1s, L1s, M1s = KLMs(Hs[i], Xs[i], Ys[i], beta0, kappa, h_0, tau_m, t)
            K2s, L2s, M2s = KLMs(Hs[i] + 0.5 * dt * M1s, Xs[i] + 0.5 * dt * K1s, Ys[i] + 0.5 * dt * L1s, beta0, kappa,
                                 h_0, tau_m,
                                 t)
            K3s, L3s, M3s = KLMs(Hs[i] + 0.5 * dt * M2s, Xs[i] + 0.5 * dt * K2s, Ys[i] + 0.5 * dt * L2s, beta0, kappa,
                                 h_0, tau_m,
                                 t)
            K4s, L4s, M4s = KLMs(Hs[i] + dt * M3s, Xs[i] + dt * K3s, Ys[i] + dt * L3s, beta0, kappa, h_0, tau_m, t)

            Xs[i + 1] = Xs[i] + dt / 6 * (K1s + 2 * K2s + 2 * K3s + K4s)
            Ys[i + 1] = Ys[i] + dt / 6 * (L1s + 2 * L2s + 2 * L3s + L4s)
            Hs[i + 1] = Hs[i] + dt / 6 * (M1s + 2 * M2s + 2 * M3s + M4s)

            A[i + 1] = phih(H[i + 1], beta0, kappa, h_0) / gamma * (
                1 + X[i + 1] * np.cos(2 * np.pi / gamma) - Y[i + 1] * np.sin(2 * np.pi / gamma))

            As[i + 1] = phih(Hs[i + 1], beta0, kappa, h_0) / gamma * (
                1 + Xs[i + 1] * np.cos(2 * np.pi / gamma) - Ys[i + 1] * np.sin(2 * np.pi / gamma))

        A_Matrix[i_o, i_e, :] = A.flatten()
        A_int_Matrix[i_o, i_e, :] = A_int.flatten()
        As_Matrix[i_o, i_e, :] = As.flatten()



omega_search  = {'epsilon_vec': epsilon_vec,
               'omega_vec': omega_vec,
               'mu_Matrix': mu_Matrix,
               'A_Matrix': A_Matrix,
               'A_int_Matrix': A_int_Matrix,
                'As_Matrix': As_Matrix,
               'tsim':tsim,
                'fr_vec':fr_vec,
                 'mu_1':mu_1
              }

output = open('omega_searchgammareal.pkl', 'wb')

pickle.dump(omega_search, output)

output.close()