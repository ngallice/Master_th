#! /usr/bin/env python

from pylab import *


N=100000  #number of neurons
tref=0.004
taum=0.02
mu1=0.8
mu2=5.
J=0.0
c=30.
du=1.
delay=0.005

taumax=2.#max ISI, in sec
dt=0.001
tauv=arange(0,taumax,dt)
K=len(tauv)
S=zeros(K)
rho=zeros(K)
P=zeros(K)
AAv=zeros(K)

m=zeros(K,dtype=int)   #refractory states, age vector m=q*dt*N

#random initial refractory distribution
# for i in range(N):
#     m[int(K*rand())]+=1
# print sum(m)

#synchronized initial conditions
m[0]=N

a=m[0]

#simulation time
tmax=0.3# sec
tv=arange(0,tmax,dt)

Nv=len(tv)
Av=zeros(Nv)
av=zeros(Nv)
dA=zeros(Nv)
print "time step dt=",dt," K=",K

v=zeros(K+1)
iref=int(tref/dt)
E1=exp(-dt/taum)
E2=1-E1


t=0.0
tv[0]=0.
Av[0]=(1.*m[0])/dt/N
av[0]=(1.*m[0])/dt/N
ion()
fig=figure(1)
clf()
ax1=axes((0.1,0.6,0.8,0.3))
line,=plot(tauv,m/(N*dt))
ylim((0,100.))
xlim((0,0.1))
ax2=axes((0.1,0.1,0.8,0.4))
line2,=plot(tv[:1],Av[:1])
#line3,=plot(tv[:1],av[:1])
xlim((0,tmax))
ylim((0,100.))
draw()
show()


jdelay=int(round(delay/dt))
#ion()
for j in range(1,Nv):
    if j*dt>0.1:
        mu=mu2
    else:
        mu=mu1
    Input=Av[j-1-jdelay]*J*E2
    v=roll(v,1)
    v[0]=0
    for i in range(iref,K):
        v[i+1]=mu+(v[i+1]-mu)*E1+Input

    n=0
    a=0.
    S=roll(S,1)
    S[0]=1
    for i in range(iref+1,K):
        rho[i]=c*exp((v[i]-1)/du)*dt
        pspike=1-exp(-rho[i])
        S[i]*=(1-pspike)
        P[i]=(S[i-1]-S[i])/dt
        if (m[i]>0):
            spikes=np.random.binomial(m[i],pspike)
        else:
            spikes=0
#        print m[i],spikes
        m[i]=m[i]-spikes
        n+=spikes

    m=roll(m,1)
    m[0]+=n
    AAv=roll(AAv,1)
    AAv[0]=(1.*m[0])/dt/N
    
    av[j]=sum(P*AAv)*dt
    Av[j]=(1.*m[0])/dt/N
    line.set_ydata(m/(N*dt))
    line2.set_xdata(tv[:j+1])
    line2.set_ydata(Av[:j+1])
#    line3.set_xdata(tv[:j+1])
#    line3.set_ydata(av[:j+1])
    fig.canvas.draw()

