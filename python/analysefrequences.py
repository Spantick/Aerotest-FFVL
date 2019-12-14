#!/usr/bin/python
#coding:utf-8
import sys
from matplotlib import pyplot as plt
import numpy as np, scipy
from numpy import asarray, zeros, linspace, absolute, arange, zeros_like, ones
from math import pi, sqrt, sin
import scipy.integrate as integrate
# import scipy.signal.convolve as convolve
from matplotlib import pyplot as plt
from path import Path
from scipy.fftpack import *
from numpy.fft import rfft
from scipy import signal
from math import sqrt, exp,log
"""
Pour filtrer un signal et plein d'autres choses : a lire absolument
http://www.f-legrand.fr/scidoc/docimg/numerique/filtre/introfiltres/introfiltres.html
"""
filter = signal.convolve

# class A(object):
#     def __init__():
#         pass
# def Kg2N(x) : return x*9.81
# def ms2Kmh(x) : return x*3.6

def filtreGaussien(P_gauss, epsilon):
    #filtre gaussien, largeur 2*P_gauss+1
    sigma=P_gauss/sqrt(-2.0*log(epsilon))
    som = 0.0
    b_gauss = zeros(2*P_gauss+1)
    for k in range(2*P_gauss+1):
        b_gauss[k] = exp(-(k-P_gauss)**2/(2*sigma**2))
        som += b_gauss[k]
    return b_gauss/som
###Lecture data
fin = Path('/Users/puiseux/Google Drive/Aerotest-FFVL/Aerotest/NervuresConfidentiel/2016-08-20-Stantick2.txt')
with open(fin) as f : 
    lines = f.readlines()
# print (lines)
# print (len(lines))
R00 = zeros(len(lines),float)
for k,line in enumerate(lines) : 
    # print (line)
    if line[0] == '#':continue
    try : 
        w1,w2 = line.split(';')
        R00[k] = float(w1) + float(w2)
    except ValueError : 
        # print(line)
        pass
##Fin lecture
#########################################
## suppression très hautes fréquences
#On bouffe 10*dt=1 seconde à chaque bout
#########################################
P_gauss = 3
plt.plot(R00[P_gauss:-P_gauss], label='original')
#moyenne sur dt*(2*P_gauss+1) secondes => 0.5 secondes si P_gauss=2
gaussien = filtreGaussien(P_gauss=P_gauss, epsilon=0.01)
R00 = filter(R00, gaussien, mode='valid')#[P_gauss:-P_gauss]
plt.plot(R00, label='filtré (gauss(%d)'%P_gauss)
plt.legend()
# plt.show()
# exit()
# Rh = R0[P_gauss:-P_gauss]-Rb#Rh = Les hautes frequences
# T = T[P_gauss:-P_gauss]

## Extraction tranche utile
N = len(R00)
dt = 0.1#secondes, pulsation 0.1 s
freq = 1/dt#frequence = 10Hz
T00 = linspace(0, N*dt,N)
T0, T1 = T00[0], T00[-1]
###########################################
fenetre = (t0, t1) = (50,65)#secondes
###########################################
n0, n1 = int(round(t0/dt,0)), int(round(t1/dt,0))
R0 = R00[n0:1+n1]#les data à analyser
T = T00[n0:1+n1]
print('\nN, R0.shape, n0,n1', str((N, R0.shape, n0,n1)))
# exit()

R = R0#.copy()
P_gauss = 15
gaussien = filtreGaussien(P_gauss=P_gauss, epsilon=0.01)
Rb = filter(R, gaussien, mode='valid')
Rh = R0[P_gauss:-P_gauss]-Rb#Rh = Les hautes frequences
T = T[P_gauss:-P_gauss]
#On re-filtre les très hautes fréquences (bruit)
# gaussien = filtreGaussien(P_gauss=P_gauss, epsilon=0.01)
# R = filter(Rh, gaussien, mode='valid')
# Rh = Rh[P_gauss:-P_gauss]-R
# Rb = Rb[P_gauss:-P_gauss]
# Rh = Rh[n0:n1]
# Rb = Rb[n0:n1]
Ne = len(Rh)
# T = t0+arange(Ne)*dt
# x, y = T[:,0], T[:,1]
# u = x+y
#Te = T[1]-T[0]
print('dt',dt)
hspectre = absolute(fft(Rh))/Ne
frequences = arange(Ne,)/(dt*Ne)
bspectre = absolute(fft(Rb))/Ne
# bfrequences = arange(Ne,)/(dt*Ne)
# FFT=fft(c)
# print(FFT)
fenetre = []
nf0, nf1 = 0,0
fig, ((sbax,shax),(fbax,fhax)) = plt.subplots(2, 2)
sbax.set_title("Basses fréquences")
sbax.grid(True)
sbax.plot(T,Rb)
sbax.set_xlabel("temps (s)")
sbax.set_ylabel("Force (N)")
fbax.plot(frequences, bspectre)
fbax.set_xlabel("f (Hz)")
fbax.set_ylabel("A")
fbax.grid(True)

shax.set_title("Hautes fréquences")
shax.grid(True)
shax.plot(T,Rh)
shax.set_xlabel("temps (s)")
shax.set_ylabel("Force (N)")
fhax.plot(frequences, hspectre)
fhax.set_xlabel("f (Hz)")
fhax.set_ylabel("A")
fhax.grid(True)
plt.show()
