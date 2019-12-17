#!/usr/bin/python
#coding:utf-8
#Pierre Puiseux, 14 décembre 2019
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
filter = signal.fftconvolve
filter = scipy.ndimage.convolve1d

"""
1 Lire le signal,
2 régler le paramètre freq, la fréquence d'echantillonnage en Hz(constant)
    = nb d'enregistrements par seconde. chez Aérotest freq=10 à ce jour (Dec. 2019)
Le signal S est une somme S = THF + HF + BF
    THF=bruit, très hautes fréquences
    BF=basses fréquences, 
    HF=hautes fréquences

3 filtrer le bruit (très hautes fréquences), pour cela régler le paramètre P_gauss, nb de points de convolution.
    Il faut que 2*P+1 points couvrent les parasites, mais pas plus. Tester plusieurs fois
4 régler la fenètre de temps 
5 Extraire les BF : elle sont obtenues par lissage du signal, avec un filtre gaussien. 
    Régler ce filtre avec P_gauss, a nouveau

That's all
"""
def filtreGaussien(P_gauss, epsilon):
    #filtre gaussien, largeur 2*P_gauss+1
    #filtre passe bas
    sigma=P_gauss/sqrt(-2.0*log(epsilon))
    som = 0.0
    b_gauss = zeros(2*P_gauss+1)
    for k in range(2*P_gauss+1):
        b_gauss[k] = exp(-(k-P_gauss)**2/(2*sigma**2))
        som += b_gauss[k]
    return b_gauss/som

def analyse(R, T, fenetre, pgauss, titre, show=True):
    if len(R) != len(T) :
        raise IndexError("R et T doivent être de même taille %d, != %d"%(len(R),len(T)))
    R0 = R.copy()#les data à analyser
    T0 = T.copy()
    if fenetre == 'all' : 
        t0, t1 = tmax, tmin = [T[0], T[-1]]
        n0, n1 = 0, len(T)
    else : 
        tmax, tmin = T[-1], T[0]
        t0, t1 = fenetre
        t0, t1 = max(tmin,t0), min(tmax,t1)
        n0, n1 = int(round(t0*freq,0))-1, int(round(t1*freq,0))-1
        R = R[n0:1+n1]#les data à analyser
        T = T[n0:1+n1]
    print('\nN=%d, len(R0)=%d, len(R)=%d, n0=%d,n1=%d'%(N, len(R0), len(R), n0, n1))
    #############################################################
    ## décomposition du signal filtré
    ## par extraction de basses fréquences avec le filtre gaussien 
    ##############################################################
    gaussien = filtreGaussien(P_gauss=pgauss, epsilon=0.01)
    mode = 'mirror'
    mode = 'reflect'
    # Rb = filter(R, gaussien, mode='full')
    # print('full', R.size, Rb.size)
    Rb = filter(R, gaussien, mode=mode)
    print(mode, R.size, Rb.size)
    if mode == 'valid' : 
        R  = R[pgauss:-pgauss]
        T  = T[pgauss:-pgauss]

    # exit()
    ##Les hautes fréquences Rh sont obtenues par soustraction des basses 
    ###fréquences au signal d'origine : Rh = R-Rb
    Rh = R - Rb#Rh = Les hautes frequences
    Ne = len(Rh)#nb pts echantillon
    Nes2 = int(1+len(Rh)/2)#nb pts echantillon
    ############################################################
    ## calcul du spectre HF et BF
    ############################################################
    hspectre = absolute(fft(Rh))/Ne
    # bspectre = absolute(fft(Rb))/Ne
    frequences = arange(Ne,)/(dt*Ne)

    ############################################################
    ## tracé
    ############################################################
    fig, (ax1,ax2) = plt.subplots(2, 1)
    fig.set_size_inches(20,10)
    ax1.set_title("La fenêtre")
    ax1.grid(True)
    ax1.plot(T0, R0, 'k-.', linewidth=0.2, label='signal complet')
    ax1.plot(T,Rh, label='HF')
    ax1.plot(T,Rb, label='BF')
    ax1.set_xlabel("temps (s)")
    ax1.set_ylabel("Force (N)")
    ax1.legend()

    ax2.set_title("Analyse spectrale HF")
    ax2.plot(frequences[:Nes2], hspectre[:Nes2])
    ax2.set_xlabel("f (Hz)")
    ax2.set_ylabel("A")
    ax2.grid(True)

    if titre =='' : 
        titre = 'Test en charge %s'%fin.name + '\n fenêtre [%.0f s, %.0f s],  durée=%.0f s, '%(t0,t1,t1-t0)
    plt.suptitle(titre)

    if show : plt.show()
    return R, T

def lireData(fin):
    with open(fin) as f : 
        lines = f.readlines()
    R = []#zeros(len(lines),float)
    N = len(lines)#nb d'enregistrements lus
    invalides, comm = [],[]
    for k,line in enumerate(lines) : 
        # print (line)
        if line[0] == '#':
            comm.append(line[1:])
            continue
        try : 
            w1, w2 = line.split(';')
            R.append(float(w1) + float(w2))
        except ValueError : 
            invalides.append[k]
            N -= 1
            continue
    R = asarray(R)
    return R, invalides, comm

if __name__ == '__main__' :
    ############################################
    ###Lecture data
    ############################################
    #freq = nb points par secondes
    freq, fin, pgauss, fenetres = (10, 
        Path('/Users/puiseux/Google Drive/Aerotest-FFVL/Aerotest/NervuresConfidentiel/2016-08-20-Stantick2.txt'),
        3, 
        (((15,30),50), ((33,43),100), ((48,66),100)))

    freq, fin, pgauss, fenetres = (50, 
        Path('/Users/puiseux/Google Drive/Aerotest-FFVL/Aerotest/ITV-Daytona-32/2019-12-05-ITV-Daytona-32.TXT'),
        20, 
        (((20,58),100), ((58, 78),200), ((78,85),100)))
    
    
    R00, invalides, comm = lireData(fin)
    print('\n\nFin lecture, %d enregistrements dont %d invalides.'%(len(R00)+len(invalides), len(invalides)))
    if comm :
        print('Commentaires : %s'%('\n'.join(comm)))
    ############################################
    ## freq = période d'échantillonnage eg en Hz 
    ##      = nb de points par sec.
    ############################################
    dt = 1/freq#secondes, période 0.1 s
    # freq = 1/dt#frequence = 10Hz
    N = len(R00)
    T00 = linspace(0, N*dt, N)
    print ('Durée du run : %.2f secondes'%T00[-1])
    print('len(R) = %d'%len(R00))
    ##Fin lecture
    # exit()
    #########################################
    ## suppression très hautes fréquences
    #On bouffe 10*dt=1 seconde à chaque bout
    #########################################
    titre0 = 'Test en charge %s'%fin.name + '\n %d valeurs, dt=%.2f s => durée=%.1f s, '%(N, dt, N/freq)
    titre = titre0 + '\n Filtrage du bruit : les HF doivent couvrir un large spectre (graphique du bas)...'
    # P_gauss = 3
    R0, T0 = analyse(R00, T00, fenetre='all', pgauss=pgauss,show=False, titre=titre)
    # recalage des temps
    #T0 -= T0[0]
    # P_gauss = 100
    print('tmin = %.2f, tmax=%.2f'%(T0[0],T0[-1]))
    for (t0, t1), pg in fenetres :
    # for (t0, t1), pgauss in ():
        ###########################################
        ##Fenetre de temps ((t0,t1),p)
        #   tmin < t0 < t1 <tmax 
        #   p = nb points de gauss
        ###########################################
        titre = fin.name + " : analyse fréquencielle de la fenêtre [%.1f, %.1f], pts de Gauss : %d"%(t0,t1, pgauss)
        try : 
            R1, T1 = analyse(R0, T0, fenetre=[t0,t1], pgauss=pg, show=False, titre=titre)
            print('tmin = %.2f, tmax=%.2f'%(T1[0],T1[-1]))
        except ValueError : 
            pass
    plt.show()
    exit()
