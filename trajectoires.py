#!/usr/bin/python
#coding:utf-8
import sys
from matplotlib import pyplot as plt
import numpy as np
from numpy import asarray
from math import pi, sqrt, sin
import scipy.integrate as integrate
class A(object):
    def __init__():
        pass


def Kg2N(x) : return x*9.81
def ms2Kmh(x) : return x*3.6

def toLaTeX(boite, bis) :
    hline = "\\hline"
    def dataLine(bi) : 
        name, fmax, dt, L, comm = bi
        gamma = L/(dt**2)
        V = ms2Kmh(gamma*dt)
        return "%s & %.0f & %.1f  & %.1f & %.1f & %s \\tabularnewline"%(name,fmax,dt,gamma,V,comm)
    
    lines = []
    lines = ["\\begin{tabular}{|l|l|l|l|l|l|}"]
    lines.append(hline)
    lines.append("\\noalign{\\vskip0.1cm}")

    lines.append("Modèle & $\\Delta F\\ \\left(N\\right)$ & $\\Delta t\\ \\left(s\\right)$ & $\\gamma \\left(m/s^2\\right)$ & $V_{max} (km/h)$ & Remarque \\tabularnewline")
    lines.append("\\noalign{\\vskip0.1cm}")
    lines.append(hline) 
    for bi in bis :
        # lines.append(hline) 
        lines.append(dataLine(bi))
    lines.append(hline) 
    
    lines.append("\\end{tabular}")
    return '\n'.join(lines)



def run4() :
    """
    """
    def Vmax(t1,t2,t3,L) : 
        return ms2Kmh(3*L/(2*t3+t2-t1))
    
    def toText(boite, L, bis) :
        print ("\n"+(len(boite)*'=') + '\n' + "%s"%boite+'\n' + len(boite)*'=' + '\n')

        hline = 4*' '+65*'-'
        print(hline)
        print('    Modèle               : t1   : t2   : t3   : Vmax    : Fmax    : 𝛾1    :')
        print(hline)
        for bi in bis : 
            name, t0, t1, t2, t3, fmax = bi
            t1 -= t0
            t2 -= t0
            t3 -= t0
            vmax = Vmax(t1,t2,t3, L)
            gam1 = vmax/t1
            # gam3 = -vmax/(t3-t2)
            print('    %-20s :'%name, end='')#---------------------
            print(" %-4.0f :"%(t1),end='')
            print(" %-4.0f :"%(t2),end='')
            print(" %-4.0f :"%(t3),end='')
            print(" %-7.1f :"%(vmax),end='')
            print(" %-7.1f :"%(fmax),end='')
            print(" %-5.1f :"%(gam1),end='')
            # print(" %-5.1f :"%(gam3),end='')
            print()
        print(hline)
    
    for boite, (L,bis) in D.items() : 
        toText(boite, L,    bis)
    print()
    print()
    # for boite, bis in D.items() :
    #     print() 
    #     print() 
    #     print(toLaTeX(boite, bis))

class Trajectoire():
    """
    - A est l'ampitude relative des oscillations : F_max*A~4000N, doit être à peu pres constant, 
        à lire sur les diagrammes
    - les k_i sont des entiers de sorte que au bout des intervalles les F_i  et les V_i se raccordent
    """
    params = dict(options='sqrt+sin',om=2, A=0.05, util=0.75)
    # def __init__(self,operator,L,name,t0,t1,t2,t3,Fmax)
    def __init__(self,**kargs):
        for arg, value in self.params.items():
            setattr(self, arg, value)
        
        for arg, value in kargs.items():
            setattr(self, arg, value)
        self.t1 -= self.t0
        self.t2 -= self.t0
        self.t3 -= self.t0
        if self.A==0.0 or '0' in self.options : 
            self._Cv = self._Cv0
        elif 'sin' in self.options:
            self._Cv = self._CvSin  
        else : 
            raise NotImplementedError
        self.racine = 'sqrt' in self.options
        self.k1 = self.om*self.t1
        self.k2 = self.om*(self.t2-self.t1)
        self.k3 = self.om*(self.t3-self.t2)
        self.L *= self.util

    def __str__(self):
        print('type(self.k1)',type(self.k1))
        msgs = [
        ('Opérateur', self.operator),
        ('Modèle',self.name),
        ('Fmax (N)',self.Fmax),
        ('Vmax (km/h)',"%.2f"%(3.6*self.Vmax)),
        ('longueur piste (m)',self.L),
        ('longueur utile (m)',self.L*self.util),
        ('A, om', (self.A,self.om)),
        ('(k1, k2, k3)', (round(self.k1/(2*pi),1),round(self.k2/(2*pi),1),round(self.k3/(2*pi),1))),
        ('(t1, t2, t3)(s)', (self.t1,self.t2,self.t3)),
        ]
        return "<Trajectoire> :\n"+'\n'.join(['%20s = %s'%(key,value) for key,value in msgs])

    def Cv(self,t):
        """
        Coef de vitesse = V/Vmax;
        """
        Cv2 = self._Cv(t) 
        if self.racine : 
            return sqrt(Cv2) if Cv2>0 else 0.0
        else :
            return Cv2

    def F(self,t):
        return self.Fmax*self._Cv(t)

    def _CvSin2(self,t):
        if 0.0 <= t <= self.t1 :
            return t/self.t1+self.A*sin(self.k1*(t/self.t1))
        elif self.t1 <= t <= self.t2:
            return 1+self.A*sin(2*self.k2*((t-self.t1)/(self.t2-self.t1)))
        elif self.t2 <= t <= self.t3:
            return (self.t3-t)/(self.t3-self.t2)+self.A*sin(self.k3*((self.t3-t)/(self.t3-self.t2)))
        else : 
            return 0.0

    def _CvSin(self,t):
        if 0.0 <= t <= self.t1 :
            return t/self.t1+self.A*sin(self.k1*(t/self.t1))
        elif self.t1 <= t <= self.t2:
            return 1+self.A*sin(self.k2*((t-self.t1)/(self.t2-self.t1)))
        elif self.t2 <= t <= self.t3:
            return (self.t3-t)/(self.t3-self.t2)+self.A*sin(self.k3*((self.t3-t)/(self.t3-self.t2)))
        else : 
            return 0.0
   
    def _Cv0(self,t):
        if 0.0 <= t <= self.t1 :
            return t/self.t1
        elif self.t1 <= t <= self.t2:
            return 1
        elif self.t2 <= t <= self.t3:
            return (self.t3-t)/(self.t3-self.t2)
        else : 
            return 0.0
    
    @property
    def Vmax(self) :
        if not hasattr(self,'_Vmax') : 
            self._Vmax = self.L/self.__T(self.t3)
        return self._Vmax
    
    def X(self,t) :
        return self.Vmax*self.__T(t)
    
    def __T(self,t):
        """La trajectoire X(t) divisée par Vmax, homogene à un temps"""
        if not 0.0 <= t <= self.t3 : 
            return np.nan
        if t < self.t1 :
            return integrate.quad(self.Cv,0,t,limit=100)[0]
        elif t < self.t2:
            T1 = integrate.quad(self.Cv,0,self.t1,limit=100)[0]
            T2 = integrate.quad(self.Cv,self.t1,t,limit=100)[0]
            return T1+T2
        else:
            T1 = integrate.quad(self.Cv,0,self.t1,limit=100)[0]
            T2 = integrate.quad(self.Cv,self.t1,self.t2,limit=100)[0]
            T3 = integrate.quad(self.Cv,self.t2,t,limit=100)[0]
            return T1+T2+T3

    def plotF(self):
        n = 1000
        T = np.linspace(0.0, self.t3, n)
        F = asarray([self.F(t) for t in T])
        fig, ax = plt.subplots()
        plt.plot(T,F)
        ax.set_title("%s, modélisation de la force appliquée"%self.name)
        plt.plot([0.0, self.t1,self.t2,self.t3],[0.0,self.Fmax,self.Fmax,0],
            'k--',linewidth=0.4)
        # plt.plot([self.t1],[0.0,self.Fmax])
        plt.plot([self.t1, self.t1],[0.0,self.Fmax])
        plt.plot([self.t2, self.t2],[0.0,self.Fmax])
        ax.grid(True)
        plt.xlabel('$time (s)$')
        plt.ylabel('$force (N)$')
        plt.show()
    
    def plotV(self):
        n = 1000
        T = np.linspace(0.0, self.t3, n)
        V = asarray([self.Cv(t) for t in T])
        # fig, ax = plt.subplots()
        plt.title("%s, modélisation de la vitesse"%self.name)
        plt.grid(True)
        plt.plot(T,V)
        self.Vmax = 1.0
        plt.plot([0.0, self.t1,self.t2,self.t3],[0.0,self.Vmax,self.Vmax,0],
            'k--',linewidth=0.4)
        # plt.plot([self.t1],[0.0,self.Fmax])
        plt.plot([self.t1, self.t1],[0.0,self.Vmax])
        plt.plot([self.t2, self.t2],[0.0,self.Vmax])
        # ax.grid(True)
        plt.xlabel('$time (s)$')
        plt.ylabel('coefficient de vitesse (sans dimension)')
        plt.show()
    
    def plotX(self):
        n = 100
        T = np.linspace(0.0, self.t3, n)
        X = asarray([self.X(t) for t in T])
        print('X=',X.shape)
        # fig, ax = plt.subplots()
        plt.title("%s, modélisation de la trajectoire ($V_{max}=%.1f km/h$)"%(self.name,3.6*self.Vmax))
        plt.grid(True)
        plt.plot(T,X)
        # self.Vmax = 1.0
        # plt.plot([0.0, self.t1,self.t2,self.t3],[0.0,self.Vmax,self.Vmax,0],
            # 'k--',linewidth=0.4)
        # plt.plot([self.t1],[0.0,self.Fmax])
        plt.plot([self.t1, self.t1],[0.0,self.X(self.t1)],'k--')
        plt.plot([self.t2, self.t2],[0.0,self.X(self.t2)],'k--')
        plt.plot([self.t3, self.t3],[0.0,self.X(self.t3)],'k--')
        # ax.grid(True)
        plt.xlabel('$time (s)$')
        plt.ylabel('distance (m)')
        plt.show()
    
    def plot(self):
        n = 1000
        T = np.linspace(0.0, self.t3, n)
        fig, axs = plt.subplots()
        plt.subplot(221)
        plt.grid(True)
        F = asarray([self.F(t) for t in T])
        plt.plot(T,F)
        plt.plot([0.0, self.t1,self.t2,self.t3],[0.0,self.Fmax,self.Fmax,0],
            'k--',linewidth=0.4)
        plt.plot([self.t1, self.t1],[0.0,self.Fmax],'r--',linewidth=0.6)
        plt.plot([self.t2, self.t2],[0.0,self.Fmax],'r--',linewidth=0.6)
        plt.ylabel('$force (N)$')

        plt.subplot(222)
        plt.grid(True)
        V = asarray([self.Vmax*self.Cv(t) for t in T])
        plt.plot(T,V)
        plt.plot([0.0, self.t1,self.t2,self.t3],[0.0,self.Vmax,self.Vmax,0],
            'k--',linewidth=0.4)
        plt.plot([self.t1, self.t1],[0.0,self.Vmax],'r--',linewidth=0.6)
        plt.plot([self.t2, self.t2],[0.0,self.Vmax],'r--',linewidth=0.6)
        plt.xlabel('$time (s)$')
        plt.ylabel('coefficient de vitesse $C_V$')

        plt.subplot(223)
        plt.grid(True)
        X = asarray([self.X(t) for t in T])
        plt.plot(T,X)
        plt.plot([self.t1, self.t1],[0.0,self.X(self.t1)],'r--',linewidth=0.6)
        plt.plot([self.t2, self.t2],[0.0,self.X(self.t2)],'r--',linewidth=0.6)
        plt.plot([self.t3, self.t3],[0.0,self.X(self.t3)],'r--',linewidth=0.6)
        plt.ylabel('distance (m)')
        plt.xlabel('$time (s)$')
        plt.ylabel('$X (m)$')
        fig.suptitle('%s : modélisation force, vitesse, trajectoire'%self.name, fontsize=16)


        plt.show()
if __name__ == '__main__' :
    D = dict(
        Turquoise = [900.00,#longueur de la piste
            #(nom, t0,t1,t2,t3,fmax (N))
            [('BiGolden 4 light', 43, 72, 77, 83, 16797.),
            ('Windtech Ru-Bi 2',  48, 78, 82, 92, 17739.59, ),
            ('Supair Sora 2',     18, 43, 51, 57, 17360.56),
            ('Supair Sora 1-41',  10, 55, 58, 62, 1739*9.81),]
                    ],
        Aérotest = [1700.00,
            [('Hercules',        6, 54, 61, 63, 22840),
            ('MacPara Trike 42', 7, 51, 56, 58, 20850),
            ('Stewart',          0, 36, 43, 45, 1760*9.81),
            ('Stewart-DGAC',     7, 27, 34, 35, 2500*9.81),
            ('Shuttle',         15, 49, 59, 60, 2500*9.81),]
                    ]
            )
    P = dict(options='sin+sqrt',om=2, A=0.08, util=0.75)
    AeroTest = [
                Trajectoire(operator='Aérotest', L=1700, name='Stewart-DGAC',
                    t0=7, t1=27, t2=34, t3=45, Fmax=2500*9.81,**P),
                Trajectoire(operator='Aérotest', L=1700, name='Hercules',
                    t0=6, t1=54, t2=61, t3=70, Fmax=22840, **P),
                Trajectoire(operator='Aérotest', L=1700, name='MacPara Trike 42',
                    t0=7, t1=51, t2=56, t3=70, Fmax=20850, **P),
                Trajectoire(operator='Aérotest', L=1700, name='Stewart',
                    t0=0, t1=36, t2=43, t3=55, Fmax=1760*9.81, **P),
                Trajectoire(operator='Aérotest', L=1700, name='Shuttle',
                    t0=15, t1=49, t2=59, t3=70, Fmax=2500*9.81, **P)
                ]

    P = dict(options='sin+sqrt',om=2, A=0.08, util=0.90)
    Turquoise = [
                Trajectoire(operator='Turquoise',L=900, name='Supair Sora 2',
                    t0=18, t1=43, t2=51, t3=57, Fmax=17360.56, **P),
                Trajectoire(operator='Turquoise',L=900, name='Supair Sora 1-41',
                    t0=10, t1=55, t2=58, t3=62, Fmax=1739*9.81, **P),
                Trajectoire(operator='Turquoise',L=900, name='BiGolden 4 light',
                    t0=43, t1=72, t2=77, t3=83, Fmax=16797., **P),
                Trajectoire(operator='Turquoise',L=900, name='Windtech Ru-Bi 2',
                    t0=48, t1=78, t2=82, t3=92, Fmax=17739.59, **P),
                ]
    for T in AeroTest :
        print("%s"%T)
        # print(T.X(T.t3))
        if 0 : T.plot()
    for T in Turquoise :
        print("%s"%T)
        # print(T.X(T.t3))
        if 1 : T.plot()


