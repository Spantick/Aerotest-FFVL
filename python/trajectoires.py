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
        lines.append(dataLine(bi))
    lines.append(hline) 
    
    lines.append("\\end{tabular}")
    return '\n'.join(lines)

class Trajectoire():
    #Valeurs par défaut
    # params = dict(options='sin', om=2, A=0.05)
    params = dict(options='sin', om=2, A=0.08)
    # def __init__(self,operator,L,name,t0,t1,t2,t3,Fmax)
    def __init__(self,**kargs):
        for arg, value in self.params.items():
            setattr(self, arg, value)
        
        for arg, value in kargs.items():
            setattr(self, arg, value)
        #on se ramène à t0 = 0
        self.t1 -= self.t0
        self.t2 -= self.t0
        self.t3 -= self.t0
        if self.A==0.0 or 'lin' in self.options : 
            #modelisation linéaire
            self._Cv = self._CvLin
        elif 'sin' in self.options:
            #modelisation sinusoide+droite
            self._Cv = self._CvSin  
        else : 
            raise NotImplementedError
        #le nombre de pic par palier
        self.k1 = self.om*self.t1
        self.k2 = self.om*(self.t2-self.t1)
        self.k3 = self.om*(self.t3-self.t2)

    def Cv(self,t):
        """
        Coef de vitesse = V/Vmax;
        """
        Cv =  self._Cv(t)
        return sqrt(Cv) if Cv>0 else 0.0

    def F(self,t):
        return self.Fmax*self._Cv(t)

    def _CvSin(self,t):
        if 0.0 <= t <= self.t1 :
            return t/self.t1 + self.A*sin(self.k1*(t/self.t1))
        elif self.t1 <= t <= self.t2:
            return 1 + self.A*sin(self.k2*((t-self.t1)/(self.t2-self.t1)))
        elif self.t2 <= t <= self.t3:
            return (self.t3-t)/(self.t3-self.t2) + self.A*sin(self.k3*((self.t3-t)/(self.t3-self.t2)))
        else : 
            return 0.0
   
    def _CvLin(self,t):
        """Vitesse linéaire """
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

        plt.subplot(224)
        plt.grid(False)
        plt.axis([0, 10, 0, 10])
        msg = '\n'.join(self.LaTeX())
        print(msg)
        plt.text(0.3,0.5,msg,{'color': 'black', 'fontsize': 12}, horizontalalignment='left',clip_on=True, multialignment="left")
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        fig.suptitle('%s : modélisation force, vitesse, trajectoire'%self.name, fontsize=16)


        plt.show()

    def LaTeX(self) : 
        msgs = [
        ('Opérateur', self.operator,''),
        ('Modèle',self.name,''),
        ('$F_{max}$',self.Fmax,'$N$'),
        ('$V_{max}$',"%.2f"%(3.6*self.Vmax),'$km/h$'),
        ('$L_u$',self.L,'$m$'),
        ('$(t_0, t_1, t_2, t_3)$', (self.t0, self.t1+self.t0,self.t2+self.t0,self.t3+self.t0),'$s$'),
        ]
        if 'lin' in self.options or self.A == 0:
            msgs.insert(3,('Modélisation', 'linéaire par palier', ''))
        else:#if 'sin' in self.options and self.A != 0:
            msgs.insert(3,('Modélisation', 'sinusoïdale par palier', ''))
            msgs.insert(6,('Coeff d\'Amplitude $\\sin$', self.A, ''))
            msgs.insert(8,('$\\omega$', self.om, '$s^{-1}$'))
        # ('$(k_1, k2, k3)$', (round(self.k1/(2*pi),1),round(self.k2/(2*pi),1),round(self.k3/(2*pi),1))),
        return ['%s : %s %s'%(key,value,unit) for key,value,unit in msgs]

    def __str__(self):
        msgs = [
        ('Opérateur', self.operator),
        ('Modèle',self.name),
        ('Fmax (N)',self.Fmax),
        ('Vmax (km/h)',"%.2f"%(3.6*self.Vmax)),
        ('longueur utile (m)',self.L),
        ('(t1, t2, t3)(s)', (self.t1,self.t2,self.t3)),
        ]
        if 'lin' in self.options or self.A == 0:
            msgs.insert(2,('Modélisation', 'Linéaire par palier'))
        if 'sin' in self.options and self.A != 0:
            msgs.insert(2,('Modélisation', 'Sinusoïdale par palier'))
            msgs.insert(6,('A, om', (self.A,self.om)))
            msgs.insert(7,('(k1, k2, k3)', (round(self.k1/(2*pi),1),round(self.k2/(2*pi),1),round(self.k3/(2*pi),1))))

        return "<Trajectoire> :\n"+'\n'.join(['%20s = %s'%(key,value) for key,value in msgs])


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
    """
    - option est le type de modélisation de F(t) : linéaire ou linéaire+sinusoïde
    - om la fréquence ? pulsation ? pour les sinusoïdes sin(om*t)
    - A est l'ampitude relative des oscillations : F_max*A~4000N, doit être à peu pres constant, 
        à lire sur les diagrammes.
    - operateur = 'Aérotest' ou 'Turquoise' ou autre
    - L est la longueur utile de la piste pour les 3 paliers "montée en charge", "plateau", "décélération"
    - name est le nom du modèle
    - t0, t1, t2, t3 sont les temps limites des trois paliers, à lire sur les diagrammes. t0 peut être non nul
    - Fmax est la charge exigée par la norme
    """
    #valable pour tous
    P0 = dict(options='sin', om=4, A=0.125)

    #Pour Aérotest uniquement
    P = dict(operator='Aérotest', L=1600)
    AeroTest = [
                Trajectoire(name='Hercules',         t0=6,  t1=54, t2=61, t3=70, Fmax=round(22840,1), **P, **P0),
                Trajectoire(name='MacPara Trike 42', t0=7,  t1=51, t2=56, t3=70, Fmax=round(20850,1), **P, **P0),
                Trajectoire(name='Stewart-DGAC',     t0=7,  t1=27, t2=34, t3=45, Fmax=round(2500*9.81,1), **P, **P0),
                Trajectoire(name='Stewart',          t0=0,  t1=36, t2=43, t3=55, Fmax=round(1760*9.81,1), **P, **P0),
                Trajectoire(name='Shuttle',          t0=15, t1=49, t2=59, t3=70, Fmax=round(2500*9.81,1), **P, **P0)
                ]

    #Pour Turquoise uniquement
    P = dict(L=800, operator='Turquoise')
    Turquoise = [
                Trajectoire(name='Supair Sora 2',    t0=18, t1=43, t2=51, t3=57, Fmax=round(17360.56,1), **P, **P0),
                Trajectoire(name='Supair Sora 1-41', t0=10, t1=55, t2=58, t3=62, Fmax=round(1739*9.81,1), **P, **P0),
                Trajectoire(name='BiGolden 4 light', t0=43, t1=72, t2=77, t3=83, Fmax=round(16797.,1), **P, **P0),
                Trajectoire(name='Windtech Ru-Bi 2', t0=48, t1=78, t2=82, t3=92, Fmax=round(17739.59,1), **P, **P0),
                ]
    
    for i, T in enumerate(AeroTest) :
        pass
        print("%s"%T)
        if i>=0 : T.plot()
    for i, T in enumerate(Turquoise) :
        print("%s"%T)
        if i>=0 : T.plot()


