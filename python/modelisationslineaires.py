#!/usr/parapn/python
#coding:utf-8
import sys
class A(object):
    def __init__():
        pass

def Kg2N(x) : return x*9.81
def ms2Kmh(x) : return x*3.6

def run0():
    """Modélisation de la charge : linéaire à un seul palier"""
    D = dict(
        Turquoise = [
            (450,600),#longueur de piste utile
            #(nom, fmax (N), dt(s))
            ('BiGolden 4 light', 16797.33, 72-42),
            ('Windtech Ru-Bi 2', 17739.59, 77-48),
            ('Supair Sora 2',    17360.56, 46-15),
                    ],
        Aérotest = [(1150,1300),
            ('Hercules',         22840, 56-6),
            ('MacPara Trike 42', 20850, 52-6),
                    ]
            )

    def Vmax(tm,Lu) :
        return 1.5*(Lu/tm)

    def toLaTeX(operateur, valeurs) :
        Lu0,Lu1 = valeurs[0]#longueurs utiles de piste
        hline = "\\hline"
        def dataLine(parap) :
            name, fmax, tmax = parap
            vmax0 = Vmax(tmax,Lu0)
            vmax1 = Vmax(tmax,Lu1)
            Vs = "%3.0f - %-3.0f"%(ms2Kmh(vmax0),ms2Kmh(vmax1))
            Ls = "%4.0f - %-4.0f"%(Lu0,Lu1)
            gammas = "%2.1f - %2.1f"%(vmax0/tmax, vmax1/tmax)
            return "%s & %.0f & %.0f  & %s & %s & %s \\tabularnewline"%(name,fmax,tmax,Ls,Vs,gammas)
        
        lines = []
        lines = ["\\begin{tabular}{|l|l|l|l|l|l|}"]
        lines.append(hline)
        lines.append("\\noalign{\\vskip0.1cm}")

        lines.append("Modèle & $F_{max} \\left(N\\right)$ & $t_{max} \\left(s\\right)$ & $L \\left(m\\right)$ & $V_{max} \\left(km/h\\right)$ & $\\gamma \\left(m/s^2\\right)$ \\tabularnewline")
        #             Modèle & $F_{max}  \left(N\right)$  & $t_{max} \left(s\right)$ & $ L \left(m\right)$ & $V_{max} (km/h)$ & $\gamma (m/s^2)$ \tabularnewline
        lines.append("\\noalign{\\vskip0.1cm}")
        lines.append(hline) 
        for parap in valeurs[1:] :
            # lines.append(hline) 
            lines.append(dataLine(parap))
        lines.append(hline) 
        
        lines.append("\\end{tabular}")
        return '\n'.join(lines)

    def toText(operateur, valeurs) :
        Lu0,Lu1 = valeurs[0]#longueurs utiles de piste
        print ("\n"+(len(operateur)*'=') + '\n' + "%s"%operateur+'\n' + len(operateur)*'=' + '\n')

        separateur = 4*' '+81*'-'
        print(separateur)
        print('    : Modèle               : Fmax (N) : tmax (s) :    L (m)    :Vmax (km/h): 𝛾 (m/s2) :')
        print(separateur)
        # print(valeurs[2:])
        for parap in valeurs[1:] : 
            name, fmax, tmax = parap
            print('    : %-20s :'%name, end='')
            print(" %-8.0f :"%(fmax),end='')
            print(" %-8.1f :"%(tmax),end='')
            
            # print(" %-7.2f :"%(gamma),end='')
            vmax0 = Vmax(tmax,Lu0)
            vmax1 = Vmax(tmax,Lu1)
            print(" [%-4.0f;%-4.0f] :"%(Lu0,Lu1),end='')
            print(" [%-3.0f;%-3.0f] :"%(ms2Kmh(vmax0),ms2Kmh(vmax1)),end='')
            print(" [%2.1f;%2.1f] :"%(vmax0/tmax, vmax1/tmax),end='')

            print()
        print(separateur)
    
    for operateur, valeurs in D.items() : 
        # print (operateur)
        # print(valeurs)
        # toText(operateur, valeurs)
        print(toLaTeX(operateur, valeurs))
    # print()
    # for operateur, valeurs in D.items() :
    #     print() 
    #     print() 
    #     print(toLaTeX(operateur, valeurs))

def run1():
    """Modélisation de l'accélération : linéaire à un seul palier => faux ?"""
    D = dict(
        Turquoise = [
            #(nom, fmax (N), dt(s), Longueur roulage(m), commentaire)
            ('BiGolden 4 light', 16797.33, 72-42, 750,'requis à 10 g : 20976.39 N'),
            ('Windtech Ru-Bi 2', 17739.59, 77-48, 750,'requis à 10 g : 22154.09 N'),
            ('Supair Sora 2',    17360.56, 46-15, 750,'requis à 10 g : 21676.96 N'),
                    ],
        Aérotest = [
            ('Hercules',         22840, 56-6 , 1600, '5 pics à 2487 N'),
            ('MacPara Trike 42', 20850, 52-4 , 1600, '5 pics à 2233 N'),
            ('??',               Kg2N(2150), 49-14, 1600, '??'),
                    ]
            )

    
    def toText(operateur, paraps) :
        print ("\n"+(len(operateur)*'=') + '\n' + "%s"%operateur+'\n' + len(operateur)*'=' + '\n')

        separateur = 4*' '+93*'-'
        print(separateur)
        print('    Modèle               : ẟF      :    ẟt   :  Vmax   : remarque                               :')
        print('    '+93*'-')
        for parap in paraps : 
            name, fmax, dt, L, comm = parap#76.39 N) :
            print('    %-20s :'%name, end='')
            parap = list(parap)
            print(" %-7.0f :"%(fmax),end='')
            print(" %-7.1f :"%(dt),end='')
            gamma = L/(dt**2)
            # print(" %-7.2f :"%(gamma),end='')
            Vmax = ms2Kmh(gamma*dt)
            print(" %-7.2f :"%(Vmax),end='')
            print(" %-38s :"%(comm),end='')
            print()
        print('    '+93*'-')
    
    for operateur, valeurs in D.items() : 
        toText(operateur, valeurs)
    print()
    for operateur, valeurs in D.items() :
        print() 
        print() 
        print(toLaTeX(operateur, valeurs))
    
    def toLaTeX(operateur, paraps) :
        hline = "\\hline"
        def dataLine(parap) : 
            name, fmax, dt, L, comm = parap
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
        for parap in paraps :
            # lines.append(hline) 
            lines.append(dataLine(parap))
        lines.append(hline) 
        
        lines.append("\\end{tabular}")
        return '\n'.join(lines)
    # return D

"""
=========
Turquoise
=========

    ---------------------------------------------------------------------------------------------
    Modèle               : ẟF      : ẟF/ẟt   :  𝛾      : remarque                               :
    ---------------------------------------------------------------------------------------------
    BiGolden 4 light     : 16797   : 559.9   : 0.83    : Pas de pic à 10g (requis : 20976.39 N) :
    Windtech Ru-Bi 2     : 17740   : 611.7   : 0.89    : Pas de pic à 10g (requis : 22154.09 N) :
    Supair Sora 2        : 17361   : 560.0   : 0.78    : Pas de pic à 10g (requis : 21676.96 N) :
    ---------------------------------------------------------------------------------------------

========
Aérotest
========

    ---------------------------------------------------------------------------------------------
    Modèle               : ẟF      : ẟF/ẟt   :  𝛾      : remarque                               :
    ---------------------------------------------------------------------------------------------
    Hercules             : 22840   : 456.8   : 0.64    : 5 pics à 2487 N                        :
    MacPara Trike 42     : 20850   : 434.4   : 0.69    : 5 pics à 2233 N                        :
    ??                   : 21092   : 602.6   : 1.31    : ??                                     :
    ---------------------------------------------------------------------------------------------
"""


def run2() :
    """Modélisation charge : et vitesse linéaire à 3 paliers : archi faux !!!"""
    D = dict(
        Turquoise = [900.00,
            #(nom, t0,t1,t2,t3,fmax (N))
            [('BiGolden 4 light', 43,72,77,83,16797.),
            ('Windtech Ru-Bi 2', 48,78,82,92,17739.59, ),
            ('Supair Sora 2',    18,43,51,57,17360.56),]
                    ],
        Aérotest = [1700.00,
            [('Hercules',         6,54,61,63,22840),
            ('MacPara Trike 42', 7,51,56,58,20850),
            ('Stewart',          0,27,34,35,1760*9.81),
            ('Stewart-DGAC',     7,27,34,35,2500*9.81),
            ('Shuttle',         15,49,59,60,2500*9.81),]
                    ]
            )
    def Vmax(t1,t2,t3,L) : 
        return ms2Kmh(2*L/(-t1+t2+t3))
    
    def toText(operateur, L, paraps) :
        print ("\n"+(len(operateur)*'=') + '\n' + "%s"%operateur+'\n' + len(operateur)*'=' + '\n')

        hline = 4*' '+79*'-'
        print(hline)
        print('    Modèle               : t1   : t2   : t3   : Vmax    : Fmax    : 𝛾1    : 𝛾3    :')
        print(hline)
        for parap in paraps : 
            name, t0, t1, t2, t3, fmax = parap
            t1 -= t0
            t2 -= t0
            t3 -= t0
            vmax = Vmax(t1,t2,t3,L)
            gam1 = vmax/t1
            gam3 = -vmax/(t3-t2)
            print('    %-20s :'%name, end='')#---------------------
            print(" %-4.0f :"%(t1),end='')
            print(" %-4.0f :"%(t2),end='')
            print(" %-4.0f :"%(t3),end='')
            print(" %-7.1f :"%(vmax),end='')
            print(" %-7.1f :"%(fmax),end='')
            print(" %-5.1f :"%(gam1),end='')
            print(" %-5.1f :"%(gam3),end='')
            print()
        print(hline)
    
    for operateur, (L,paraps) in D.items() : 
        toText(operateur, L,    paraps)
    print()
    print()
    # for operateur, paraps in D.items() :
    #     print() 
    #     print() 
    #     print(toLaTeX(operateur, paraps))

def run3() :
    """Modélisation charge : linéaire à 3 paliers !!!"""

    """
        Turquoise = [900.00,
            #(nom, t0,t1,t2,t3,fmax (N))
            [('BiGolden 4 light', 43,72,77,83,16797.),
            ('Windtech Ru-Bi 2', 48,78,82,92,17739.59, ),
            ('Supair Sora 2',    18,43,51,57,17360.56),]
                    ],
        Aérotest = [1700.00,
            [('Hercules',         6,54,61,63,22840),
            ('MacPara Trike 42', 7,51,56,58,20850),
            ('Stewart',          0,27,34,35,1760*9.81),
            ('Stewart-DGAC',     7,27,34,35,2500*9.81),
            ('Shuttle',         15,49,59,60,2500*9.81),]
                    ]

    """
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

    def Vmax(t1,t2,t3,L) : 
        return ms2Kmh(3*L/(2*t3+t2-t1))
        # return ms2Kmh(L/(t2-t1/2))/2.5
    
    def toText(operateur, L, paraps) :
        print ("\n"+(len(operateur)*'=') + '\n' + "%s"%operateur+'\n' + len(operateur)*'=' + '\n')

        hline = 4*' '+65*'-'
        print(hline)
        print('    Modèle               : t1   : t2   : t3   : Vmax    : Fmax    : 𝛾1    :')
        print(hline)
        for parap in paraps : 
            name, t0, t1, t2, t3, fmax = parap
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
    
    for operateur, (L,paraps) in D.items() : 
        toText(operateur, L,    paraps)
    print()
    print()
    # for operateur, paraps in D.items() :
    #     print() 
    #     print() 
    #     print(toLaTeX(operateur, paraps))

if __name__ == '__main__' :
    run0()
    # run1()
    # run2()
    # run3()


