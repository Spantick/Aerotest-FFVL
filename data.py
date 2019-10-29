#!/usr/bin/python
#coding:utf-8
class A(object):
    def __init__():
        pass



def ruin():
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

    
    def toText(boite, bis) :
        print ("\n"+(len(boite)*'=') + '\n' + "%s"%boite+'\n' + len(boite)*'=' + '\n')

        separateur = 4*' '+93*'-'
        print(separateur)
        print('    Modèle               : ẟF      :    ẟt   :  Vmax   : remarque                               :')
        print('    '+93*'-')
        for bi in bis : 
            name, fmax, dt, L, comm = bi#76.39 N) :
            print('    %-20s :'%name, end='')
            bi = list(bi)
            print(" %-7.0f :"%(fmax),end='')
            print(" %-7.1f :"%(dt),end='')
            gamma = L/(dt**2)
            # print(" %-7.2f :"%(gamma),end='')
            Vmax = ms2Kmh(gamma*dt)
            print(" %-7.2f :"%(Vmax),end='')
            print(" %-38s :"%(comm),end='')
            print()
        print('    '+93*'-')
    
    for boite, valeurs in D.items() : 
        toText(boite, valeurs)
    print()
    for boite, valeurs in D.items() :
        print() 
        print() 
        print(toLaTeX(boite, valeurs))
    # return D

def Kg2N(x) : return x*9.81
def ms2Kmh(x) : return x*3.6

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

if __name__ == '__main__' :
    pass
    # ruin()


