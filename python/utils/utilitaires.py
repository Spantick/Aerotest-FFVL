#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Created on 11 mai 2012
__updated__ = "2019-10-29"
@author: puiseux
utilitaires chaines de caractères,
fonctions geometriques qui marchent en 2d ET 3D
'''
import os
import shutil
import sys,string

import gc#garbage colector
import numpy as np
from numpy import asarray, zeros#, ones, array
from path import Path
import pickle
import datetime
from aperoconfig import RUNS_DIR, DATA_DIR
from utils.debog import debug, rdebug, mexit
from numpy.linalg import norm
from numpy import delete, where, nan, argmax, logical_not, arange
import math
from random import random, choice

def sapero(n=None):
    fn = Path(DATA_DIR)/'art'/'titres.pkl'
    with open(str(fn),'br') as f:
        titres = pickle.load(f)

    if n is None :
        n = choice(range(len(titres)))
    return titres[n]

def _sapero():
    """définition et sauvegarde des titre en ascii-art"""
    sapero = []
    sapero.append(\
"""
   ___,
  /   |
 |    |        _      _      ,_       __
 |    |      |/ \_   |/     /  |     /  \_
  \__/\_/    |__/    |__/      |_/   \__/
            /|
            \|
""")
    sapero.append(\
"""
  __     ___    ___    ___      __
 (  )   (  ,\  (  _)  (  ,)    /  \\
 /__\    ) _/   ) _)   )  \   ( () )
(_)(_)  (_)    (___)  (_)\_)   \__/
""")
    sapero.append(\
"""
    ,.
   / |   ,-. ,-. ,-. ,-.
  /~~|-. | | |-' |   | |
,'   `-' |-' `-' '   `-'
         |
         '
""")
    sapero.append(\
"""
 /\   _   _  _  _
/--\ |_) (- |  (_)
     |
""")
    sapero.append(\
"""
,---.
|---|,---.,---.,---.,---.
|   ||   ||---'|    |   |
`   '|---'`---'`    `---'
     |
""")
    sapero.append(\
"""


               AAA                    PPPPPPPPPPPPPPPPP        EEEEEEEEEEEEEEEEEEEEEE     RRRRRRRRRRRRRRRRR             OOOOOOOOO
              A:::A                   P::::::::::::::::P       E::::::::::::::::::::E     R::::::::::::::::R          OO:::::::::OO
             A:::::A                  P::::::PPPPPP:::::P      E::::::::::::::::::::E     R::::::RRRRRR:::::R       OO:::::::::::::OO
            A:::::::A                 PP:::::P     P:::::P     EE::::::EEEEEEEEE::::E     RR:::::R     R:::::R     O:::::::OOO:::::::O
           A:::::::::A                  P::::P     P:::::P       E:::::E       EEEEEE       R::::R     R:::::R     O::::::O   O::::::O
          A:::::A:::::A                 P::::P     P:::::P       E:::::E                    R::::R     R:::::R     O:::::O     O:::::O
         A:::::A A:::::A                P::::PPPPPP:::::P        E::::::EEEEEEEEEE          R::::RRRRRR:::::R      O:::::O     O:::::O
        A:::::A   A:::::A               P:::::::::::::PP         E:::::::::::::::E          R:::::::::::::RR       O:::::O     O:::::O
       A:::::A     A:::::A              P::::PPPPPPPPP           E:::::::::::::::E          R::::RRRRRR:::::R      O:::::O     O:::::O
      A:::::AAAAAAAAA:::::A             P::::P                   E::::::EEEEEEEEEE          R::::R     R:::::R     O:::::O     O:::::O
     A:::::::::::::::::::::A            P::::P                   E:::::E                    R::::R     R:::::R     O:::::O     O:::::O
    A:::::AAAAAAAAAAAAA:::::A           P::::P                   E:::::E       EEEEEE       R::::R     R:::::R     O::::::O   O::::::O
   A:::::A             A:::::A        PP::::::PP               EE::::::EEEEEEEE:::::E     RR:::::R     R:::::R     O:::::::OOO:::::::O
  A:::::A               A:::::A       P::::::::P               E::::::::::::::::::::E     R::::::R     R:::::R      OO:::::::::::::OO
 A:::::A                 A:::::A      P::::::::P               E::::::::::::::::::::E     R::::::R     R:::::R        OO:::::::::OO
AAAAAAA                   AAAAAAA     PPPPPPPPPP               EEEEEEEEEEEEEEEEEEEEEE     RRRRRRRR     RRRRRRR          OOOOOOOOO







""")
    sapero.append(\
"""

  _____      _____       ______     _____       _____
 (_____)    (_____)     (______)   (_____)     (_____)
(_)___(_)   (_)__(_)    (_)__      (_)__(_)   (_)   (_)
(_______)   (_____)     (____)     (_____)    (_)   (_)
(_)   (_)   (_)         (_)____    ( ) ( )    (_)___(_)
(_)   (_)   (_)         (______)   (_)  (_)    (_____)


""")
    sapero.append(\
"""
  AAA    PPPPPP   EEEEEEE  RRRRRR    OOOOO
 AAAAA   PP   PP  EE       RR   RR  OO   OO
AA   AA  PPPPPP   EEEEE    RRRRRR   OO   OO
AAAAAAA  PP       EE       RR  RR   OO   OO
AA   AA  PP       EEEEEEE  RR   RR   OOOO0
          """)
    fn = Path(DATA_DIR)/'art'/'titres.pkl'

    with open(fn,'wb') as f :
        pickle.dump(sapero, f)

def loadDict(self, filename, type_=''):
    with open(filename,'r'+type_) as f :
        D = eval(f.read())
    for k,v in D.items():
        setattr(self, k, v)
    return D

def suppressionPoints(P, S, C=[], debog=False):
    """
    obsolete ??
    Suppression points de P referencés par S, et mise à jour des connexions C
    :param P: ndarray((nbp,dim),dtype=float ou int) tableau de points à nettoyer.
    :param S: ndarray(ns, dtype=int), les n° de points à supprimer
        max(S)<nbp.
    :param C: ndarray(nc,dtype=int) tableau de n° de points à mettre à jour
        max(C)<nbp.
        - si C∩S = ∅, on ne supprime que des points non référencés dans C, RAS
        - si C∩S ≠ ∅, pour i∈C∩S, dans C on remplace i par -9999
    :return P, C: le tableau de points modifié et les connexions mises à jour.
    """
    K = sorted(list(set(range(len(P))) - set(S)))#Keep
    debug(K=K)
    P = P[K]
    S = list(S)
    S.sort(reverse=True)
    for ks in S :
        if ks in K :
            i = K.index(ks)
            print('suppression de %d impossible'%ks)
        else :
            C[where(C>ks)] -= 1
#     print("C=", C)
    return P, C

def pointsUniques(A, eps=1.0e-3):
    """
    #http://scipy-user.10969.n7.nabble.com/remove-duplicate-points-td4636.html
    :param A: ndarray((nbp,dim),dtype=float ou int) tableau de points contenant
        possiblement des points doubles
        Il vaut mieux que dim<<nbp, pour tri.
    :returns j_unique : le tableau des n° de points uniques
    :returns j_sorted : la permutation telle que A[j_sorted] est trié
    :returns unique_mask : un masque (tableau booleens) t.q
        j_unique == j_sorted[unique_mask]
    NB la liste des n° de points doubles donc j_sorted[logical_not(unique_mask)]
    """
    j_sorted = np.lexsort(A.T)
    v = None
    for q in A.T:
        q = q[j_sorted] # q = le tableau X des x_i, idem Y, idem Z
#         print('q=',q.tolist())
        w = (np.abs(q[1:] - q[:-1])>eps)# liste de True ou False : w[i] vrai si q[i] == q[i+1]
#         print('w=',w.tolist())
        if v is None:
            v = w
        else:
            v |= w #v[i] = v[i] ou w[i]
#         print('v=',v.tolist())
    unique_mask = np.hstack([True, v])
    j_unique = j_sorted[unique_mask]
#     j_double = j_sorted[logical_not(unique_mask)]
    return j_sorted, unique_mask, j_unique

# def supprimerPointsOld(A, S, C=asarray([],dtype=int), replace=True):
#     """
#     suppression dans A des points numeros j∈S et mise à jour de C.
#     :param replace:bool
#         - Si replace = True, la mise à jour de C consiste à remplacer le point
#             n° j∈S par le point n°j-1
#         - Si replace = False, la mise à jour de C consiste à supprimer toutes les
#             connexions de C qui se réfèrent aux points n° j∈S
#     :param A: ndarray((nbp,dim),dtype=float ou int) tableau de points a nettoyer.
#         Il vaut mieux que dim<<nbp, pour tri.
#     :param S: ndarray(n,dtype=int) numéros des points à supprimer max(S)∈[0, nbp[
#         - le point n° j∈S (à supprimer) est remplacé par le point n°j-1
#         - toutes les connexions faisant référence à un j∈S sont supprimées
#     :param C: ndarray((nc,n),dtype=int) tableau de connexions entre les points
#         de A, avec max(C)<nbp et nc=1,2,3 ou plus.
#
#     :return A, J: le tableau de points modifié et les connexions mises à jour.
#
#     """
# #     U.sort()
# #     S.sort()
#
#     As = delete(A, S, axis=0)
# #     Au = A[U]
#     print('A =',A.tolist())
# #     print('U =',U)
# #     print('Au =',Au.tolist())
#     print('S =',S)
#     print('As =',As.tolist())
#     A = As
#     if replace is True :
#         """ne marche pas"""
#         J = C.view()
#         J.shape = -1#on fait de C un tableau a une seule dimension
#         debug(J=J, S=S)
#         S.sort()
#         for js in S[::-1] :
#             print ('js=%d'%js, end=' ; ')
#             w = where(J>=js)
#             print ('where(J>=js) = %-20s'%w[0], end=' ; ')
#             J[w] -= 1
#             print (' => J=',J)
#         J.shape = C.shape
#         return A, J
#     else:#replace=False,
# #         raise NotImplementedError
#         #On supprime les connexions CONTENANT les points j∈S
#
#         tokeep = []
#         for k, c in enumerate(C) :
#             keep = True
#             for i in c :
#                 if i in S :
#                     keep = False
#                     continue
#             if keep : tokeep.append(k)
# #         print('C=',C.tolist(), 'tokeep=',tokeep)
#         C = C[tokeep]
#         #Décalage des n° j∈C.flat() tq j>js, js∈S
# #         print('C=',C.tolist())
#         for js in S :
#             for c in C : c[where(c>js)]-=1
#         return A, C

def supprimerPoints(A, S, C=asarray([],dtype=int)):
    """
    ATTENTION, supprime également toutes les connexions de C qui contiennent un j∈S
    suppression dans A des points numeros j∈S et des connexions contenant j.
    :param A: ndarray((nbp,dim),dtype=float ou int) tableau de points a nettoyer.
        Il vaut mieux que dim<<nbp, pour tri.
    :param S: ndarray(n,dtype=int) numéros des points à supprimer max(S)∈[0, nbp[
        toutes les connexions faisant référence à un j∈S sont supprimées
    :param C: ndarray((n,coord),dtype=int) tableau de connexions entre les points
        de A, avec max(C)<nbp et coord=1,2,3 ou plus.

    :return A, C: le tableau de points modifié et les connexions mises à jour.

    """
    #On supprime les points j∈S
    As = delete(A, S, axis=0)
#     print('A =',A.tolist())
#     print('S =',S)
#     print('As =',As.tolist())
    A = As
    S = sorted(list(set(S)))
    #On supprime les connexions CONTENANT les points j∈S
    tokeep = []
    for k, c in enumerate(C) :
#         print ('(k,c) =', (k,c.tolist()))
        keep = True
        for i in c :
            if i in S :
                keep = False
                continue
#         print('keep=',keep)
        if keep : tokeep.append(k)
#     print('C=',C.tolist(), 'tokeep=',tokeep)
    C = C[tokeep]
    #Décalage des n° j∈C.flat() tq j>js, js∈S
    for js in S[::-1] :
        for c in C : c[where(c>js)] -= 1
    return A, C

# def supprimerPointsDoublesNeMarchePas(A, C=asarray([], dtype=int), eps=1.0e-3):
#     #calcul des points doubles (Pd)
#     _, Pd = pointsUniques(A, eps)
#     """C'est supprimerPoints qui ne marche pas"""
#     A, C = supprimerPoints(A, Pd, C, replace=True)
#     debug('Apres_supprimerPoints',A=A, C=C.tolist(), )
#     #suppression des connexions doubles [i,j], [i,j]
#     C.sort()
#     _, Cd = pointsUniques(C,eps=0)
#     if len(Cd)>0 :
#         C, _ = supprimerPoints(C, Cd, replace=False)
#     #suppression des connexions de type [i,i]
#     C = supprimerAutoConnexions(C)
# #     tokeep = []
# #     for k,c in enumerate(C) :
# #         if not np.all(c==c[0]) :
# #             tokeep.append(k)
# #     debug(tokeep=tokeep)
# #     C = C[tokeep]
#     return A,C

def supprimerAutoConnexions(C):
    """
    :param C: ndarray((nc,coord),dtype=int), un tableau de connexions i.e. une
        liste de nc cellules [i_1,i_2,...,i_coord] coord est la coordinance,
        i.e le nb de points des cellules (line=>2, triangle=>3 etc...)
    une auto-connexion est une cellule [i,i,...,i]. On les supprime
    """
    tokeep = []
    for k,c in enumerate(C) :
        if not np.all(c==c[0]) :
            tokeep.append(k)
#     debug(tokeep=tokeep)
    return C[tokeep]

def supprimerPointsDoubles(A, C=asarray([], dtype=int), eps=1.0e-3):
    """
    suppression points doubles de A et mise à jour de C, le cas échéant
    :param A: ndarray((nbp,dim),dtype=float ou int) tableau de points a nettoyer.
        Il vaut mieux que dim<<nbp, pour tri.
    :param C: ndarray((nc,coord),dtype=int) tableau de connexions entre les points
        de A, avec max(C.flat)<nbp. et coord=1,2,3,... est la coordinance,
        i.e. le nb de (n° de) points de chaque connexion.
    :returns A1, C1: le tableau de points et les connexions mises à jour.
    :Attention : A et C sont modifiés, si on veut les conserver,
        il faut en faire une copie avant l'appel à supprimerPointsDoubles()
    """
    j_sorted, unique_mask, j_unique = pointsUniques(A, eps)
#     debug(j_sorted=j_sorted, A=A.shape)
    if len(C)==0 :
        #Cas ou il n'y a pas de tableau de connexions C
        j_unique.sort()#On veut les points uniques, mais dans l'ordre initial
        return A[j_unique], []
    #Cas ou il y a un tableau de connexions C à traiter
    #les n° des points doubles de A, triés
    avirer = sorted(j_sorted[logical_not(unique_mask)])
#     debug(avirer=avirer)
    J = C.view()
    J.shape = -1#on fait de C un tableau a une seule dimension
#     avirer = sorted([i for i,v in zip(j_sorted, unique_mask) if not v])
#     newJ = np.arange(0, 1+np.max(J))#il peut y avoir des n° de points>max(J)
    newJ = np.arange(0, 1+max(j_sorted))
    #2. Nouveaux pointeurs sur points
    # Si je vire le point j, je le remplace par le point newJ[j]=ju
    # newJ ne contient que les (numéros de) points à garder.
    for k, (v, j) in enumerate(zip(unique_mask, j_sorted)) :
        if v :
            ju = j#v=True => on garde le point j
        else :
            newJ[j] = ju #le point j est remplacé par ju
    #3. Je supprime effectivement les n° de points avirer, les points agarder ont
    # un nouveau numéro que je met dans newJ
    # si je supprime le point j, tous les points k>j reculent d'une place k-=1
    for j in avirer[::-1] :
        #il faut partir du numero le +élevé
        #sinon les n° avirer doivent être aussi décalés de -1
        newJ[newJ>j] -= 1
    for k, j in enumerate(J) :
        J[k] = newJ[j]
    j_unique.sort()#On les veut dans l'ordre initial
    C,_ = supprimerPointsDoubles(C, eps=0)
#     C = supprimerAutoConnexions(C)
    if False : C.sort(axis=1)#transforme les connexions (i,j) en (min(i,j), max(i,j))
    return A[j_unique], C

def supprimerPointsDoubles1(A, eps=1.0e-3, newnum=True):
    """
    suppression points doubles de A
    Dans A soit un point multiple que l'on retrouve en position i1<i2<i3...
    Alors
    - A[i2], A[i3],... sont supprimés.
    - les points sont renumérotés la nouvelle numérotation est retournée dans newJ
        - tous les indices i2, i3,... sont transformés en i1
        - les indices sont mis à jour pour tenir compte du chgt de
            position des points d'indice ≥ i2
            (si on supprime le point i, les n° de points≥i sont décrémentés de 1)

    :param A: ndarray((nbp,dim),dtype=float ou int) tableau de points a nettoyer.
        Il vaut mieux que dim<<nbp, pour tri.
    :param eps: float, deux points sont considérés comme identiques si leur distance est <eps
    :param newnum : bool, True si l'on veut la nouvelle numérotation, False sinon
    :returns A1: le tableau de points mis à jour.
    :returns newJ: ndarray(nbp, dtype=int) nouvelle numérotation
        si J est un(e liste de) n° de point(s) de l'ancienne numérotation, alors
        newJ[J] est le (la liste des) n° de points dans la nouvelle numérotation
    """
    #1. calcul des points doubles
    j_sorted, unique_mask, j_unique = pointsUniques(A, eps)
    if not newnum :
        j_unique.sort()#On les veut dans l'ordre initial
        return A[j_unique]
    #
    #2. Nouvelle numérotation
    #les n° des points doubles de A, triés
    avirer = sorted(j_sorted[logical_not(unique_mask)])
    # Si je vire un point n°j, je le remplace par le point n° newJ[j]=ju
    newJ = np.arange(0, 1+max(j_sorted))#nouvelle numérotation
    for (v, j) in zip(unique_mask, j_sorted) :
        if v :
            ju = j#v=True => on garde le point j
        else :
            newJ[j] = ju #le point j est remplacé par ju
    #si je supprime le point j, tous les points k>j reculent d'une place k-=1
    for j in avirer[::-1] :
        #il faut partir du numero le +élevé
        #sinon les n° avirer doivent être aussi décalés de -1
        newJ[newJ>j] -= 1
    j_unique.sort()#On les veut dans l'ordre initial
    return A[j_unique], newJ

def className(obj):
    try : return obj.__class__.__name__
    except : return 'Inconnu'

def locate(t, T, eps = 1.0e-5):
    """(2D et 3D)
    Localise t dans T=ndarray(n).
    On doit avoir T[0]<=t<=T[-1] et T croissant, len(T)>=2
    retourne (True_ou_False, index),
    - (True, k) si t=T[k] à eps près (t est une des valeurs T[k])
    - (False, k) si T[k] < t < T[k+1] (t n'est pas un T[k])
    """
    g = where(T<=t+eps)[0]
    d = where(T>=t-eps)[0]
    return (True, g[-1]) if g[-1]==d[0] else (False, g[-1])

def computeCordeAndNBA(points):
    """
    :param points:ndarray((n,dim), dtype=float)
    retourne la distance et le numéro du point le plus éloigné de BF=points[0]
    Dés que la distance décroit, on y est.
    """
    p1 = points - points[0]
    n = norm(p1,axis=1)
    i = argmax(n)
    return n[i],i
#     corde, nba = -1.0, np.nan
#     bf = points[0]
#     for k, point in enumerate(points) :
#         d = dist2(bf, point)
#         if d > corde :
#             corde, nba = d, k
#     return math.sqrt(corde), nba

def computeCordeAndNBAOld(points):
    """
    :param points:ndarray((n,dim), dtype=float)
    retourne la distance et le numéro du point le plus éloigné de BF=points[0]
    Dés que la distance décroit, on y est.
    """
    corde, nba = -1.0, np.nan
    bf = points[0]
    for k, point in enumerate(points) :
        d = dist2(bf, point)
        if d > corde :
            corde, nba = d, k
    return math.sqrt(corde), nba

def dist2(p1,p2):
    """retourne le carré de la distance de p1 à p2 en norme n=2"""
    try :
        return sum(v**2 for v in p2-p1)
    except TypeError : # Si p1 ou p2=None
        return nan

def dist(p1,p2):
    """retourne la distance de p1 à p2 en norme n=2"""
    try :
        return math.sqrt(sum(v**2 for v in p2-p1))
    except TypeError : # Si p1 ou p2=None
        return nan

def moyenneMobileClosed(points,molecule=[1.,2.,1.],nnodes=None):
    """
    Lissage moyenne mobile pour les polygones fermés self[-1]==self[0]
    nnodes est la liste des numeros de points à lisser. On suppose que 0 n'y est pas
    """
    molecule=np.asarray(molecule)#np.asarray([1., 2., 1.])/4.
    molecule/=np.sum(np.absolute(molecule))
    new=points
    if nnodes is None :
        nnodes=list(range(len(points)))
    else :
        nnodes.sort()
        if nnodes[0]==0 and nnodes[-1]==len(points) :
            debug('non implémenté')
            return new
        else : pass
    old=new[:-1]
    n=len(old)
#         deb, fin = 0, n
    deb,fin=nnodes[0],nnodes[-1]
    for k in range(deb,fin) :
        pm=old[k-1]
        p=old[k]
        if k==n-1 :
            pp=old[0]
        else:
            pp=old[k+1]
        new[k]=molecule[0]*pm+molecule[1]*p+molecule[2]*pp
    new[-1]=new[0]
    return new

def moyenneMobile(X,n=1):#,molecule=[1.,2.,1.]):
    """
    Lissage moyenne mobile pour fonction simple i -> X[i]
    Y[i] = (X[i-1]+2X[i]+X[i+1])/4, 0<i<n
    Y[0] = (       2X[0]+X[1])/3, i=0
    Y[n] = (X[n-1]+2X[n]       )/3, i=n
    """
    Y = X.copy()
    for k in range(1, len(X)-1):
        Y[k] = 0.25*(X[k-1] + 2*X[k] + X[k+1])
    Y[0]  = (        2*X[0] + X[1])/3.0
    Y[-1] = (X[-2] + 2*X[-1]      )/3.0
    if n==1 :
        return Y
    else :
        return moyenneMobile(Y, n-1)

def moyenneMobile1(X, n=1):#,molecule=[1.,2.,1.]):
    """
    n lissages par moyenne mobile pour fonction simple i -> X[i]
    Y[i] = (X[i-1]+2X[i]+X[i+1])/4, 0<i<n
    Y[0] = X[0]                   , i=0
    Y[n] = X[n]                   , i=n
    """
    Y = X.copy()
    for k in range(1, len(X)-1):
        Y[k] = 0.25*(X[k-1] + 2*X[k] + X[k+1])
    Y[0]  = X[0]
    Y[-1] = X[-1]
    if n==1 :
        return Y
    else :
        return moyenneMobile1(Y, n-1)

def maintenant():
    """La date et l'heure formatées "human readable\""""
    return str(datetime.datetime.now().strftime("%A, %d. %B %Y %I:%M%p"))

def toDict(cles,valeurs):
    """
    - cles = "cle1 cle2 cle3 ..."
    - valeurs = "val1 val2 val3...", les valeurs sont des entiers ou des reels
    retourne un dictionnaire cle,valeurs
    """
    d={}
    for key,value in zip(cles.split(),valeurs.split()) :
        try: w=int(value)
        except ValueError : w=float(value)
        d[key.lower()]=w
    return d

def findAll(tag,lines,first=0,last=-1):
    """
    Retourne une liste triée des numeros de ligne de
    TOUTES les occurences de tag dans lines[first:last]
    """
#     debug(len(lines),tag=repr(tag), first_last=(first,last))
    n=first-1
    n0=findRel(tag,lines,n+1,last)
    if n0 is None : return ()
    N=[]
    while n0 is not None and n0>=0 :
        n=n0+n+1
        N.append(n)
        n0=findRel(tag,lines,n+1,last)

    return tuple(N)

def findAllLines(tag,lines,first=0,last=-1):
    """
    Comme findAll(), mais il faut que la ligne complete (nettoyée) soit égale à tag, pas seulement une partie de la ligne.
    Par exemple : line = 'TOTO_EST_CONTENT' ne match pas avec tag='TOTO'
    """
    liste=findAll(tag,lines,first,last)
    newlist=[]
    for n in liste :
        if lines[n].strip()==tag :
            newlist.append(n)
    return tuple(newlist)

def findRel(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i-first de la premiere occurence trouvee,
          c'est à dire le numéro de ligne RELATIF : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    found=find(tag,lines,first,last)
    if found is None : return None
    else : return found-first

def find0(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i-first de la premiere occurence trouvee,
          c'est à dire le numéro de ligne RELATIF : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    found=findRel(tag,lines,first=first,last=last)
    if found is not None : return  found+first
    else : return None

def find(tag,lines,first=0,last=-1):
    """
    Cherche 'tag' dans lines[first:last] et retourne
        * le numero de ligne i de la premiere occurence trouvee,
          c'est à dire le numéro de ligne ABSOLU : dans self.lines[0:]
        * None si le tag n'est pas trouvé
    """
    if last<0 : last=len(lines)+last
    if first is None : first=0
    elif first<0 : first=len(lines)+first
    if not tag : return None
    i=first-1
    while i<last:
        i+=1
        try :
            if lines[i].find(tag)>=0 :
                return i
        except IndexError :
            return None
    return None

def rreplace(chaine,sub,replace):
    '''
    Remplace récursivement dans chaine, toutes les occurences de sub  par replace
    >>> rreplace("aaaaab","aa","a")
    ... "ab"
    '''
    while 1:
        oc=chaine
        chaine=oc.replace(sub,replace)
        if oc==chaine : break
    return chaine

def find_names(obj):
    """Renvoit les noms de variable de l'objet 'obj' passé en paramètre"""
    frame=sys._getframe()
    for frame in iter(lambda: frame.f_back,None):
        frame.f_locals
    result=[]
    for referrer in gc.get_referrers(obj):
        if isinstance(referrer,dict):
            for k,v in referrer.items():
                if v is obj:
                    result.append(k)
    return result

def load(filinname):
    try :
        filin=open(filinname,'r')
        return pickle.load(filin)
    except IOError :
        print("Impossible d'ouvrir le fichier dump %s, pas de lecture."%filinname, file=sys.stderr)

def goodName(name):
    allowed=string.ascii_letters+string.digits+'_-+=#'
    for letter in name :
        if letter not in allowed :
            return False
    return True

def diff(A):
    """
    J'en ai marre de taper tout le temps 'dA = A[1:]-A[:-1]'
    :param A: un tableau de points np.ndarray de shape (n, dim)
    :return dA: ndarray((n-1,dim)), les différences premières de A:
        dA[i] = A[i+1]-A[i]
    """
    return A[1:]-A[:-1]

def absCurv(X, normalise=False):
    """Abscisse curviligne des points de X (2d/3d)
    :param X: ndarray((n,dim))"""
    l = len(X)
    if l==0 : return []
    elif l==1 : return [0]

    T = zeros(l)
    for k in range(l-1) :
        T[k+1] = T[k] + norm(X[k]-X[k+1])
    if normalise and T[-1] != 0.0 :
        T /= T[-1]
        #sinon T=[0,0,...] On peut très bien entrer ici avec 3 points = 0
    return T

# def pointsDoubles(points, eps=1.0e-10):
#     """retourne la liste des numéros de points doubles CONSECUTIFS à eps pres.
#     points est de shape(n,2), n>=2"""
# #     if not isinstance(points, (np.ndarray,list,tuple)) :
# #         raise NotImplementedError('points doit etre de type numpy.ndarray')
#     if len(points)<=1 : return []
#     dac = diff(absCurv(points))
# #     dac = ac[1:]-ac[:-1]
#     return [k for (k,d) in enumerate(dac) if abs(d)<=eps]

def pointsDoublesNew(points, eps=1.0e-10):
    """Calcule uniquement les (n° de) points doubles consécutifs
    Cf fonction pointsUniques()"""
    return where(norm(diff(points), axis=1)<=eps)[0].tolist()

def eliminerPointsDoublesConsecutifs(points, eps=0.0, vires=False):
    """a refaire en utilisant pointsUniques()"""
    avirer = pointsDoublesNew(points, eps)
#     debug(len_avirer=len(avirer))
    if avirer : points = delete(points, avirer, 0)
#     rdebug(points)
    return (points, avirer) if vires else points
#         msg1 = 'Les points de controle %s sont des doublons consecutifs.'%avirer
#         msg1 += 'Suppression des doublons'
#         rdebug(msg1,'')
#                 rdebug(avirer=avirer, delta_abscurv=dac.tolist(), controle=sum(dac))
#         for k in reversed(avirer) :
#             self.removePoint(k, update=False)

def safeRemoveTree(rep:Path) :
    if rep.isdir() :
        contents = sorted(os.listdir(rep))
    elif rep.isfile() :
        contents = [rep]
    if contents :
        debug (paragraphe="Suppression de %s. Contenu :"%rep.name)
        for c in contents :
            c = Path(c)
            if c.isdir() : print ('%30s (D)'%c.name)
            elif c.isfile() : print ('%30s (F)'%c.name)
            else :  print ('%30s (?)'%c.name)
    else :
        debug(paragraphe="%s est vide"%rep.name)
    ans = 'w'
    while ans not in ('y','n','') :
        ans = input('\nOK pour suppression ? y/n (defaut : n)').lower()
    if ans == 'y' :
        rdebug("JE REMOVE %s"%rep)
        return shutil.rmtree(rep)
    else :
        rdebug("CANCELLED REMOVE %s"%rep)
        return None

def dinput(prompt, defaut='',type_=str) :
    """input avec valeur par defaut et type de résultat
    :param defaut: valeur par defaut
    :param type_: class
    """
    if defaut == '':
        return type_(input(prompt))
    else :
        value = input(prompt+'[%s]'%defaut)
        return type_(defaut) if value == '' else  type_(value)

def ordered(objet) :
    """
    :type objet: un iterable (list, tuple, ) ou un dict
    :returns un tuple trié
    """
    if isinstance(objet, (dict, )) :
        return sorted([(key, objet[key]) for key in objet])
    else :
        return sorted(objet)

if __name__=="__main__":
    _sapero()
    for a in range(10) :
        print(sapero(a))
    mexit()
    from tests.utils.testsutilitaires import TestUtils
    t = TestUtils()
    if 0 : t.testSafeRemoveTree()
    if 1 : t.testsDivers()
    if 0 : t.testPointsDoubles()
    if 0 : t.testSuppressionPointsDoubles()
    if 0 : t.testSuppressionPointsDoubles1()
    if 0 : t.testSuppressionPoints()


