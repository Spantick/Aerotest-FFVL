#!/usr/bin/env python
# -*- coding: utf8 -*-
# import Torseur3D
# from utils.geometry3d import ProdVect3d
from utils.debog import debug, mexit, smallStack, souligne
from utils.utilitaires import className
from operator import __matmul__
from numbers import Number
from random import random
from inout.format import fmt
from _ast import If
# from simulation.aerodynamique.postaero import self

## This file in part of Torseur3D
#############################################################################
#############################################################################
##                                                                         ##
##                                   torseur                                 ##
##                                                                         ##
#############################################################################
#############################################################################

## Copyright (C) 2009-2012 Cédrick FAURY

#    Torseur3D is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.

#    Torseur3D is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Torseur3D; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""
torseur.py
Modèle mathématique d'un torseur d'action mécanique
Copyright (C) 2009 Cédrick FAURY

"""
# from widgets import strSc, strRound
import wx
import os
import numpy
from numpy import sqrt, zeros
# from widgets import mathtext_to_wxbitmap

rapportUniteLong = {'m' : 1,
                    'mm': 0.001}
class AlgBase(object):
    """Pour obliger les classes Point, Vecteur, ...
    a implémenter ces méthodes"""
    def __init__(self,*args,**kargs):
        raise NotImplementedError(className(self))

    def copy(self):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def __iadd__(self, x):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def __add__(self,x):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def __isub__(self, x):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def __sub__(self, x):
        raise NotImplementedError(className(self))
        raise NotImplementedError

    def setCoord(self, x):
        raise NotImplementedError(className(self))
        raise NotImplementedError


################################################################################
################################################################################
class Point(AlgBase):
    def __init__(self, x = 0, y = 0, z = 0, nom = "O", unite = 'm'):
        self.x = x
        self.y = y
        self.z = z
        self.nom = nom
        self.unite = unite

    def __eq__(self, other):
        return self.x==other.x and self.y==other.y and self.z==other.z

    def __repr__(self):
        return '%s'%className(self)+"("+str(self.x)+", "+str(self.y)+", "+str(self.z)+")"

    def __str__(self):
        return "("+str(fmt(self.x))+", "+str(fmt(self.y))+", "+str(fmt(self.z))+")"

    def __iadd__(self, V):
        """:returns : Point self, auquel on rajoute V
        >>> self += Vecteur(1,1,1) <=> self.__iadd__(Vecteur(1,1,1))"""
        if not isinstance(V,Vecteur) :
            raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(V)))
        self.x += V.x
        self.y += V.y
        self.z += V.z
        return self

    def __add__(self, V):
        """:returns : nouvelle instance de Point, le point self+V
        >>> self + v <=> self.__add__(v)"""
#         if not isinstance(V,Vecteur) :
#             raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(V)))
        return self.copy().__iadd__(V)

    def __isub__(self,v):
        if not isinstance(V,Vecteur) :
            raise TypeError('On ne sait pas faire %s - %s'%(className(self),className(V)))
        self.x -= V.x
        self.y -= V.y
        self.z -= V.z
        return self

    def __sub__(self,v):
#         if not isinstance(V,Vecteur) :
#             raise TypeError('On ne sait pas faire %s - %s'%(className(self),className(V)))
        return self.copy().__isub__(v)
#         C = self.copy()
#         C.x -= V.x
#         C.y -= V.y
#         C.z -= V.z
#         return C

    def copy(self):
        return Point(self.x, self.y, self.z)

    def setCoord(self, pt):
        """ Appliquer les composantes du Point <pt> a self
        """
        if type(pt) == list and len(pt) == 3:
            self.x = pt[0]
            self.y = pt[1]
            self.z = pt[2]
        elif isinstance(pt, Point):
            self.x = pt.x
            self.y = pt.y
            self.z = pt.z

#     def changeUnite(self, unite):
#         r = rapportUniteLong[self.unite]/rapportUniteLong[unite]
#         self.x = self.x * r
#         self.y = self.y * r
#         self.z = self.z * r
#         self.unite = unite

################################################################################
################################################################################
class Droite():
    def __init__(self, pt, vd, nom = "", unite = 'm'):
        # Point et Vecteur directeur
        self.P = pt
        self.V = vd
        self.nom = nom
        self.unite = unite

    def __repr__(self):
        return "Droite(%s, %s)"%(self.P, self.V)

    def __str__(self):
        return "Δ = {%s + λ.%s, λ∈ℝ}"%(self.P, self.V)

    def intersectionPlan(self, Q, N):
#        print "intersectionPlan"
#        print self
#        print Q, N

        k = self.V@N
        if k == 0.0:
            return None
        PQ = Vecteur(P1 = self.P, P2 = Q)
        PI = self.V*(PQ@N/k)
#        print PI
        I = Point(PI.x+self.P.x, PI.y+self.P.y, PI.z+self.P.z)
        return I
    def __eq__(self, other):
        """2 droites D1(A,u) et D2(B,v) sont les mêmes ssi :
        det(u,v)=0 et det(B-A,u)=0"""
    def setCoord(self, pt):
        """ Appliquer les composantes du Point <pt> à self
        """
        if type(pt) == list and len(pt) == 3:
            self.x = pt[0]
            self.y = pt[1]
            self.z = pt[2]
        elif isinstance(pt, Point):
            self.x = pt.x
            self.y = pt.y
            self.z = pt.z

#     def changeUnite(self, unite):
#         self.P.changeUnite(unite)
#         self.unite = unite

################################################################################
################################################################################
class Vecteur(AlgBase):
    def __init__(self, x = 0, y = 0, z = 0, P1 = None, P2 = None, nom = "V"):
        if P1 is not None and P2 is not None:
            self.x = P2.x - P1.x
            self.y = P2.y - P1.y
            self.z = P2.z - P1.z
        else:
            self.x = x
            self.y = y
            self.z = z
        self.nom = nom

    def copy(self):
        return Vecteur(self.x, self.y, self.z)

    def __repr__(self):
        return className(self)+"("+str(self.x)+", "+str(self.y)+", "+str(self.z)+")"

    def __str__(self):
        return "("+str(fmt(self.x))+", "+str(fmt(self.y))+", "+str(fmt(self.z))+")"
#         return "("+str(self.x)+", "+str(self.y)+", "+str(self.z)+")"

    def __mul__(self, v):
        """ Produit : self*v
            - par un scalaire si v est réel
            - vectoriel si v est Vecteur
        """
#         debug('mul',self=self, v=v)
        if isinstance(v, Vecteur):
            return Vecteur(self.y * v.z - self.z * v.y,
                           self.z * v.x - self.x * v.z,
                           self.x * v.y - self.y * v.x,
                           )
        else:
            return Vecteur(self.x*v, self.y*v, self.z*v)

        return

    def __rmul__(self, v):
        """ Produit : v*self
            - par un scalaire si v est réel
            - vectoriel si v est Vecteur
        """
        debug('rmult')
        if isinstance(v, Vecteur):
            return Vecteur(-self.y * v.z + self.z * v.y,
                           -self.z * v.x + self.x * v.z,
                           -self.x * v.y + self.y * v.x,
                           )
        else:
            return Vecteur(self.x*v, self.y*v, self.z*v)

    def __truediv__(self, k):
        return Vecteur(self.x/k, self.y/k, self.z/k)


    def __matmul__(self, v):
        """
        >>> self@b #<=> self.__matmult__(v), retourne le produit scalaire"""
        return self.x*v.x + self.y*v.y + self.z*v.z

    prod_scal = __matmul__

    def __abs__(self):
        """
        >>> abs(self) <=> self.__abs__()
        """
        return sqrt(self.x**2+self.y**2+self.z**2)

    def __iadd__(self, v):
        """
        Addition
        >>> self += v <=> self.__iadd__(v)
        """
        if isinstance(v, Vecteur):
            self.x += v.x
            self.y += v.y
            self.z += v.z
        else:
            raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(v)))
        return self

    def __add__(self, v):
        """
        Addition
        >>> self + v <=> self.__add__(v)
        """
        if isinstance(v, Vecteur):
            return self.copy().__iadd__(v)
#
#             return Vecteur(self.x + v.x,
#                            self.y + v.y,
#                            self.z + v.z)
        else:
            raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(v)))

    def __sub__(self, v):
        """
        soustraction
        >>> self - v <=> self.__sub__(v)
        """
        if isinstance(v, Vecteur):
            return self.copy().__isub__(v)
#             return Vecteur(self.x - v.x,
#                            self.y - v.y,
#                            self.z - v.z)
        else:
            raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(v)))

    def __isub__(self, v):
        """
        soustraction
        >>> self -= v <=> self.__isub__(v)
        """
        if isinstance(v, Vecteur):
            self.x -= v.x
            self.y -= v.y
            self.z -= v.z
            return self
        else:
            raise TypeError('On ne sait pas faire %s + %s'%(className(self),className(v)))

    def __neg__(self):
        """
        >>> -self <=> self.__neg__()
        """
        return Vecteur(-self.x,-self.y,-self.z)

    def setCoord(self, vect):
        """
        Appliquer les composantes du Vecteur <vect> à self
        """
        if type(vect) == list and len(vect) == 3:
            self.x = vect[0]
            self.y = vect[1]
            self.z = vect[2]
        elif isinstance(vect, Vecteur):
            self.x = vect.x
            self.y = vect.y
            self.z = vect.z

################################################################################
################################################################################
class Torseur(AlgBase):
    def __init__(self, A = Point(), R = Vecteur(), M = Vecteur(), nom = ""):
        """
        Les "éléments de réduction" d'un Torseur T sont A, R et M
            A = origine (espace affine)
            R = résultante
            M = moment en A
        On est certain qu'ils gardent les même id() tout au long de leur vie.
        Donc, excepté à l'initialisation, ils ne peuvent être modifiés que par coordonnée
        >>> T.A = Point(0,0,0)
        ne change pas l'id() de T.A mais modifie son contenu, et celui de T.M
        >>> T.M = Vecteur(1,2,3) #lève une exception
        >>> T.R = Vecteur(5,5,5) # ne modifie pas l'id() de T.R mais son contenu

        L'ev des torseurs est de dim 6. Il suffit donc de 6 parametres pour définir
        un torseur. Ici, on en a 9...
        Oui mais comme un torseur est défini sur un espace affine
        il faut fixer une origine, c'est le point A
        """
        self.A = A
        self.R = R
        self.M = M
        self.nom = nom

    def __repr__(self):
        return '%s(A=%r, R=%r, M=%r)'%(className(self),self.A,self.R,self.M)

    def __str__(self):
        D = self.getAxeCentral()
        inv = self.invariant()
        info = [
            ('Origine (O)', self.A),
            ('Résultante(R)',self.R),
            ('Moment en O (M)',self.M),
            ('Axe central',D),
            ('invariant',inv)
            ]
        if inv==0 : info.append(('glisseur', True))
        if 1: info.append(('ids','self:%d ; A:%d ; R:%d ; M:%d'%(id(self),id(self.A),id(self.R),id(self.M),)))
        return '<%s>\n'%className(self)+'\n'.join(['    %15s : %s'%(k,v) for (k,v) in info])

    def __call__(self, B:Point):
        """Retourne le moment en B (BABAR)"""
        BA = Vecteur(P1=B, P2=self.A)
        return self.M + BA*self.R

    def __iadd__(self, other):
        """self += other <==> self.__iadd__(other)"""
        if not isinstance(other, Torseur):
            raise TypeError('Torseur + %s : opération interdite'%className(other))
#         if other.A != self.A :

#         debug(other,'\n',self)
        cother = other.copy(self.A)#T += T1 ne doit pas modifier ni T1
#         debug(cother,'\n',self)
#             mexit()
        self._R += cother._R
        self.__M += cother.__M
        return self

    def __add__(self, other):
        """self += other <==> self.__iadd__(other)"""
        if not isinstance(other, Torseur):
            raise TypeError('Torseur + %s : opération interdite'%className(other))
        T = self.copy()
#         debug(T,'\n',self)
        T += other
#         debug(T,'\n',self)
#         mexit()
        return T
#         return T.__iadd__(other)
    def __radd__(self, other):
        raise NotImplementedError('--radd__')

    def __isub__(self, other):
        """self -= other <==> self.__isub__(other)"""
        if not isinstance(other, Torseur):
            raise TypeError('Torseur + %s : opération interdite'%className(other))
        if other.A != self.A :
            other = other.copy(self.A)#T += T1 ne doit pas modifier ni T1
        self._R -= other.R
        self.__M -= other.M
        return self

    def __sub__(self, other):
        """self += other <==> self.__iadd__(other)"""
        if not isinstance(other, Torseur):
            raise TypeError('Torseur + %s : opération interdite'%className(other))
        return self.copy(self.A).__isub__(other)

    @property
    def R(self):
        return self._R

    @R.setter
    def R(self, R):
        if hasattr(self, "_R") :
#             raise TypeError('Le résultante R ne peut être modifié que via self.R.setCoord()')
            self._R.setCoord(R)
        else :
            self._R = R

    @property
    def M(self):
        return self.__M

    @M.setter
    def M(self, M):
        if hasattr(self, "__M") :
            raise TypeError('Le moment M d\'un torseur T ne peut pas être modifié')
#             self._M.setCoord(M)
        else :
            self.__M = M

    @property
    def A(self):
        return self._A

    @A.setter
    def A(self,A):
        if hasattr(self, "_A") :
#             raise TypeError('Le point de référence A d\'un torseur T ne peut être modifié que via T.setPtRed(Point)')
            self._setPtRed(A)
        else :
            self._A = A

    def _setPtRed(self, A):
        self.__M += Vecteur(P1 = A, P2 = self._A)*self.R
        self._A.setCoord(A)
#         debug('%r'%self)

    def isCouple(self):
        return self.R@self.R == 0

    def invariant(self, P=None):
        """Pour test et debug : le produit scalaire <R,self(P)> ne doit pas
            dépendre de P"""
        if P is None :
            return self.R@self.M
        else :
            return self.R@self(P)

    def isGlisseur(self, eps=1.0e-10):
        return not self.isCouple() and abs(self.invariant)<eps

    def getAxeCentral(self):
        k = self._R@self._R
        if k == 0.0:
            return None
#         debug('R = %r, M = %r'%(self.R,self.M))
        PI = self._R*self.__M/k
        I = Point(self._A.x+PI.x, self._A.y+PI.y, self._A.z+PI.z)
        return Droite(I, self._R)

#     def setElementsRed(self, tors_R, M = None):
#         if isinstance(tors_R, Torseur):
#             self.R.setCoord(tors_R.R)
#             self.M.setCoord(tors_R.M)
#             self.A.setCoord(tors_R.P)
#         elif isinstance(tors_R, Vecteur):
#             self.R.setCoord(tors_R)
#             self.M.setCoord(M)

    def copy(self, A = None):
        t = Torseur(A=self.A.copy(), R=self.R.copy(), M=self.M.copy())
#         t = Torseur()
#         t.R.setCoord(self.R)
#         t.M.setCoord(self.M)
#         t.A.setCoord(self.A)
        if A is not None and t.A != A:
            t.A = A
        return t
if __name__ == '__main__' :
    from pytictoc import TicToc
    if 0 :
        debug(titre='test Point')
        P = Point(1,0,0);idP=id(P)
        print("P = %r"%P)
        P1 = Point(1,1,1);idP1=id(P1)
        cP = P.copy()
        print('P == P.copy() :',P==P.copy())
        print('id(P) == id(P.copy()) :',id(P)==id(P.copy()))#, 'P==P.copy() :',P==P.copy())
        # P1 = P
        # print('id(P)==id(P1)',id(P)==id(P1))#, 'P==P.copy()',P==P.copy())
        V = Vecteur(2,2,3)
        print('V = %r'%V)
        P += V

        print('P += V : P =', P,'id(P)==idP :',idP==id(P))#, 'P==P.copy()',P==P.copy())
        P = P+V
        print('P = P+V : P =',P,', id(P)==idP',idP==id(P))#, 'P==P.copy()',P==P.copy())
    # mexit()
    if 0:
        debug(titre='test Vecteur')
        I = Vecteur(1,0,0)
        J = Vecteur(0,1,0)
        K = Vecteur(0,0,1)
        debug('produit vectoriel *')
        print('0,0,1 ?=',I*J)
        print('1,0,0 ?=',J*K)
        print('0,1,0 ?=',K*I)
        P = Point(1,2,3)
        print('%r + %r =  %r'%(P,I,P+I))
        try :
            print('%r + %r =  %r'%(P, P, P+P))
        except TypeError as msg :
            print(msg, 'c\'est normal, juste un test')
        try :
            print('%r + %r =  %r'%(I, P, I+P))
        except TypeError as msg :
            print(msg, '. C\'est normal, juste un test')
    #################################################################
    ## VECTEUR TEST, NE PAS MODIFIER
    T  = Torseur(A=Point(0,0,0), R=Vecteur(0,1,0), M=Vecteur(0,1,1))
    ## VECTEUR TEST, NE PAS MODIFIER
    #################################################################
    if 0 :
        debug(titre='méthode T.__call__(Point(x,y,z)) <=> T(Point(x,y,z))')
        print('T=\n',repr(T))
        print('T=\n',str(T))
        print(souligne('Pour B = (x,y,z) on doit trouver T(B) = (z,1,1-x)'))
        for B in (Point(1,0,0),Point(0,0,1),Point(0,1,0),Point(1,1,0),Point(0,1,1),
                  Point(1,0,1), Point(1,1,1), Point(1,2,3), Point(random(),random(),random())):
            print('B = %s'%B,'=> T(B) = %s'%T(B))
        print('T = %r'%T)
    if 0:
        debug(paragraphe="l'invariant est-il invariant ?")
        T1 = Torseur(A=Point  (random(),random(),random()),
                     R=Vecteur(random(),random(),random()),
                     M=Vecteur(random(),random(),random()))
        for i in range(10):
            P = Point(random(),random(),random())
            print("T1.invariant(P)", fmt(T1.invariant(P),10))
    #################################################################
    ## VECTEUR TEST, NE PAS MODIFIER
    T1 = Torseur(Point(1,1,1),Vecteur(1,2,-1), Vecteur(-1,1,-2))
    #################################################################
    ## VECTEUR TEST, NE PAS MODIFIER
    if 0 :
        debug(titre='somme de Torseurs T+=T1')
        print('T  = %r'%(T), 'invariant=%s'%fmt(T.invariant()))
        print('T1 = %r'%(T1),'invariant=%s'%fmt(T1.invariant()))
        T += T1
        print('T += T1')# =>\nT  = %r'%(T))
        print('T  après : %r'%T, 'invariant=%s'%fmt(T.invariant()))# =>\nT  = %r'%(T))
        print('T1 après : %r'%T1, 'invariant=%s'%fmt(T1.invariant()))# =>\nT  = %r'%(T))
        T -= T1
        print('T -= T1')# =>\nT  = %r'%(T))
        print('T  = %r'%(T), 'invariant=%s'%fmt(T.invariant()))
        print('T  après : %r'%T, 'invariant=%s'%fmt(T.invariant()))# =>\nT  = %r'%(T))
        print('T1 après : %r'%T1, 'invariant=%s'%fmt(T1.invariant()))# =>\nT  = %r'%(T))
    # T -= T
    # print('T  après T-=T: %r'%T, 'invariant=%s'%fmt(T.invariant()))# =>\nT  = %r'%(T))
    if 1:
        debug(titre='somme de Torseurs T2=T±T1')
        print('T  : %s'%T)
        print('T1 : %s'%T1)
        T2 = T + T1
    #     print('T  : %s'%T)
    #     print('T1 : %s'%T1)
        print('T2 : %s'%T2)
    #     T += T1
    #     print('T  : %s'%T)
    #     print('T1 : %s'%T1)
    #     print('T2 : %s'%T2)
        T3 = T2-T1
    #     print('T  : %s'%T)
    #     print('T1 : %s'%T1)
    #     print('T2 : %s'%T2)
        print('T3 : %s'%T3)
    #     mexit()
    if 0:
        t = TicToc()
        n=100000
        t.tic()
        for i in range(n):
            Torseur(Point(1,1,1),Vecteur(1,2,-1), Vecteur(-1,1,-2))
        t.toc(restart=True)
        for i in range(n):
            zeros(1000)
        t.toc()
    
