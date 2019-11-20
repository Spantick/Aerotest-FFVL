#!/usr/local/bin/python
#-*-coding: utf-8 -*-

'''
Created on 11 mai 2012

@author: puiseux
'''
import sys, math, os
from utils.debog import debug, stack, mexit
from plotly.offline import plot
import plotly.graph_objs as go
from shapely import speedups
from scipy.interpolate.fitpack2 import InterpolatedUnivariateSpline
from utils.utilitaires import (absCurv, className)
from utils.debog import rdebug, rstack
from preferences import GraphicsPrefs
from array import array
from numpy.linalg.linalg import norm, det
from scipy import random
from pytictoc import TicToc
speedups.enable()
import numpy as np
import numpy.linalg as LA
from numpy import (sin, cos, ones, linspace, zeros, average, matrix, pi,
                   ndarray,asarray, cross)
from aperoconfig import RUNS_DIR#, DATA_DIR
from scipy.interpolate import (SmoothBivariateSpline,)
from scipy.optimize import  minimize_scalar
from vtk import vtkXMLPolyDataWriter
from vtk_nervures.vtkUtilitaires import toVtkPolyData
from path import Path
from utils.utilitaires import dist
OUTPUT = RUNS_DIR/'trash'

def computeSurfaces(points, closed):
    """
    Calcul des surfaces des panels définis par points
    :param points:ndarray((nbn,nbpn,3),float) des points 3d d'un maillage régulier
    Les panels sont délimités par points[i,j] points[i+1,j], points[i+1,j+1], p[i,j+1]
    :param closed: bool
        True si le bord de fuite est en double (points[k,-1]==points[k,0])
            c'est le cas d'un maillage aerodynamique
        False si non, c'est le cas du maillage elastique.
            Dans ce cas, il faut rajouter un panel au bout de chaque caisson
    """
    nbn, nbpn,_ = points.shape #ici les nervures sont fermées en BF... ou pas
    debug(shape_points=(nbn,nbpn))
    nbcais = nbn-1
    points = points.view().reshape((-1,3))
    nbpanels = nbcais*(nbpn-1) if closed else nbcais*nbpn
    debug(points=points.shape, nbpanels=nbpanels)
    surf = zeros(nbpanels)
#         debug(nbc=nbc, nbpn=nbpn, points=shape, surf=surf.shape)
#         debug(self_nbp=nbp, self_nbn=nbn, self_nbpn=nbpn)
#         exit()
    for ip, panel in enumerate(panels((nbn, nbpn), closed)):
#         if ip>=nbpanels-1 : debug(ip=ip, panel=panel)
        sommets = [points[k] for k in panel]
        surf[ip] = airePanel(*sommets)
#             debug(p=p, panel=panel, surf=surf[p])
#         celldata['surfaces'] = surf
#         mesh.celldata = celldata
    return surf


def doublons(points,voisinages=[],eps=1.0e-8,sym=True):
    """
    Verifier s'il y a des points doubles dans points.
    Paramètres :
    ----------
    :param points : ndarray(nbp,3) nbp points 3d. nbp peut être une
        shape (nl, nc, 3) et dans ce cas nbp = nl*nc
    :param voisinages: ndarray(shape=(np,1+nv), dtype=int) np numéros de points,
        faisant référence au tableau 'points', chacun des np points à
        vérifier possède nv voisins. voisinage est donc de la forme
        [
         [k_0,  v_0_1,  v_0_2, ...,  v_0_nv],
         [k_1,  v_1_1,  v_1_2, ...,  v_1_nv],
         ...
         [k_np, v_np_1, v_np_2, ..., v_np_nv]
        ]
        les valeurs k_i, v_i_j sont des entiers entre 0 et nbp-1
    :param eps : float. deux points p,q sont considérée comme confondus si dist(p,q)<eps
    :param sym : bool. True si i voisin de j <==> j voisin de i. Divise par deux
        le travail.

    Retourne :
    --------
    :return : une liste de doublons (d_0,dd_0), ... (d_n, dd_n)
    """
    def cond(v,k):
        if sym : return v>k
        else : return v!=k

    qoints=points.view()
    qoints.shape=(qoints.size//3,3)
    nbp=qoints.shape[0]
    if voisinages!=[] :
        voisinages=asarray(voisinages)
#        debog(whoami(), voisinages.shape)
        nbp=voisinages.shape[-1]
        voisinages.shape=(-1,nbp)
    else :
        if sym :
            voisinages=asarray([[i]+list(range(i,nbp)) for i in range(nbp)])
        else :
            voisinages=asarray([[i]+list(range(nbp)) for i in range(nbp)])

    doublons=[]
    for voisinage in voisinages:
        k,voisins=voisinage[0],voisinage[1:]
        point=qoints[k]
        for voisin in voisins :
            if cond(voisin,k) :
                if np.linalg.norm(point-qoints[voisin])<eps :
                    doublons.append((k,voisin))
    return doublons

def quadToPoint(shape, celldata):
    '''
    Transformer 'celldata', de type CellData(QUAD), en PointData,
        en répartissant celldata[k] sur les 4 sommets du panel k
    :param shape : (nbn, nbpn, dim)= nb nervures, nb points par nervure,dim=1 ou 3
        (doit provenir de l'elasticité, bf n'est pas en double, nbpn=85)
    :param celldata : np.ndarray((N,dim)) ou CellData(), avec N=(nbn-1)*nbpn
    ATTENTION:les numeros de points dans les panels sont des numéros de la geométrie aero,
    ie dans le tableau de points aero, de shape (nbn_aero, nbpn_aero,3)
    ou nbpn_aero est le nb de points par nervure de la géométrie aero, (ou le bf est en double)
    nbpn_aero = nbpn_elasto+1

    '''
    vshape = celldata.shape
    if len(shape)==2 : nbn, nbpn_e, _ = shape, 1
    elif len(shape)==3: nbn, nbpn_e, _ = shape #nb nerv, nb pt par nerv, 3, maillage elastique
#     nbpn_a = nbpn_e+1
    N = (nbn-1)*nbpn_e#nb panels QUAD
#     rdebug(vshape=vshape, shape_e=shape, nbn=nbn, nbpn_e=nbpn_e, nbpanels=N)
    if vshape[0] != N :
        msg = "Erreur de dimension : celldata doit etre de shape (%d) ou (%d,%d) ou (%d,%d). Or il est de shape %s"%(N,N,1,N,3,str(vshape))
        raise ValueError(msg)
    pointdata = zeros(shape).reshape((nbn*nbpn_e,-1))
#     debug(pointdata=pointdata.shape, celldata=celldata.shape, dim=dim)
    for p, panel in enumerate(panels(shape, closed=False)):
#         if p>=255 : continue
        v = 0.25*celldata[p]
#         if p <= 4000 :
#             debug('panel %d'%p,panel=panel, dF=v)
        for k_e in panel : #ka=indice aero ke=indice elasto
#             k_a = e2a(k_e)
            pointdata[k_e] += v
#             if p <= 4 :
#                 print 'pointdata[%d]='%k_e,pointdata[k_e]
    return pointdata

# def splineInterpolation(points, methode='c cubic', tension=5, degre=3):
#     u"""
#     Une paire de spline cubiques paramétriques qui interpole ou ajuste le polygone points
#     sx(t), sy(t) sont deux splines.
#     Voir la doc Scipy.interpolate...
#     - si methode dans {'x cubic', 'ius','periodic','interp1d',} c'est de l'interpolation
#     - si methode est {'us','univariatespline'} c'est de l'ajustement, le poids est 1 pour tous les points
#     Retourne:
#     --------
#     T = les valeurs du paramètre t, abscisse curviligne NORMALISEE, entre 0 et 1.
#     sx, sy = les deux splines
#     """
#     if methode in ('periodic', 'p cubic', ) :
#         if all(points[0] == points[-1]) : pass
#         else : #On rajoute le premier point a la fin
#             points = np.vstack((points, points[0]))
# #             points.resize((1+len(points),2)) #moins cher que vstack mais marche pas !!
# #             points[-1] = points[0]
#     if tension == 0.0 :
#         eps = 0.0
#     else :
#         eps = 10.0**(-tension)
#
#     N = len(points)
#     T = absCurv(points)
#     if len(points)<2 : return T, None, None
#     T /= T[-1]
#     X = points[:,0]
#     Y = points[:,1]
#     try : methode = methode.lower()
#     except AttributeError : pass
# #     debug(methode=methode, tension=tension, degre=degre)
#     if methode in ('ius','interpolatedunivariatespline') :
#         try :
#             sx = InterpolatedUnivariateSpline(T, X, k=degre)#s=la précision de l'ajustement s=0 <=> interpolation
#             sy = InterpolatedUnivariateSpline(T, Y, k=degre)
#         except Exception as msg:
#             debug()
#             print unicode (msg)
# #             debug(u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
#             sx = sy = None
#     elif methode in ('us','univariatespline') :
#         try :
#             weights = np.ones(N)
#             W = 1000.0
#             # en supposant que tous les termes erreur di^2=wi*(xi-f(ti))^2 sont egaux
#             # le choix de s suivant implique
#             # abs(xi-f(ti))<eps et
#             # abs(x1-f(t1))<eps/(N*W) et abs(xN-f(tN))<eps/(N*W)
# #             eps = 10.0**(-tension)
#             weights[0] = weights[-1] = W
#             weights /= np.sum(weights)
#             s = eps/(N*W)
#             sx = UnivariateSpline(T, X, w=weights, k=degre, s=s)#s=la précision de l'ajustement s=0 <=> interpolation
#             sy = UnivariateSpline(T, Y, w=weights, k=degre, s=s)
#         except Exception as msg:
#             debug()
#             print unicode(msg)
# #             debug(u'Impossible de calculer la testspline (pas assez de points ?, degré trop élévé ?)')
#             sx = sy = None
#     elif methode in ('interp1d',) :
#         try :
#             sx = interp1d(T, X, kind=degre)
#             sy = interp1d(T, Y, kind=degre)
#         except ValueError as msg:
#             debug()
#             print unicode(msg)
#             sx = sy = None
# #     elif methode in ('periodic',) :
# #         try :
# #             sx = PeriodicSpline(T, X, k=degre, s=eps)
# #             sy = PeriodicSpline(T, Y, k=degre, s=eps)
# #         except ValueError as msg:
# #             debug()
# #             print unicode(msg)
# #             sx = sy = None
#     elif 'cubic' in methode :#or isinstance(methode, (tuple, list, np.ndarray)):
#         if methode == 'p cubic' : bc_type='periodic'
#         elif methode == 'c cubic' : bc_type='clamped'
#         elif methode == 'n cubic' : bc_type='natural'
#         else : bc_type = 'not-a-knot'
#
#         try :
# #             debug(T)
#             sx = CubicSpline(T, X, bc_type=bc_type)
#             sy = CubicSpline(T, Y, bc_type=bc_type)
#         except ValueError as msg:
#             print unicode(msg)
#             sx = sy = None
#
#     elif isinstance(methode, (tuple, list, np.ndarray)):
#         bc_type = methode
#         try :
# #             debug(T)
#             sx = CubicSpline(T, X, bc_type=bc_type)
#             sy = CubicSpline(T, Y, bc_type=bc_type)
#         except ValueError as msg:
#             print unicode(msg)
#             sx = sy = None
#     return T, sx, sy
#
# # def rotate(points,alfa,centre):
# #     u'''
# #     rotation (2d)
# #     de centre 'centre' et d'angle 'alfa' en radians
# #     Retourne une COPIE de 'points' rotationnée(!).
# #     'points' est supposé stocké par ligne (shape=(n,2)), chaque point est de shape (1,2),
# #     il les faudrait en colonne (shape=(2,1)) pour faire le produit matriciel.
# #     Donc on transpose tout et on ecrit Xi' = C' + (Xi'-C')*A' au lieu de
# #     Xi = C + A*(Xi-C), pour i= 0, 1,...
# #     '''
# #     Ct=np.asarray(centre).reshape((1,2))
# #     cosa,sina=np.cos(alfa),np.sin(alfa)
# #     At=np.matrix([[cosa,-sina],[sina,cosa]]).transpose()
# #     Xt=points-Ct
# #     Xt=Xt*At+Ct
# # #     debug(Xt.shape)
# #     return np.asarray(Xt)

def symetriser(points, sym):
    """
    Symetrie in situ de 'points' par rapport à un des plan Oxy, Oyz ou Oxz
    Parametres:
    ----------
        - points = np.ndarray de shape (m,n,...,3) ou (m,n,...,2)
        - sym = un chaine de caracteres pour indiquer l'orientation.
            sym contient au plus 3 caracteres parmi {'x','y','z'}
            si sym contient 'x' les points sont transformés [x,y,z] -> [-x,y,z] (points est symétrisé par rapport au plan Oyz)
            si sym contient 'y' les points sont transformés [x,y,z] -> [x,-y,z] (points est symétrisé par rapport au plan Oxz)
            si sym contient 'z' les points sont transformés [x,y,z] -> [x,y,-z] (points est symétrisé par rapport au plan Oxy)
    """

    shape = points.shape
    dim = shape[-1]#2d ou 3d
    if not dim in (2,3) : return points
    vpoints = points.view()
    vpoints.shape = -1,dim # on remet a plat sauf sur la dimension. vpoints est une liste (simple) de points 2D ou 3D
    sym = sym.lower()
    if 'x' in sym : vpoints[:,0] *= -1
    if 'y' in sym : vpoints[:,1] *= -1
    if dim==3 and 'z' in sym : vpoints[:,2] *= -1
    return points

def symetrieAxe(points,axe,centre=0.0):
    '''
    symetrie
    - d'axe 'axe' = 0 (vertical) ou 1 (horizontal)
    - de centre 'centre', réel
    X = 2*centre - X pour axe vertical
    Y = 2*centre - Y pour axe horizontal
    points est modifié in situ.
    '''
#     debug(points.shape, axe, centre)
    points.shape=(-1,2)
    points[:,axe]*=-1
    if centre !=0.0 : points[:,axe]+=2*centre
    return points

def normale(a, b, c, d=None):
    """DIMENSION 3 : calcule la normale à un panel, un plan, un triangle, sens direct"""
    if d is None:#triangle, pas de pb
        v = np.cross(b-a, c-a)

    #les milieux des 4 segments sont coplanaires, on prend la normale à ce plan
    #c'est egalement le produit vectoriel (c-a)vect(b-d)
#     m1, m2, m3 = 0.5*(a+b),0.5*(b+c),0.5*(c+d)#,0.5*(d+a)
    else :
        v = np.cross(a-c, b-d)
    return v/np.linalg.norm(v)

def panels(shape, closed=False):
    """
    un generateur qui parcourt les panels d'un tableau de points de shape 'shape' dans le sens direct
    un generateur renvoit les éléments les uns après les autres SANS créer explicitement la liste des éléments
    => economie memoire
    'shape' represente une famille de nervures consécutives.
    Tous les tableaux 'points((nbn, nbpn,3))' de Axile ont cette shape.
    Le tableau peut representer des nervures fermées (closed=True) en BF (points[k,0] == points[k,-1])
    ou bien ouvertes (closed=False)  (points[k,0] != points[k,-1])
    :param shape : tuple = (nbn, nbpn,3)
    :param closed : bool True si le BF est en double, False sinon
    :return : les panels de 'points' sous la forme de 4 numeros de points A,B,C,D
    --------
    Utilisation :
    -----------
    >>> shape = (10,20,3)
    >>> for panel in panels(nervures, False) :
    >>>     print panel
    """
    nbn = shape[0]#nb nervures
    nbpn = shape[1]#nb points par nervure
    kcais = 0
    while kcais < nbn-1:
        N1 = list(range(kcais*nbpn, (kcais+1)*nbpn))#Des numeros de points
        N2 = list(range((kcais+1)*nbpn, (kcais+2)*nbpn))
        kpan = 0
#         debug(caisson=kcais, kpan=kpan)
        while kpan < nbpn-1 :
            yield (N1[kpan], N2[kpan], N2[kpan+1], N1[kpan+1])
            kpan += 1

        if not closed :#Cas particulier : le dernier panel du caisson kcais
            yield (N1[-1], N2[-1], N2[0], N1[0])
        kcais += 1

def airePanel(a,b,c,d, precis=True):
    """
    Calcul de la surface du panel en DIMENSION 3 (quadrilatere abcd) comme la somme des aires de deux triangles.
        Il y a deux manières de decouper le panel en deux triangles : (BAD et BCD) ou bien (ABC et ACD)
        Les surfaces correspondantes étant S1 et S2, et en supposant que les points sont
        A:[0,0,0]$B:[1,0,0]$C:[0,1,0]$D:[1,1,h],
        Si h<<1 (quadrilatere presque plat), on a
            S1 ~ 1+0.75 h^2, et
            S2 ~ 1+0.25 h^2 et
            S1-S2 ~ 0.5 h^2 et
            (S1-S2)/S2 ~ 0.5 h^2.
    j'utilise S2
    Si precis=True, je calcule l'aire du quadrilatere plan constitué par
    les 4 milieux (ils sont coplanaires) et j'y rajoute les 4 petits triangles.
    """
    if precis :
        m1 = 0.5*(a + b)
        m2 = 0.5*(b + c)
        m3 = 0.5*(c + d)
        m4 = 0.5*(d + a)
        return airePanel(m1, m2, m3, m4, precis=False)\
            + aireTriangle(a, m1, m4) + aireTriangle(b, m2, m1)\
            + aireTriangle(c, m3, m2) + aireTriangle(d, m4, m3)
    else :
        return aireTriangle(a, b, c) + aireTriangle(a, c, d)

def aireTriangle(a,b,c):
    """
    Aire du triangle abc dans l'espace.
    C'est la moitié de la norme du produit vectoriel ab vect ac
    beaucoup plus rapide que son equivalent numpy :
    >>> return norm(cross(b-a, c-a))*0.5
    (sauf vectorisation, bien sûr)
    """
#     print(a.shape)
    u, v = b-a, c-a
    r = u[2]*v[0]-u[0]*v[2]
    s = u[0]*v[1]-u[1]*v[0]
    t = u[1]*v[2]-u[2]*v[1]
    return 0.5*np.sqrt(r*r + s*s + t*t)

def XYZ(points):
    """Pour ne pas avoir à écrire sans arrêt
    X,Y,Z=points[:,0],points[:,1],points[:,2]
    on exrit X,Y,Z=XYZ(points)
    """
    return points[:,0],points[:,1],points[:,2]

def structuredGrid(nx,ny,nz,delta=-1):
    """
    retourne un tableau np.ndarray de shape(nx,ny,nz,3)

    produit tensoriel de X=(x_1,x_2,...x_nx), Y=(y_1,y_2,...y_ny)
    et Z=(z_1,z_2,...z_nz), aléatoirement espacés, mais triés.
    ou bien régulièrement espacés
    Contrainte : nx>1, ny>1, nz>1
    """
#    nx,ny,nz = shape
    if delta==-1 :
        X=np.random.random(nx)
        Y=np.random.random(ny)
        Z=np.random.random(nz)
        X.sort()
        Y.sort()
        Z.sort()
    else :#TODO: delta=dx,dy,dz
        dx,dy=1.0/(nx-1),1.0/(ny-1)
        X=np.arange(0.0,1.0+dx,dx)
        Y=np.arange(0.0,1.0+dy,dy)
    if nz==1 :
        Z=np.asarray([0.0])
    else :
        dz=1.0/(nz-1)
        Z=np.arange(0.0,1.0+dz,dz)
    G=np.ndarray((nx,ny,nz,3),dtype=float)
    for i,x in enumerate(X) :
        for j,y in enumerate(Y) :
            for k,z in enumerate(Z) :
                G[i,j,k,:]=[x,y,z]
    return G

class Rotation(object) :
    """"
    Une classe qui modelise les rotations dans l'espace.
    """
    def __init__(self, alfa, centre, axe) :
        """
        - alfa : float, angle de la rotation (en radians) ou bien (float, float) = (cos(a), sin(a))
        - centre : np.ndarray de shape (3,), le centre de la rotation
        - axe : np.ndarray de shape (3,), l'axe de la rotation (peut être de norme != 1)
        """
        if isinstance(alfa, tuple) :
            self.cosa, self.sina = alfa
            self.alfa = None
        elif isinstance(alfa, (float, int)) :
            self.cosa = cos(alfa)
            self.sina = sin(alfa)
            self.alfa = alfa
#         if centre is None  : rstack()
#         debug(alfa=className(alfa),centre=className(centre),axe=className(axe))
        self._centre = asarray(centre).reshape((1,3))
        self._axe = asarray(axe).reshape((1,3))
        self.computeMatrice()

    def __str__(self) :
        infos = []
        infos.append('cos(a)=%.2g, sin(a) = %.2g'%(self.cosa, self.sina))
        if self.alfa is not None :
            infos.append('alfa (deg) = %.2f'%(self.alfa*180/math.pi))
        infos.append('centre     = %s'%self.centre)
        infos.append('axe        = %s'%self.axe)
        infos.append('Matrice    = \n%s'%self.matrice)
        return '\n'.join(infos)

    def computeMatrice(self) :
        axe = self.axe
#         centre = self.centre
#         debug(axe=axe)
        u = asarray(axe).reshape((3,))
        u = u/LA.norm(u)
        c, s = self.cosa, self.sina
        c1 = 1-c
        c1u0 = c1*u[0]
        c1u1 = c1*u[1]
        c1u01 = c1u0*u[1]
        c1u02 = c1u0*u[2]
        c1u12 = c1u1*u[2]
        su0 = s*u[0]
        su1 = s*u[1]
        su2 = s*u[2]
        R = np.matrix([[c+u[0]*c1u0, c1u01-su2    , c1u02+su1],
                       [c1u01+su2  , c+u[1]*c1u1  , c1u12-su0],
                       [c1u02-su1  , c1u12+su0    , c+u[2]**2*c1]]).transpose()
        self.matrice = R
#         print 'rotation : alfa, cos, sin, axe, centre\n',self.alfa, c, s, axe, centre
#         print 'matrice rotation : =\n',R
#         print 'matrice rotation : inverse-transpose=\n',R.I-R.T
#         print 'det R=',LA.det(R)
        return R

    def apply(self, X):
        """
        Applique la rotation à X. Modifie X lui-même.
        Si R est une instance de Rotation, et X un vecteur (1,3)
        R.apply(X) <==> R(X)
        """
#         debug(centre_X=(self.centre,X))
        X -= self.centre
#         debug(centre_X=(self.centre,X))
        X[:,:] = X*self.matrice
#         debug(centre_X=(self.centre,X))
        X += self.centre
        return X
    __call__=apply
    @property
    def centre(self):
        return self._centre
    @centre.setter
    def centre(self, c):
        raise AttributeError("On ne change pas une rotation, on la reconstruit")
    @property
    def axe(self):
        return self._axe
    @axe.setter
    def axe(self, c):
        raise AttributeError("On ne change pas une rotation, on la reconstruit")

class Polyline3D(object):
    '''Polyline3D. Un spline paramétrique 3d sommaire.
    a vocation à être initialisée, mais  modifiée avec précaution car
    Le calcul des centre, encombrement, périmetre, spline etc.
    N'EST PAS MIS A JOUR
    marker et mode sont des préférences pour le plot()
    parametres autorisés dans __init__() et load()
    '''
#     ALLOWED_KEYS = ['points', 'name','role','parent']
    DEFAULT = dict(points=ndarray((0,3), dtype=float),
                   name='', role='none', parent=None)
    prefs = GraphicsPrefs
    def __init__(self, *args, **kargs):
#     def __init__(self, points=None, role=None, parent=None, name=''):
        self.spline = None
        dump = self.DEFAULT.copy()
        """Constructeur vide Polyline3D()"""
        """Constructeur recopie Polyline3D(une_instance_de_Polyline3D, name=..)"""
        if len(args)>1:
            msg = """
    Un seul argument non nommé autorisé, doit être de type Polyline3D.
        pour constructeur recopie.
        Ici, on a %d arguments non nommés de types : %s"""\
            %(len(args),[className(a) for a in args])
            raise TypeError(msg)
        elif len(args)==1 :
            #Un seul args, de type Polyline3D => constructeur recopie
            if isinstance(args[0], Polyline3D):
                dump.update(args[0].toDump())
            else :
                msg = """
    Un seul argument non nomme autorisé, doit etre de type Polyline3D (pour constructeur recopie.)
        Ici, 1 argument non nommé de type %s"""%(className(args[0]))
                debug(msg)
                raise TypeError
#         else :#len(args)==0
        dump.update(kargs)
        self.load(dump)
#         if not kargs and not args:

#         else :
#         self.load(dump)
#         self.name = name
#         self.parent = parent
        if self.points.shape[0]>2 : self.update()

    def __call__(self, t, k=0):
#         debug(t=t)
#         debug('dx%d(%s) = %s'%(k,t,self.x(t,k)))
        if isinstance(t, (list, tuple, ndarray, array))  :
            return asarray(list(zip(self.x(t,k), self.y(t,k), self.z(t,k))))
        elif isinstance(t,(int, float)) :
            return asarray((self.x(t,k), self.y(t,k), self.z(t,k)))

    @property
    def dpoints(self):
        return self(linspace(0,1,500))

    def toDump(self):
        return dict(points=self.points, role=self.role, name=self.name)

    def load(self, dump):
        """Normalement self.load(self.toDump()) ne modifie pas self"""
        for key, value in list(dump.items()):
            setattr(self, key, value)

    def update(self):
        self._centre = self.isoBarycentre.reshape(1,3)
        self.computeAbsCurv()
        self.computeSpline()

    @property
    def isoBarycentre(self):
        """
        Calcule et retourne l'isobarycentre des masses ponctuelles égales, placées aux points de 'X'
        X est un np.ndarray de shape n,3
        """
        return average(self.points,0)#.reshape(1,3)
    centre=isoBarycentre

    @property
    def perimetre(self):
        """Longueur du contour"""
        return self.abscurv[-1]

    @property
    def x(self):
        """spline x(t)"""
        return self.spline[0]

    @property
    def y(self):
        """spline y(t)"""
        return self.spline[1]

    @property
    def z(self):
        """spline z(t)"""
        return self.spline[2]

    def info(self):
        infos = ['<'+self.__class__.__name__+'>']
        infos.append('    shape = '+str(self.points.shape))
        if self.points.shape[0] <=3 :
            infos.append('    Pas de spline')
            return infos
        (mx,my,mz), (Mx, My, Mz) = self.encombrement()#(mx,my,mz), (Mx, My, Mz)
        infos.append('points:')
        infos.append('    x dans [%.3g, %.3g]'%(mx,Mx))
        infos.append('    y dans [%.3g, %.3g]'%(my,My))
        infos.append('    z dans [%.3g, %.3g]'%(mz,Mz))
        ((mx,my,mz), (Mx, My, Mz)), ((tmx,tmy,tmz), (tMx, tMy, tMz)) = self.encombrement('spline')
        infos.append('splines:')
        infos.append('    x dans [%.3g, %.3g], t = [%.3g, %.3g]'%(mx,Mx, tmx, tMx))
        infos.append('    y dans [%.3g, %.3g], t = [%.3g, %.3g]'%(my,My, tmy, tMy))
        infos.append('    z dans [%.3g, %.3g], t = [%.3g, %.3g]'%(mz,Mz, tmz, tMz))
        infos.append('dx, dy, dz=%s'%(str(self.encombrement1())))

#         infos.append('platitude = %g'%(self.platitude()))#très cher !
        return infos

    def __str__(self):
        return '\n'.join(self.info())

    def encombrement(self, mode='points'):
        if mode=='points' : return encombrement(self.points)
        elif mode=='spline' : return self._encombrementSpline()

    def _encombrementSpline(self):
        """Encombrement sous forme V,T
        V = [(min(x),min(y),min(z)), (max(x),max(y),max(z))]
        T = [(tmx, tmy, tmz), (tMx, tMy, tMz)]
        où
            - tmx est la valeur du parametre t pour lequel x est minimale, tmx dans [0,1],
            - tMx est la valeur du parametre t pour lequel x est maximale, tMx dans [0,1],
        etc...
        """
        VE, TE = np.ndarray((2,3)), np.ndarray((2,3))
        for i in (0,1,2):#,1,2) :
            s = self.spline[i]
            result = minimize_scalar(s,bounds=(0.0,1.0),
                                     method='Bounded',
                                     options={'disp':False, 'maxiter':200})
            VE[0,i], TE[0,i] = result['fun'], result['x']
            result = minimize_scalar(lambda t:-s(t),bounds=(0.0,1.0),
                                     method='Bounded',
                                     options={'disp':False, 'maxiter':200})
            VE[1,i], TE[1,i] = -result['fun'], result['x']
#         self.TE = TE
        return VE, TE

    def encombrement1(self):
        """encombrement sous forme dx, dy, dz"""
        return encombrement1(self.points)

    def translate(self, u):
        """
        Effectue, IN-SITU, une translation de vecteur u et met à jour les
        differents éléments
        :param u': list(float,float,float) [x,y,z]
            ou bien ndarray de shape (1,3) ou (3,1),
            est le vecteur de la translation
        :return u.
        """
        u1 = asarray(u).reshape((3,))
        self.points += u1
        self.update()
        return u1

    def rotate(self, alfa=None, centre=None, axe=None, rot=None):
        """
        Effectue, IN-SITU, une rotation (de centre C='centre', d'axe u='axe', d'angle a='alfa'),
            pour tous les points de X='points'.
            ou bien la rotation rot
        :param alfa: float, l'angle en radians
        :param centre : list de 3 reels [x,y,z] ou ndarray((1,3)) ou ndarray((3,1))
            c'est le centre de la rotation
        :param axe: list de 3 reels [x,y,z] ou ndarray de shape (1,3) ou (3,1),
            est l'axe de la rotation, non normalisé.
        Retourne la Rotation.
        --------
        """
        if rot is None :
            rot = Rotation(alfa, centre, axe)
        elif not isinstance(rot, Rotation) :
            msg = self.__doc__+"""
            le parametre 'rot' est de type %s, il doit etre de type Rotation.
            sinon, precisez les trois parametres 'alfa' , 'centre' et 'axe'
            de la rotation"""%(className(rot))
            raise TypeError(msg)
        self.points = rot(self.points)
        self.update()
        return rot
    nbappelsplatitude=0

    def platitude(self):
        '''
        ATTENTION, coute cher !!
        Détermine un coefficient de platitude du polyligne.
        (platitude=0.0 si les points sont coplanaires)
        On parcours tous les tetrahedres T(i,j,k)=A,Xi,Xj,Xk. où A=X[0]
        On calcule le volume Vt(i,j,k) de T(i,j,k).
        On calcule la sphere circonscrite S(i,j,k) et son volume Vs(i,j,k)

        Si T est grand (4 aretes de meme longueur) alors Vs=Vt.3.pi.sqrt(3).
        T est 'plat' si 3.pi.sqrt(3).Vt/Vs est petit devant 1.
        On fait la moyenne des coefficients a(i,j,k) = 3.pi.sqrt(3).Vt(i,j,k)/Vs(i,j,k)
        '''
        self.nbappelsplatitude += 1
        if self.nbappelsplatitude <10 :
            stack('nbappelsplatitude=%d'%self.nbappelsplatitude)
        n = len(self.points)
#         print n
        X = self.points
        if n<=3 : return 0.0
        A = X[0]
        moy = 0.0
        nt = 0
        cs, cv = 4.0*math.pi/3.0, 3*math.pi*math.sqrt(3)
        for i in range(1, n) :
            B = X[i]
            AB = B-A
            for j in range(i+1, n):
                C = X[j]
                AC = C-A
                for k in range(j+1, n) :
                    D = X[k]
                    AD = D-A
                    nt += 1
                    D = X[k]
                    [I, r] = sphereCirconscriteA4Points(A, B, C, D)
                    if I is not None :
                        Vs = cs*r**3
                        Vt = abs(det(AB, AC, AD))/6.0
                        moy += cv*Vt/Vs
        self._platitude = moy/nt
        return moy/nt
#scipy.interpolate.SmoothBivariateSpline(x, y, z, w=None, bbox=[None, None, None, None], kx=3, ky=3, s=None, eps=None)

    def compute2DSpline(self, interpolated=[0,-1]) :
        "TODO : Calcul spline 2d ???Je ne sais pas à quoi ca sert, ni ce que j'ai voulu essayer..."
        X = self.points[:,0]
        Y = self.points[:,1]
        Z = self.points[:,2]
        maxw = 100.0
        weight = ones(len(self.points))
        for index in interpolated :
            weight[index] = maxw
        spl2d = SmoothBivariateSpline(X, Z, Y, w=weight, kx=3, ky=3, s=len(X)*1.0e-8)
        self.spline2d = spl2d

    def computeSpline(self) :
        """
        Calcule d'une spline cubique X(t),Y(t),Z(t) avec t dans [0,1].
        qui interpolte self.points
        """
        t = self.abscurv/self.abscurv[-1]
        X = self.points[:,0]
        Y = self.points[:,1]
        Z = self.points[:,2]
        x = InterpolatedUnivariateSpline(t, X, k=3)
        y = InterpolatedUnivariateSpline(t, Y, k=3)
        z = InterpolatedUnivariateSpline(t, Z, k=3)
        self.spline = [x, y, z]

    def computeAbsCurv(self ):
        self.abscurv = absCurv(self.points)
#         debug(len_abscurv=len(self.abscurv), len_points=len(self.points))

    def isFlat(self, eps=1.0e-8):
        return self.platitude()<eps

    def writeVTK(self, filename=sys.stdout):
#         return
        if isinstance(filename, str) :
            filename = Path(filename).abspath()
#             debug('sauvegarde %s'%filename)
        nbp = len(self.points)
        cells = list(zip(list(range(nbp-1)), list(range(1,nbp))))
        vtkpd = toVtkPolyData(self.points, cells, type_='lines')
        writer = vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()
        writer.SetFileName(filename)
        writer.SetInputData(vtkpd)
        writer.Update()
        writer.Write()

    def plot(self, show=True, filename=None):
        data = []
        points = self.points.view()
        points.shape=(-1,3)
        X, Y, Z = XYZ(points)
        mode = self.prefs.mode
        marker = self.prefs.marker
        line = self.prefs.line

        goP = go.Scatter3d(dict(x=X, y=Y, z=Z, text=list(range(len(X)))),
                           mode   = mode[self.role],
                           marker = marker[self.role],
                           line   = line[self.role],
                           name   = self.role,
                           hoverinfo='name',
                           showlegend=True
                           )
        data.append(goP)
        if show :
            layout = go.Layout(width = 1800, height = 1200,
                               title = self.name,
                               legend=dict(x=0.2,y=0.2),
                               hovermode='closest')
            if filename is None :
                try : filename = self.pjdir/self.session+className(self)+'.html'
                except : filename = RUNS_DIR/'trash'/className(self)+'.html'
            figure = go.Figure(data=data, layout=layout)
            plot(figure, show_link=True,
                 link_text='Exporter vers plot.ly (pour édition)', filename=filename)
            try : self.files.add(filename)
            except : pass
            return filename
        else :
            return data

def plotConnexions3D(A, C, session='default', show=True, data=[],
                     key='connexions', filename=None, withnodes=True):
    noeuds = A
    connexions = C
    mode = GraphicsPrefs.mode
    marker = GraphicsPrefs.marker
    data = data
    if withnodes :
        data.extend(plotNodes3D(noeuds, data=data, show=False))
    #construction des cellules vraies (avec des points3D)
    _, nbp = connexions.shape#nb de cellules, nb de points par cellule
    if nbp>2 :
        #si nbp>2  (e.g. tri ou quad), un point suppl. pour refermer le polyline
        nbp += 1
    for k, I in enumerate(connexions) :
        cell = ndarray((nbp, 3), dtype=float)
        for kp, i in enumerate(I) : #kp-ieme point de cell
            cell[kp] = noeuds[i]
        #fermeture cell si nbp>2
        if nbp>2 : cell[1+kp] = noeuds[I[0]]

        goc = go.Scatter3d(x=cell[:,0], y=cell[:,1], z=cell[:,2],
                           text        = 'Cell %d => '%k+repr(I.tolist()),
                           mode        = mode[key],
                           marker      = marker[key],#{'size':3, 'color':'rgb(255,100,100)'},
                           name        = 'connexions' if k==0 else str(k),
                           hoverinfo   = 'text',
                           showlegend  = k==0,)
#                            legendgroup = 'connexions')
        data.append(goc)
    if show :
        if filename is None :
            filename = RUNS_DIR/'trash'/'trash.html'
        #Pour aspect ratio
        layout = go.Layout(width = 1800, height = 1200,
                           title = 'connexions ' + session,
                           scene = dict(xaxis = dict(title = 'x'),
                                        yaxis = dict(title = 'y'),
                                        zaxis = dict(title = 'z'),
    #                        legend = dict(x=0.9,y=0.9),
                                        hovermode = 'closest'))

        figure = go.Figure(data=data, layout=layout)
        plot(figure, show_link=True,
             link_text='Exporter vers plot.ly (pour édition)', filename=filename)
        return filename
    else :
        return data

def plotNodes3D(nodes, session='default', show=True, data=[],
                key='nodes', filename=None):
    mode = GraphicsPrefs.mode
    marker = GraphicsPrefs.marker
    data = data
    X,Y,Z = XYZ(nodes)
#     text = ["noeud %d"%i for i in range(len(nodes))]
    goc = go.Scatter3d(x=X, y=Y, z=Z,
                       mode        = mode[key],
                       marker      = marker[key],#{'size':3, 'color':'rgb(255,100,100)'},
#                        text = text,
#                        hoverinfo   = 'text',
                       showlegend  = True,
                       legendgroup = 'noeuds')
    data.append(goc)
    if show :
        if filename is None :
            filename = RUNS_DIR/'trash'/'trash.html'
        #Pour aspect ratio
        layout = go.Layout(width = 1800, height = 1200,
                           title = 'connexions ' + session,
                           scene = dict(xaxis = dict(title = 'x'),
                                        yaxis = dict(title = 'y'),
                                        zaxis = dict(title = 'z'),
    #                        legend = dict(x=0.9,y=0.9),
                                        hovermode = 'closest'))

        figure = go.Figure(data=data, layout=layout)
        plot(figure, show_link=True,
             link_text='Exporter vers plot.ly (pour édition)', filename=filename)
        return filename
    else :
        return data

def det(X,Y,Z):
    """Determinant dans R^3"""
    return X[0]*Y[1]*Z[2] + X[2]*Y[0]*Z[1] + X[1]*Y[2]*Z[0]\
         - X[2]*Y[1]*Z[0] - X[1]*Y[0]*Z[2] - X[0]*Y[2]*Z[1]

###########################################################################################
###########################################################################################
###### Fonctions
###########################################################################################
###########################################################################################

# def locate(t, T, eps = 1.0e-5):
#     """
#     Localise t dans T=ndarray(n).
#     On doit avoir T[0]<=t<=T[-1] et T croissant, len(T)>=2
#     retourne (True_ou_False, index),
#     - (True, k) si t=T[k] à eps près (t est une des valeurs T[k])
#     - (False, k) si T[k] < t < T[k+1] (t n'est pas un T[k])
#     """
#     g = where(T<=t+eps)[0]
#     d = where(T>=t-eps)[0]
#     return (True, g[-1]) if g[-1]==d[0] else (False, g[-1])

def prodVect(u,v):
    ''' Produit vectoriel des vecteurs u et v '''
    u1,u2,u3 = u[0],u[1],u[2]
    v1,v2,v3 = v[0],v[1],v[2]
    w1 = u2*v3-u3*v2
    w2 = u3*v1-u1*v3
    w3 = u1*v2-u2*v1
    return asarray([w1,w2,w3]).reshape((1,3))

def ProdVect3d(u, v):
    ''' Produit vectoriel de deux vecteurs 3D '''
    u1,u2,u3 = u[0],u[1],u[2]
    v1,v2,v3 = v[0],v[1],v[2]
    return asarray([u2*v3-u3*v2, u3*v1-u1*v3, u1*v2-u2*v1])

def encombrement(piece, dim=3):
    """
    piece doit être un np.ndarray de shape (N, dim) quelconque,
    retourne le parallelepipede d'encombrement du nuage de points
    """
    if isinstance(piece, (list, tuple )) :
        return encombrement(np.asarray(piece),dim)#recursif
    elif isinstance(piece,np.ndarray):
        Max, Min = np.max, np.min #ca reste local
        if dim == 1 :
            return Min(piece), Max(piece)
        points = piece.view()
        points.shape = -1, dim
        if dim == 3 :
            M = [Max(points[:,0]), Max(points[:,1]), Max(points[:,2])]
            m = [Min(points[:,0]), Min(points[:,1]), Min(points[:,2])]
#             M = np.asarray([Max(points[:,0]),Max(points[:,1]),Max(points[:,2])])
#             m = np.asarray([Min(points[:,0]),Min(points[:,1]),Min(points[:,2])])
            return(m,M)
        elif dim == 2 :
            M = [Max(points[:,0]), Max(points[:,1])]
            m = [Min(points[:,0]), Min(points[:,1])]
            return(m,M)
    else :
        raise NotImplementedError

def encombrement1(piece, dim=3):
    """retourne l'encombrement de la piece sous forme dx[,dy[,dz]]"""
    m, M = encombrement(piece, dim)
    return np.asarray(M) - np.asarray(m)

def prodScal(u,v):
    ''' Produit scalaire des vecteurs u et v '''
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def rotated(points, alfa, centre, axe):
    '''
    Retourne une COPIE de points rotationnée(!).
    Arguments
    ---------
    :param alfa : reel, l'angle en radians
    :param centre : liste de 3 reels [x,y,z] ou ndarray de shape (1,3) ou (3,1),
        c'est le centre de la rotation
    :param axe : liste de 3 reels [x,y,z] ou ndarray de shape (1,3) ou (3,1),
        est l'axe de la rotation, non normalisé.
    Fonctionnement:
    --------------
    X est une copie de points
    La rotation est définie, pour M=X[i], i=1,2,... par
        M <- C + R(u,a).(M-C), ou R(u,a) est la matrice de rotation
    X est supposé stocké par ligne (shape=(n,3)), chaque point est de shape (1,3),
    il les faudrait en colonne (shape=(2,1)) pour faire le produit matriciel.
    Donc on transpose tout et au lieu de M <- C + R(u,a)*(M-C), on ecrit (' designe la transposition)
        M' <- C' + (M'-C')*R(u,a)' , pour M=X[i], i= 0, 1,...
    La matrice de la rotation est R(u,a) =
    \begin{array}{ccc}
    \cos\theta+u_{x}^{2}\left(1-\cos\theta\right) & u_{x}u_{y}\left(1-\cos\theta\right)-u_{z}\sin\theta & u_{x}u_{z}\left(1-\cos\theta\right)+u_{y}\sin\theta\\
    u_{y}u_{x}\left(1-\cos\theta\right)+u_{z}\sin\theta & \cos\theta+u_{y}^{2}\left(1-\cos\theta\right) & u_{y}u_{z}\left(1-\cos\theta\right)-u_{x}\sin\theta\\
    u_{z}u_{x}\left(1-\cos\theta\right)-u_{y}\sin\theta & u_{z}u_{y}\left(1-\cos\theta\right)+u_{x}\sin\theta & \cos\theta+u_{z}^{2}\left(1-\cos\theta\right)
    \end{array}
    '''
    X = points.copy()
    u = asarray(axe).reshape((3,))
    u = u/LA.norm(u)
    Ct = asarray(centre).reshape((1,3))
    c, s = cos(alfa), sin(alfa)
    c1 = 1-c
    c1u0 = c1*u[0]
    c1u1 = c1*u[1]
    c1u01 = c1u0*u[1]
    c1u02 = c1u0*u[2]
    c1u12 = c1u1*u[2]
    su0 = s*u[0]
    su1 = s*u[1]
    su2 = s*u[2]
    Rt = matrix([[c+u[0]*c1u0, c1u01-su2    , c1u02+su1],
                    [c1u01+su2  , c+u[1]*c1u1  , c1u12-su0],
                    [c1u02-su1  , c1u12+su0    , c+u[2]**2*c1]]).transpose()
    X -= Ct
    X[:,:] = X*Rt
    X += Ct
    return X

# def dist(p, q, n=2):
#     '''retourne la distance de p1 à p2 en norme n=2'''
#     return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)
#
# def dist2(p, q, n=2):
#     '''retourne la distance de p1 à p2 en norme n=2'''
#     return (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2

def sphereCirconscriteA4Points(A,B,C,D):
    '''I = centre, R = rayon'''
    AB, AC, AD = B-A, C-A, D-A
    M = asarray([[AB[0], AB[1], AB[2]],
                  [AC[0], AC[1], AC[2]],
                  [AD[0], AD[1], AD[2]]])
    b = asarray([prodScal(AB, A+B), prodScal(AC, A+C), prodScal(AD, A+D)])/2.0
    try :
        I = LA.solve(M,b)
        return I, LA.norm(I-A)
    except LA.LinAlgError :
        return (None, None)

def ellipse(centre, axes, nbp=20):
    if len(centre)==2 : centre = [centre[0], centre[1], 0.0]
    a,  b= axes
    T = np.linspace(0.0, 2*pi, nbp)
    T[-1] = 0.0
    return asarray(list(zip(a*cos(T), b*sin(T), len(T)*[0.0])))+centre

def cercle(centre, r, nbp=20):
    return ellipse(centre, (r,r), nbp)

def volumeTetraedre(A,B,C,D):
    """C'est det(AB,AC,AD)/6"""
    u, v, w = B-A, C-A, D-A
    return abs(u[0]*v[1]*w[2]+u[2]*v[0]*w[1]+u[1]*v[2]*w[0]\
            - (u[2]*v[1]*w[0]+u[1]*v[0]*w[2]+u[0]*v[2]*w[1]))/6

def volumeTetraedre1(A,B,C,D):
    """C'est det(AB,AC,AD)/6"""
    return det(B-A,C-A,D-A)

def surfaceTriangle(a,b,c):
    raise NotImplementedError('surfaceTriangle() : Utiliser plutôt aireTriangle(), 2 à 10 fois plus rapide')
    return norm(cross(b-a, c-a),axis=1)*0.5

if __name__=="__main__":
    A, B, C, D = random.rand(4, 3)
    print(A)
#     A,B,C = asarray([[[0,0,0]],[[1,0,0]],[[0,1,0]]])
#     print (A.shape)
#     mexit()
    t = TicToc()
#     t.tic()
#     for i in range(1000) : NN = surfaceTriangle(A, B, C)
#     t.toc(restart=True)
#     t.tic()
#     for i in range(1000) : AA = aireTriangle(A, B, C)
#     t.toc()
#     print('NN =',norm(NN-AA))
#     print('AA =',AA)
    t.tic()
    for i in range(1000) : V = volumeTetraedre(A, B, C, D)
    t.toc(restart=True)
    for i in range(1000) : V1 = volumeTetraedre1(A, B, C, D)
    t.toc()
    mexit(V-V1)
    from tests.utils.testsgeometry3d import Tests
    t = Tests()
#     t.run()

#     if 1 : t.testPointsDoubles()
#     if 1 : t.testPointsDoubles1()
#     if 1 : t.testPlacementNervureI(None)
#     if 1 : t.testRotation()
    if 1 : t.testPanels()
#     if 1 : t.testDivers()
#     if 1 : t.testDivers1()
#     if 1 : t.testPolyline3D()
    exit()


    import cProfile, pstats
    filename = OUTPUT/'profile.stats'
    cProfile.run('main()',filename)
    stats=pstats.Stats(filename)
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.sort_stats('time')
    stats.print_stats()


