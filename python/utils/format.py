#!/usr/bin/python
#-*-coding: utf-8 -*-

'''
Axile -- Outil de conception/simulation de parapentes Nervures

    Arrondir une donnée avec n décimales

@author:     Michel Le Berre
@copyright:  2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com
@deffield    updated: 28 Mai 2013
'''
import PyQt5
import numpy as np
from numpy import ndarray, nan, Infinity, infty, array, asarray, NaN, Inf
from numbers import Number
from utils.debog import rdebug, mexit
from PyQt5.QtCore import QPointF
from math import inf
import math

# print eval('PyQt5.QtCore.QString')
# exit()

''' formatage données de m à mm avec 6 décimales '''
def formDatamm(val_in):
    fdec=1000000.0
    val_out=int(val_in*1000.0*fdec+0.5)/fdec
    return val_out

''' formatage données avec ndec décimales '''
def formData(val_in,ndec):
    '''
    Peut s'écrire:
    fmt='%%.%df'%ndec
    float('%.5f'%val_in)
    '''
    if np.isnan(val_in):
        print('Erreur: NaN in formData')
        return val_in
    else:
        if ndec>=1:
            fdec=10.0**ndec
            sig=np.sign(val_in)
            val_out=sig*int(abs(val_in)*fdec+0.5)/fdec#+0.5 pour arrondir à l'entier le plus proche
        else:#ndec==0
            sig=np.sign(val_in)
            val_out=sig*int(abs(val_in)+0.5)
        return val_out

def fmtOld(r, prec=3):
    """(récursif) Formate un tableau ou dict avec prec décimales"""
    if isinstance(r,(ndarray, list, tuple)) :
        return [fmt(x,prec) for x in r]
    elif isinstance(r, dict):
        return dict( (key,fmt(value)) for (key, value) in r.items())
    elif isinstance(r, Number) :
        prec = 10.0**prec
        return round(prec*r)/prec
    else :
        return r
def fmt(r, prec=3):
    """(récursif) Formate un tableau ou dict avec prec décimales"""
    cls = type(r)
    if cls==int :
        return r

    if issubclass(cls, Number) :
        if np.isnan(r) or math.isnan(r) or np.isinf(r) or math.isinf(r) :
            return r
        prec = 10.0**prec
        return round(prec*r)/prec
    elif issubclass(cls, (list, tuple)) :
        return cls([fmt(x,prec) for x in r])
    elif isinstance(r, ndarray) :
        return asarray([fmt(x,prec) for x in r])
    elif isinstance(r, QPointF) :
        return fmt((r.x(),r.y()),prec)
    elif issubclass(cls, dict):
        return dict((key,fmt(value,prec)) for (key, value) in r.items())
    else :
        return r

def toUnicode(r):
    """(récursif) Formate un tableau ou dict en unicode
    pour remplacer les QString('toto') par u'toto' => cf suspentage"""
    cls = type(r)
#     if cls==PyQt5.QtCore.QString :
#         return str(r)
    if cls==PyQt5.QtCore.QPointF :
        return cls([toUnicode(x) for x in r])
    elif issubclass(cls, (list, tuple)) :
        return cls([toUnicode(x) for x in r])
    elif isinstance(r, ndarray) :
        return asarray([toUnicode(x) for x in r])
    elif issubclass(cls, dict):
        return dict((key,toUnicode(value)) for (key, value) in r.items())
    else :
        return r

if __name__=='__main__' :
    X = [[1.268811468596966, 0.3365727279966774], [1.1390328407287598, 0.07332829385995865]]
    X = [[-0.5279636979103088, nan], [0.7137287259101868, -infty]]
    X = [[-0.5279636979103088, nan], [0.7137287259101868, -infty], [1.268811468596966, 0.3365727279966774], [1.1390328407287598, 0.07332829385995865], [1.2571670716747245, 0.2148051133421408], [1.2206038430660453, -0.0507648238639925], [1.5545598268508911, -0.0048885527066886425]]

    print(fmt(X,4))
#     exit()
    print(fmt(tuple(X)))
    print(fmt(asarray(X)))
    D = {'Cx0': 0.02, 'Cxp': 0.01, 'Cxs': 0.016, 'S': 28.201750906217605, 'b': 5.897692607389747, 'c': 0.28999999986302366, 'd': 0.0, 'l': 2.89304324, 'mp': 85.0, 'mv': 3.0, 's': 7.3726443586, 'xGv': 0.9539636527457533, 'zGv': 1.7292584212368416, 'zs': -2.578503039144598}
    print(fmt(D))
#     mexit()
    D = {'nbpes': [7, 11], 'ts': [0.0, 0.08280030143208625, ], 'modes': ['courbure']}
    print(fmt((D,D)))
    print(fmt(dict(dico1=D, dico2=D)))
    methode = ('cubic', ((1, 0, -0.25), (2, 0, 0)))
    print(fmt(methode))
    modech = {'QPointF':PyQt5.QtCore.QPointF(-1.0020506804808074, 0.0025114052142376124),'nbpes': [7, 11, 33, 15], 'ts': [0.0, 0.08280030143208625, 0.18859864810049215, 0.8996632881984424, 1.0], 'modes': ['courbure', 'courbure', 'courbure', 'linear']}
    print(fmt(modech))
    I = 5*[0,1]
    print(fmt(I))
    print(fmt(-9.717034075195107e+15))
#     exit()
#     base_chem = [[[['1'], [PyQt5.QtCore.QString('1.1'), PyQt5.QtCore.QString('1.2'), PyQt5.QtCore.QString('1.3')], [PyQt5.QtCore.QString('1.1.1'), PyQt5.QtCore.QString('1.1.2'), PyQt5.QtCore.QString('1.2.1'), PyQt5.QtCore.QString('1.2.2'), PyQt5.QtCore.QString('1.3.1'), PyQt5.QtCore.QString('1.3.2')], [PyQt5.QtCore.QString('1.1.1.1'), PyQt5.QtCore.QString('1.1.1.2'), PyQt5.QtCore.QString('1.1.1.3'), PyQt5.QtCore.QString('1.1.2.1'), PyQt5.QtCore.QString('1.1.2.2'), PyQt5.QtCore.QString('1.1.2.3'), PyQt5.QtCore.QString('1.2.1.1'), PyQt5.QtCore.QString('1.2.1.2'), PyQt5.QtCore.QString('1.2.2.1'), PyQt5.QtCore.QString('1.2.2.2'), PyQt5.QtCore.QString('1.3.1.1'), PyQt5.QtCore.QString('1.3.1.2'), PyQt5.QtCore.QString('1.3.2.1'), PyQt5.QtCore.QString('1.3.2.2')], [], []], [['1'], [PyQt5.QtCore.QString('1.1'), PyQt5.QtCore.QString('1.2'), PyQt5.QtCore.QString('1.3'), PyQt5.QtCore.QString('1.4')], [PyQt5.QtCore.QString('1.1.1'), PyQt5.QtCore.QString('1.1.2'), PyQt5.QtCore.QString('1.2.1'), PyQt5.QtCore.QString('1.2.2'), PyQt5.QtCore.QString('1.3.1'), PyQt5.QtCore.QString('1.3.2'), PyQt5.QtCore.QString('1.4.1'), PyQt5.QtCore.QString('1.4.2')], [PyQt5.QtCore.QString('1.1.1.1'), PyQt5.QtCore.QString('1.1.1.2'), PyQt5.QtCore.QString('1.1.1.3'), PyQt5.QtCore.QString('1.1.2.1'), PyQt5.QtCore.QString('1.1.2.2'), PyQt5.QtCore.QString('1.1.2.3'), PyQt5.QtCore.QString('1.2.1.1'), PyQt5.QtCore.QString('1.2.1.2'), PyQt5.QtCore.QString('1.2.2.1'), PyQt5.QtCore.QString('1.2.2.2'), PyQt5.QtCore.QString('1.3.1.1'), PyQt5.QtCore.QString('1.3.1.2'), PyQt5.QtCore.QString('1.3.2.1'), PyQt5.QtCore.QString('1.3.2.2')], [], []], [['1'], [PyQt5.QtCore.QString('1.1'), PyQt5.QtCore.QString('1.2'), PyQt5.QtCore.QString('1.3')], [PyQt5.QtCore.QString('1.1.1'), PyQt5.QtCore.QString('1.1.2'), PyQt5.QtCore.QString('1.2.1'), PyQt5.QtCore.QString('1.2.2'), PyQt5.QtCore.QString('1.3.1'), PyQt5.QtCore.QString('1.3.2')], [PyQt5.QtCore.QString('1.1.1.1'), PyQt5.QtCore.QString('1.1.1.2'), PyQt5.QtCore.QString('1.1.1.3'), PyQt5.QtCore.QString('1.1.1.4'), PyQt5.QtCore.QString('1.1.1.5'), PyQt5.QtCore.QString('1.1.2.1'), PyQt5.QtCore.QString('1.1.2.2'), PyQt5.QtCore.QString('1.1.2.3'), PyQt5.QtCore.QString('1.1.2.4'), PyQt5.QtCore.QString('1.1.2.5'), PyQt5.QtCore.QString('1.2.1.1'), PyQt5.QtCore.QString('1.2.1.2'), PyQt5.QtCore.QString('1.2.1.3'), PyQt5.QtCore.QString('1.2.1.4'), PyQt5.QtCore.QString('1.2.2.1'), PyQt5.QtCore.QString('1.2.2.2'), PyQt5.QtCore.QString('1.2.2.3'), PyQt5.QtCore.QString('1.2.2.4'), PyQt5.QtCore.QString('1.3.1.1'), PyQt5.QtCore.QString('1.3.1.2'), PyQt5.QtCore.QString('1.3.2.1'), PyQt5.QtCore.QString('1.3.2.2')], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [['1'], [PyQt5.QtCore.QString('1.1')], [PyQt5.QtCore.QString('1.1.1'), PyQt5.QtCore.QString('1.1.2'), PyQt5.QtCore.QString('1.1.3')], [PyQt5.QtCore.QString('1.1.1.1'), PyQt5.QtCore.QString('1.1.1.2'), PyQt5.QtCore.QString('1.1.1.3'), PyQt5.QtCore.QString('1.1.2.1'), PyQt5.QtCore.QString('1.1.2.2'), PyQt5.QtCore.QString('1.1.3.1'), PyQt5.QtCore.QString('1.1.3.2'), PyQt5.QtCore.QString('1.1.3.3')], [PyQt5.QtCore.QString('1.1.2.1.1'), PyQt5.QtCore.QString('1.1.2.1.2'), PyQt5.QtCore.QString('1.1.2.2.1'), PyQt5.QtCore.QString('1.1.2.2.2')], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []]], [[[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []], [[], [], [], [], [], []]]]
#     print(toUnicode(base_chem))
