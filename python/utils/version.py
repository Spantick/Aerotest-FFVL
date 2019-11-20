#!/usr/local/bin/python2.7
# encoding: utf-8
'''
NervuresDesignTool -- Outil de conception/simulation de parapentes Nervures

Classe xxx
Description :
@module : utils.version
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 16 mars 2013
'''
# import sys, os

import os,sys,platform, subprocess
from pprint import pprint
from aperoconfig import SOURCES_DIR
from path import Path
import subprocess, shlex
from utils.debog import debug
class Versions(object):
    def __init__(self):
#         """conda list -n Axile --export > conda-packages.txt
#             /bin/vikings -input eggs.txt -output "spam spam.txt" -cmd "echo '$MONEY'"""
# #         cmd = 'conda list -n py3 --export > "conda-packages.txt"'
# #         print (shlex.split(cmd))
#
# #         exit()
#         subprocess.run(["ls"],capture_output=True)
#         exit()
#         filename = SOURCES_DIR)/'conda-packages1.txt'
#         condapath =  os.path.expanduser('~'))/'anaconda2/condabin/conda'
#         cmd = '%s list -n py3 --export > %s'%(condapath,filename)
#         print(cmd)
#         args = shlex.split(cmd)
#         print(args)
#         subprocess.run(args)
#         exit()
#         print(cmd)
#         res = os.system(cmd)
#         print(res)
#         exit()
        exes = ('ffmpeg', 'paraview', 'cmarc'   )
        self.dexes = dict([(exe,None) for exe in exes])
        packages = self.packages = ['matplotlib', 'numpy','scipy',
                                    'vtk','dxfwrite','dxfgrabber','shapely',
                                    'cPickle','path','plotly']
        self.lpack = [('Plateform', True, '',platform.platform(),['']),
                      ('Python', True, '',sys.version.replace('\n', ' -- '),[sys.executable])
                      ]
        for p in packages:
            try :
                exec('import '+p)
                value, msg = (True, '')
            except Exception as msg :
                value, msg = (False, str(msg))
            debug(package=p,value=value, msg=str(msg))
            if value :
                if p in ('numpy', 'scipy', 'matplotlib','shapely')  :
                    result = (p, True, '',
                              eval(p+'.__version__'),
                              eval(p+'.__path__'))
                elif p in ('vtk',)  :
                    result = (p, True, '',
                              eval(p+'.vtkVersion.GetVTKVersion()'),
                              eval(p+'.__path__'))
                elif p in ('dxfgrabber'):#,'dxfwrite')  :
                    result = (p, True, '',
                              eval(p+'.version'),
                              eval(p+'.__path__'))
                elif p in ('path','cPickle')  :
                    result = (p, True, '',
                              eval(p+'.__version__'),
                              [str(eval(p+'.__file__'))])
                else :
                    try :
                        v=eval(p+'.__version__')
                    except :
                        try : v=eval(p+'.version')
                        except : v = '???'
                    try :
                        pp = eval(p+'.__path__')
                    except :
                        try : pp = [str(eval(p+'.__file__'))]
                        except : pp = '???'
                    result = (p, True, '',
                              v,
                              pp)

            else :
                result = p, value, msg, None, None
#             print 'result => ', result
            self.lpack.append(result)
#         for pack in self.lpack :
#             print '%+20s = '%pack[0], pack[1:]

#         dversions={}
#         dversions['01Plateforme  '] = sys.platform
#         dversions['02Python      '] = sys.version.replace('\n', ' -- ')#platform.python_version()
#         if pyqt_     :    dversions['03PyQt4 / Qt  '] = ' / '.join([PyQt4.Qt.PYQT_VERSION_STR,PyQt4.QtCore.QT_VERSION_STR])
#         if pyside_   :    dversions['04PySide / Qt '] = ' / '.join([PySide.__version__,PySide.QtCore.__version__])
#         if sys_        :    dversions['13sys         '] = sys.version.replace('\n', ' -- ')
#         if enthought_ : dversions['6.   Enthought   '] = enthought.version
#         self.dversions=dversions
#         dpaths={}
#        dpaths['Python        '] = os.environ['PYTHONPATH']
#         dpaths[                 '01Python        '] = [sys.executable]



    def __str__(self):

        lstring=[16*'#'+'\n####  AXile ####\n'+16*'#']

        lstring.extend([
            '\n#######   sys.path  #######',self.syspaths(),
            '\n#######   versions  #######',self.versions(),
            '\n#######    paths    #######',self.paths(),
            '\n#######   MANQUANT  #######',self.manquants(),
            '\n####### executables #######',self.exes(),])
        lstring.extend([25*'#'])
        print((lstring[7:9]))

        return '\n'.join(lstring)

    def exes(self):
        os.environ['PATH'] = os.environ['PATH']+':/usr/local/bin'
        stringlist = os.environ['PATH'].split(':')
        for key,value in sorted(self.dexes.items()) :
            stringlist.append('%+10s : %-10s'%(key,str(value)))
        return '\n'.join(stringlist)

    def manquants(self):
        stringlist=[]
        for pack, exists, msg, _, _ in self.lpack :
            if exists : continue
            stringlist.append('%+10s : %-10s'%(pack, msg))
        return '\n'.join(stringlist)

    def paths(self):
        stringlist=[]
        for pack, exists,_, _, p in self.lpack[1:] :
            if not exists : continue
            stringlist.append('%+10s : %-10s'%(pack,str(p)))
        return '\n'.join(stringlist)

    def syspaths(self,sort=True):
        paths = ['[%s]'%p for p in sys.path]
        if sort : return '\n'.join(sorted(paths))
        else  : return '\n'.join(paths)

    def versions(self):
        stringlist = []
        for pack, exists, msg, version, _ in self.lpack :
            if not exists : continue
            s = '%+10s : %-10s'%(pack,str(version))
            if msg : s = s + ' (=> Bizarre : %s)'%msg
            stringlist.append(s)
        return '\n'.join(stringlist)

if __name__=='__main__':
    print((Versions()))
