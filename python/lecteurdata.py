#!/usr/local/bin/python2.7
# encoding: utf-8
'''
NervuresDesignTool -- Outil de conception/simulation de parapentes Nervures

Classe : LecteurData
Description : lecteur générique pour fichier paramètres.
@module : inout.lecteurdata
@author:      puiseux
@copyright:   2013 Nervures. All rights reserved.
@license:    LGPL
@contact:    be@nervures.com, pierre@puiseux.name
@deffield    creation: 18 janv. 2013
'''
import sys
from inout.lecteur import Lecteur
from utils.debog import whoami
# from path import Path
from aperoconfig import DATA_DIR
from numpy import asarray

def readPyData(filename):
    """
    Lit, évalue (avec eval) une EXPRESSION contenu dans le fichier filename.
    ATTENTION une INSTRUCTION (p.ex. x=3) n'est pas une EXPRESSION
    :returns: le résultat de l'expression lue et évaluée

    :example 1 valide: le fichier toto.txt contient les lignes suivantes
    ################ debut fichier Commentaire Python

    dict(a=1, b=2) #expression Python qui démarre COLONNE 1 (=> ou pas)

    ################ Autre commentaire, fin fichier

    >>> d = readPyData('toto.txt')
    >>> d
    ... {a=1, b=2}

    :example 2 valide: le fichier toto.txt contient
        une expression à deux valeurs séparées par une ','

    dict(a=1, b=2),0

    >>> d = readPyData('toto.txt')
    >>> d

    :example 3 invalide: le fichier toto.txt contient deux valeurs
        séparées par un '\n'
    dict(a=1, b=2)
    0

    >>> d = readPyData('toto.txt')
    >>> d
    ......
    SyntaxError: invalid syntax

    :example 3: le fichier toto.txt contient une INSTRUCTION Python
    ##Le fichier toto.txt
    x=3
    # fin fichier

    >>> d = readPyData('toto.txt')
    >>> d
        x=3
         ^
    SyntaxError: invalid syntax

    """
    array = asarray;array
    with open(filename,'r') as f :
        lines = []
        for l in f : lines.append(l.strip())#nettoyage des espaces
        source = '\n'.join(lines)
        array = asarray;array
        data = eval(source)
    return data
readData = readPyData

class LecteurData(Lecteur):
    """
    Lecteur generique pour fichiers parametres
    - Pour chaque ligne, tout ce qui est après '#' est ignoré
    - Les lignes vides (ou avec que des espaces) sont ignorées
    La méthode Lire(filename) retourne un dictionnaire des variables lues sur filename.
    Exemple, pour un fichier 'machin.dat' ainsi constitué (les < xxx > ne font pas partie du fichier):
    -----------------------------------------------------
    # commentaire
    # encore commentaire
    <blank line>
    var1 1. 2. 3.                # liste de reels
    var2 2                       # un entier
    var3 'chaine de caracteres'  # une chaine
    var3 "chaine de caracteres"  # une chaine
    var5 <expression>            # toute expression Python évaluable par eval, dans le contexte
    var4 variable illisible      # ...
    <fin de fichier machin.dat>

    L'instruction de lecture est :
    >>> d = LecteurData("machin.dat").lire()
    d est un dictionnaire.

    Pour instancier ces variables, dans un objet 'toto' :
    >>> for key, value in d.iteritems() :
    ...     setattr(toto, key, value)

    C'est tout.
    """
    def __init__(self,filename):
        Lecteur.__init__(self,filename,readlines=True)

    def split(self,line):
        """Partage au premier espace trouvé"""
        line=line.strip()
        comm=line.find('#')
        if comm>=0 :
            return self.split(line[:comm])
        try : b=line.index(' ')
        except ValueError : return '' #pas de blanc, pas de valeur
        return line[:b].strip(),line[b:].strip()

#    def correctiveAction(self, word):
#        """
#        Voir self.value(), en cas de dictionnaires imbriqués, il faut instancier les valeurs du dictionnaire
#        comme des variables locales
#        exemple : on a lu un dictionnaires self.dict['d'] qui vaut {'a':1}
#        et on tombe sur une ligne
#        d1 {'aa':d}
#        word vaut donc "{'aa':d}" qui n'est pas évaluable car d est inconnu dans le contexte,
#        on ne connaît que self.dict['d’] qui vaut {'a':1}
#        on instancie self.dict['d’] comme une variable locale d
#        on réévalue word
#        """

    def value(self,word=""):
        if len(word)==0 :
            return None
        try :
            return eval(word)
        except (NameError) as msg :
            try : #peut-etre des dictionnaires imbriqués
                for key,value in self.dict.items() :
                    if isinstance(value, str) :
                        exec('%s="%s"'%(key,str(value)))
                    else :
                        exec('%s=%s'%(key,str(value)))
                return eval(word)
            except (ValueError,NameError) as msg:
                print(whoami(self)," NameError::valeur illisible : %s"%word,msg, file=sys.stderr)
        except (ValueError) as msg :
            print(whoami(self)," ValueError::valeur illisible : %s"%word,msg, file=sys.stderr)
        return None


    def lire(self):
        """
        Fichier data :
        ===============
        # commentaire
        # encore commentaire

        var1 1. 2. 3.                # liste de reels
        var2 2                       # un entier
        var3 'chaine de caracteres'  # une chaine
        var3 "chaine de caracteres"  # une chaine
        var5 <expression>            # toute expression Python évaluable par eval, dans le contexte
        var4 variable illisible      # ...

        """
        self.dict={}
        for line in self.lines :
            line=line.strip()
            if not line : continue
            try :
                key,value=self.split(line)
            except ValueError :
                continue
            self.dict[key]=self.value(value)

        return self.dict

if __name__=="__main__":
    from tests.inout.testslecteurdata import Tests
    t = Tests()
    if 0 : t.testLecteurData()
    if 1 : t.testReadPyData()
