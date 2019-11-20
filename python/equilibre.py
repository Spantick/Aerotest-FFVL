#!/usr/bin/python
#coding:utf-8
#equilibre 1
from math import sqrt, sin, cos, pi
from scipy.optimize import newton#, brentq, fsolve
from utils.debog import debug, rdebug, stack, sursouligne, mexit
# from pprint import pprint
# from scipy.interpolate import UnivariateSpline
from numpy import (radians, degrees, nan, asarray, hstack, ones_like, linspace,
                   where)
from matplotlib import pyplot as plt
# from aperoconfig import DATA_DIR
# import aperoconfig as AC
# from model.basicobjects.pathbook import PathBook
from utils.utilitaires import maintenant, sapero
# from simulation.aerodynamique.simulaero import allSimulations
from inout.format import fmt
from inout.lecteurdata import readPyData
from path import Path
from numpy.linalg.linalg import norm
from scipy.interpolate.fitpack2 import UnivariateSpline
# from simulation.aerodynamique.simulaero import allSimulations
from pprint import pprint
from utils.torseur import Torseur, Point, Vecteur
from utils import torseur
g, ρ, oswald, ν = 9.81, 1.05, 0.9, 15.6e-6
# from aperoconfig import RUNS_DIR
def reynolds(v,l):
    """nb Reynolds
    :param v:vitesse (m/s)
    :param l:longueur référence (m)
    :param nu:viscosité cinématique (m2/s)
    """
    #suspente : v=10, l dans[1,2]*10^-3 => Re~[600-1200]
    return v*l/ν

def SCxs(L, D, v=10.0, unit=['m','m','m/s']):
    """SCx du suspentage.
    :param L:list L[i]=longueur (m)
    :param D:list D[i]=diametre (m)
    :param v: float, vitesse référence (m/s)
    """
    Cx = 1.0#Cx d'un cylindre de longueur infinie, Re=1000
    return sum([l*d for l,d in zip(L,D)])*Cx

class Polaires(object) :
    """
    Calcul des polaires du parapente 'session' à partir de tous les résultats
    de simulation (les .cmo) disponibles dans le répertoire
        <RUNS_DIR>/<session>/simulation/aerodynamique/
    """
    def __init__(self, **kargs):
#         session = kargs.pop('session') if 'session' in kargs else AC.SESSION
#         self.pb = PathBook('')
#         self.session = session
        self.files = set()
        # debug("Le fichier %s existe : %s"%(self.pb.polaires.name,self.pb.polaires.isfile()))
        self.load(kargs)
        w, alfadegs, alfas, Czs, Cxs, Cms = self.setCollapseAndStall(15, None)
#         debug(alfasdeg=fmt(alfadegs))
#         mexit()
        self.Cz = UnivariateSpline(alfas, Czs, w=w, s=0.0005)
        self.Cx = UnivariateSpline(alfas, Cxs, w=w, s=0.000005)
        self.Cm = UnivariateSpline(alfas, Cms, w=w, s=0.0005)
        self.alfadegs, self.alfas, self.Czs, self.Cxs, self.Cms = alfadegs, alfas, Czs, Cxs, Cms

    def setCollapseAndStall(self, 𝛂s=15, 𝛂c=None):
        # TODO :à travailler décrochage et fermeture
        # on supprime les valeurs au dessus de l'incidence de décrochage 𝛂s~15
        # on rajoute l'incidence de fermeture 𝛂c ~ 0
        if 𝛂c is None and 𝛂s is None :
            return None, self.alfadegs, self.alfas, self.Czs, self.Cxs, self.Cms
        elif 𝛂s is not None : #and 𝛂c is not None:
            𝛂s1 = 𝛂s + 1
            alfadegs = asarray(self.alfadegs)
            w = where(alfadegs <= 𝛂s)
            #On laisse tomber ce qui est au delà de l'incidence de décrochage
            alfadegs = alfadegs[w]
            alfas = self.alfas[w]
            # On ajoute l'incidence 𝛂s+1.0°
            alfadegs = hstack((alfadegs, 𝛂s1))
            alfas = hstack((alfas, radians(𝛂s1)))
            Czs = asarray(self.Czs)[w]
            # pour éviter que la polaire Cz(𝛂) ne soit tordue
            Czs[-1] = Czs[-2]
            # on met la dernière portance à 0.0
            Czs = hstack((Czs, 0))  # Czs[-1]-0.2))

            # La dernière traînée augmente beaucoup
            Cxs = self.Cxs[0:len(Czs)]
            Cxs[-1] *= 1.5
            # le dernier coeff de moment
            Cms = self.Cms[0:len(Czs)]
            Cms[-1] = 0
            # Cms[-1] = -1
            w = 5*ones_like(Czs)
            w[-2] = 0.001
            w[-3] = 0.001
            # w[1] = 1.0
            return w, alfadegs, alfas, Czs, Cxs, Cms
        elif 𝛂c is not None :
            raise NotImplementedError(" Fermeture non prise en compte, vol impossible")
            alfadegs = asarray(self.alfadegs)
            w = where(alfadegs >= 𝛂c)
            alfadegs = alfadegs[w]
            alfas = self.alfas[w]
            alfadegs = hstack((𝛂c, alfadegs))
            alfas = hstack((radians(𝛂c), alfas))
            Czs = asarray(self.Czs)[w]
            Czs = hstack((-0.05, Czs))

            Cxs = self.Cxs[1:len(Czs)]
            Cxs[0] *= 1.5
            Cms = self.Cms[1:len(Czs)]
            Cms[0] = +1
            w = ones_like(Czs)
            w[1] = 1.0/5
            return w, alfadegs, alfas, Czs, Cxs, Cms
            # debug(w=w.shape, Czs=self.Czs.shape, alfas=self.alfas.shape)
            # debug(Czs=self.Czs.tolist())
            # debug(alfas=self.alfas.tolist())
            # exit()
            # print('knots=', self.Cz.get_knots().tolist())
            # print('coefs=', self.Cz.get_coeffs().tolist())

    def loadFromAllSimulations(self, resume=True):
        """On recharge les résultats des simulations :
            - si resume=True on relit uniquement les fichiers résumé,
            - sinon, on load les fichiers cmo et on cree les résumés
            un fichier résumé se lit avec readPyData, c'est un dict constitué ainsi :

            'dims' : dict, les dimensions nbp, nbn, nbpn, nbcells, nbsect sont les nombres de
                points, nervures, points par nervure, panels, sections (ou caissons)

            'convergence' : bool, True si les itérations Cmarc ont convergé.

            'inputdata' : dict, les paramètres passés à cmarc
                 incidence, dérapage, nb timestep etc., cf fichier .cmi

            'glob' : list, pour chaque pas de temps, un dictionnaire qui contient
                les coefficients aerodynamiques globaux, dans le repere aerodynamique,
                calcules par Cmarc.
                NOTE1: si la geometrie est decrite avec une symetrie par rapport au plan Y=0,
                (parametre RSYM=0) les coefficients ne tiennent compte que de la demi voile.
                NOTE2: on suppose que l'assemblage cmarc n'a qu'un seul composant 'Voile'

            'local' : list, pour chaque caisson, un dictionnaire des
                coefficients aerodynamiques locaux (du caisson), calcules par Cmarc,
                exprimes dans le repere aerodynamique, et calcules au DERNIER pas de temps
        """

        allsim = allSimulations(AC.SESSION, resume)
        alfadegs, yawdegs, Czs, Cxs, Cms = [],[],[],[],[]
        areas, S0s, vinfs, b0s, l0s, refms = [],[],[],[],[],[]
        for _, sim in enumerate(allsim[:]):
            print(str(sim))
#             debug('>>>> Simulation %s\n'%cmo.session)
#             sim = SimulationAero(cmo)
            if not sim.hasConverged() : continue
#             lecteur = postaero.lecteur
#             d = coefs = sim.getGlobalCoeffs()['Voile']
#             debug(coefs.keys())
            alfadeg = sim.alfa
            if alfadeg in alfadegs :
                continue
            yawdeg = sim.yaw
            vinf = sim.vinf#Vitesse ref
            S0 = sim.sref#Surface ref
            b0 = sim.sspan#1/2 envergure ref
            l0 = sim.cbar#Corde ref
            #Point de référence pour moment tangage
            refm = sim.rmpx, sim.rmpy, sim.rmpz
#             vinf = d['vinf']#Vitesse ref
#             S0 = d['sref']#Surface ref
#             b0 = d['sspan']#1/2 envergure ref
#             l0 = d['cbar']#Corde ref
#             #Point de référence pour moment tangage
#             refm = d['rmpx'],d['rmpy'],d['rmpz']
            refm
#             debug(entete=d)
    #             debug('entete')
#                 pprint(d)
#             debug('Coefs aero')
#             pprint(F)

#             coefs = F['Voile']
            Cz, Cx, Cm, area = sim.CL, sim.CD, sim.C_m, sim.area*S0
            alfadegs.append(alfadeg)
            yawdegs.append(yawdeg)
            Czs.append(Cz)
            Cxs.append(Cx)
            Cms.append(Cm)
            vinfs.append(vinf)
            areas.append(area)
            S0s.append(S0)
            b0s.append(b0)
            l0s.append(l0)
            refms.append(refm)

        alfas = radians(alfadegs)
        yaws = radians(yawdegs)
        self.alfadegs = tuple(alfadegs)
        self.yawdegs = tuple(yawdegs)
        T = list(zip(*sorted(zip(alfadegs, alfas, yawdegs, yaws, Czs, Cxs, Cms,
                            areas, S0s, vinfs, b0s, l0s))))
        (self.alfadegs, self.alfas, self.yawdegs, self.yaws, self.Czs, self.Cxs, self.Cms,
        self.areas, self.S0s, self.vinfs, self.b0s, self.l0s) = T

        self.verifications()

    def verifications(self):
        '''Verification : les surface, corde et 1/2 envergure
        des différentes simulations doivent être identiques'''
        if len(set(self.b0s))==1 :
            self.b0 = self.b0s[0]
            del self.b0s
        else :
            msg = "Les envergures de reference des differentes simulations ne sont pas identiques : \n  (alfa, b0)=%s"%list(zip(self.alfas, self.b0s))
            raise ValueError(msg)

        if len(set(self.S0s))==1 :
            self.S0 = self.S0s[0]
            del self.S0s
        else :
            msg = "Les surfaces de reference des differentes simulations ne sont pas identiques : \n  (alfa, S0)=%s"%list(zip(self.alfas, self.S0s))
            raise ValueError(msg)

        if len(set(self.l0s))==1 :
            self.l0 = self.l0s[0]
            del self.l0s
        else :
            msg = "Les cordes de reference des differentes simulations ne sont pas identiques : \n  (alfa, l0)=%s"%list(zip(self.alfas, self.l0s))
            raise ValueError(msg)

        if len(set(self.areas))==1 :
            self.area = self.areas[0]
            del self.areas
        else :
            msg = "Les aires de reference des differentes simulations ne sont pas identiques : \n  (alfa, area)=%s"%list(zip(self.alfas, self.areas))
            raise ValueError(msg)

        if len(set(self.vinfs))==1 :
            self.vinf = self.vinfs[0]
            del self.vinfs
        else :
            msg = "Les vitesses de reference des differentes simulations ne sont pas identiques : \n  (alfa, vinf)=%s"%list(zip(self.alfas, self.vinfs))
            raise ValueError(msg)

        if len(set(self.yawdegs))==1 :
            self.yawdeg = self.yawdegs[0]
            del self.yawdegs
        else :
            msg = "Les derapages des differentes simulations ne sont pas identiques : \n  (alfas, yawdegs)=%s"%list(zip(self.alfas, self.yawdegs))
            raise ValueError(msg)

    @property
    def info(self):
        infos = []
        for attr in ('alfadegs', 'Czs', 'Cxs', 'Cms',
                     'yawdeg', 'area', 'vinf', 'S0', 'b0', 'l0') :
            value = getattr(self, attr)
            infos.append('%10s : %s'%(attr, value))
        return infos

    def __str__(self):
        return '\n'.join(self.info)

    def plot(self):
        adeg = linspace(min(self.alfadegs), max(self.alfadegs), 100)
        arad = linspace(min(self.alfas), max(self.alfas), 100)
        Cxs, Czs = self.Cx(arad), self.Cz(arad)
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222, sharex=ax1)
        ax3 = fig.add_subplot(223, sharex=ax1)
        ax4 = fig.add_subplot(224)
        ax4.set_xlim(left=0.0, right=max(Cxs))
        ax4.set_ylim(bottom=0.0, top=max(Czs))
        ax1.plot(adeg, Czs, 'r-', label='$C_z$')
        ax1.plot(self.alfadegs, self.Czs, 'b.', label='data')
        ax1.legend()
        ax1.set_xlabel('$ \\alpha $ (°)')
        ax2.plot(adeg, Cxs, 'b-', label='$C_x$')
        ax2.plot(self.alfadegs, self.Cxs, 'b.', label='data')
        ax2.legend()
        ax2.set_xlabel('$ \\alpha $ (°)')
        ax3.plot(adeg, self.Cm(arad), 'g-', label='$C_m$')
        ax3.plot(self.alfadegs, self.Cms, 'b.', label='data')
        ax3.legend()
        ax3.set_xlabel('$ \\alpha $ (°)')

        ax4.plot(Cxs, Czs, 'k-', label='Polaire Eiffel')
        ax4.legend()
        ax4.set_xlabel('$ C_x $')
        ax4.set_ylabel('$ C_z $')
    #     plt.plot()
    #     plt.legend(['Linear', 'InterpolatedUnivariateSpline', 'True'])
    #     plt.axis([-0.05, 6.33, -1.05, 1.05])
        fig.suptitle('%s : Coefficients aérodynamiques\n(voile seule)'%AC.SESSION)
        plt.show()

    def toDump(self):
        return dict(alfadegs=self.alfadegs,
                    Czs=self.Czs, Cxs=self.Cxs, Cms=self.Cms,
                    area=self.area, yawdeg=self.yawdeg, vinf=self.vinf,
                    S0=self.S0, b0=self.b0, l0=self.l0)

    def dump(self, filename=None):
        if filename is None :
            filename=self.pb.aerodir/AC.SESSION+'.pol'
#             debug(filename=filename)
        with open(filename,'w') as f :
            pprint(self.toDump(), f, width=120)
        self.files.add(filename)

    def loadFromFile(self, filename=None):
        if filename is None :
            filename=self.pb.polairesfile
        return readPyData(filename)
#         self.load(D)

    def load(self, kargs):
        for key, value in kargs.items() :
            setattr(self, key, value)
        try :
            self.alfas = radians(self.alfadegs)
        except AttributeError :
            pass
#         rdebug('fin load')
        return 0

class Equilibriste(object):
    EPS = 1.0e-5
    """
        Manipulation des équations de l'equilibre longitudinale,
        a partir des donnees Cz,Cx,Cm... de CMARC, ou autre
        Les attributs modifiables sont :
        d, b, S, l, s, xGv, zGv, mp, mv, c, zs, Cx0, Cxp, Cxs,
        les autres sont calculés.
        Les coordonées sont exprimées dans le repère aérodynamique (0,(i,j,k))
        dont l'origine est O = le BA du caisson central,
        et vu par un pilote : i vers l'arrière, j vers la droite et k vers le haut
    """
    VARIABLES = dict(
        b = 'Demi envergure à plat (m).',
        S = 'Surface voilure à plat (m2).',
        l = 'Corde centrale (m).',
        s = 'Longueur du cône de suspentage (m).',
        xGv ='Position en z du centre de gravité de la voile, xGv>0. (m)',
        zGv ='Position en x du centre de gravité de la voile, zGv<0. (m)',
        zs = '''Position en z du point d'application de la trainée suspente (m). zs<0.''',
        mp = 'Masse pilote (kg)',
        mv = 'Masse voile seule (kg). NON compris la masse d\'air emprisonné.',
        ma = 'Masse masse d\'air emprisonné.',
        c = '''Calage en %de corde centrale : c = AQ/b, où
                    Q est la projection orthogonale du pilote sur la corde centrale,
                    A est le bord d'attaque de la corde centrale.''',
        d = 'Pourcentage de frein ∈ [0.0, 100.0] => non implémenté',
        Cx0 = 'Coefficient de traînée voilure seule.',
        Cxp = 'Coefficient de traînée pilote.',
        Cxs = 'Coefficient de traînée suspentage.')

    def __init__(self, **kargs):
        """
        Calcul approximatif de l'équilibre en vol droit, à vitesse constante,
        de la voile <session>, à partir des data (globales) calculées par Cmarc
        (ou un autre simulateur aérodynamique).
        Les coordonnées sont exprimées dans le repère voile,
        d'origine le bord d'attaque du caisson central, de vecteurs iv,jv,kv,
        - iv porté par la corde centrale, vers l'arrière,
        - jv vers la gauche,
        - kv vertical, vers le haut.
        Sauf précision, les unités sont en S.I. (Mètre, Kilogramme, Seconde).
        ########################################################################
        ####  ICI, ON N'UTILISE PAS LES VARIABLES GLOBALES DE APEROCONFIG   ####
        ########################################################################

        kargs['pdata'] : variables calculés par self.parapente.data2d()
        --------------------------------------------------------------
        :param b : float, demi-envergure vraie.
        :param S : float, surface vraie.
        :param l : float, corde centrale vraie.
        :param s : longueur du cône de suspentage.
            P et Q étant la position du pilote et la projection orthogonlre du
            pilote sur la corde centrale, s=dist(P,Q).
        :param xGv,zGv : float, float, position du centre de gravité de la voile
            xGv>0, zGv<0.

        Variables par défaut dans <DATA_DIR>/preferences/equilibre2d.dat
        ----------------------------------------------------------------
                puis écrasées <session>.equilibre2d.dat
                ---------------------------------------
        kargs['userdata']
        :param mp: float, masse pilote.
        :param d: pourcentage de frein ∈ [0.0, 100.0] => non implémenté
        :param Cxp : coefficients de traînée pilote.
                [TODO : a faire calculer par PARAPENTE]
                ---------------------------------------
        :param mv:float, masse voile seule, NON compris la masse 'ma' d'air emprisonné,
            dont le poids est intégralement compensé par la force d'Archimède.
            Lorsque l'on s'intéressera à la DYNAMIQUE du parapente, il faudra prendre
            'ma' en compte
        :param ma:float, masse d'air emprisonné dans la voile [non utilisé actuellement]
        :param c: float, calage, i.e. position de Q en %de corde centrale
            i.e. AQ/b, A=bord d'attaque au centre de la voile.
        :param zs: float, altitude du point d'application de la trainée suspente.
            zs<0.
        :param Cxs : coefficients de traînée suspentage.
        :param Cx0: float, coefficient de traînée voile (frottement).

        Fonctions et param. retournées par Polaires
        -------------------------------------------
        :param b0 : float, demi-envergure de référence
            (=SSPAN utilisée par Cmarc => for rolling and yawing moments).
        :param S0 : float, surface de référence
            (=SREF utilisée par Cmarc => for force and moment coefficients).
        :param l0 : float, corde de référence
            (=CBAR utilisée par Cmarc => for pitching moment).
        :param v0 : float, vitesse de référence
            (=VINF utilisée par Cmarc => free-stream velocity).
        :param Cx : fonction (quadratique ou spline) de l'incidence, le
            coefficient de trainée voile fourni par la simulation Cmarc.
        :param Cz : fonction (quadratique, linéaire ou spline) de l'incidence,
            le coefficient de  portance, fourni par la simulation Cmarc.
        :param Cm : fonction (quadratique, linéaire ou spline) de l'incidence,
            le coefficient de  moment par rapport au BA, fourni par la
            simulation Cmarc.
        NB : les trois dernières fonctions sont des fonctions f(alfa) qui
        prennent en argument 'alfa' l'incidence
        Cmarc fournit les trois fonctions aérodynamiques de 'alfa' pour la voilure
        seule alfa -> Cz(alfa), Cx(alfa), Cm(alfa)
        :TODO: l'influence des freins 'd' doit être évaluée par ailleurs...
        NB (mars 2019): pour obtenir la position de la voile en soufflerie, à vitesse V
        donnée, il suffit d'utiliser le simulateur itérativement avec une charge
        (poids pilote) variable, jusqu'a obtenir une vitesse de vol = V
        La méthode self.alfaV(v0) fait le job
        [Fixed] => NB (mars 2019) ATTENTION,
            cmarc donne le moment des forces par rapport au point
                BINP9{..., 'RMPX': 0.0, 'RMPY':0.0, 'RMPZ':0.0}
            valeur par defaut qu'on retrouve dans tous les cmi et cmo
        """

        #valeurs par defaut
        dump = readPyData(DATA_DIR/'preferences'/'equilibre2d.dict')
        dump.update(kargs['userdata'])
        dump.update(kargs['pdata'])


        #Les valeurs de b0, S0,.. et les fonctions Cx, Cz..
        #dans Polaires
        pol = kargs['polaires']
        dump.update(
            CxvCMARC = lambda alfa:pol.Cx(alfa%pi),
            Cmv      = lambda alfa:pol.Cm(alfa%pi),
            Cz       = lambda alfa:pol.Cz(alfa%pi),
            b0=pol.b0, l0=pol.l0, vinf=pol.vinf, S0=pol.S0
            )

        #On écrase tout par les paramètres d'entrée
        dump.update(kargs)

        #on se débarasse du paramètre session qui est en @property
        if 'session' in dump :
            dump.pop('session')

        #On a toutes les data dans dump, on peut instancier les attributs
        for key, value in dump.items():
            setattr(self, key, value)
        
        # ne pas oublier, pour calculer les différentes constantes utiles.
        self.load()
    
    def Cxv(self, alfa) :
        """fonction Cx voile"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return self.CxvCMARC(alfa) + self.Cx0

    def CxT(self, alfa) :
        """fonction alfa ⟼ Cx(𝛼) = coef traînée total"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))


        return self.Cxv(alfa) + self.Cxs + self.Cxp
    # CxT = Cx

    def CmT(self, alfa):
        """fonction Coef de moment global(BA), aile + pilote + suspentes"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        Cmv = self.Cmv(alfa)
        Cz = self.Cz(alfa)
        θ = self.θ(alfa)
        return self.mx + θ*self.mz + (Cmv+self.Cmz-alfa*self.Cmx)/Cz

    def γ(self, alfa) :
        """fonction alfa ⟼ γ(alfa) angle de plané"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return self.CxT(alfa)/self.Cz(alfa)

    def θ(self, alfa):
        """fonction alfa ⟼ θ(alfa) assiette"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return self.γ(alfa)-alfa

    def CP(self, alfa):
        """Centre de poussee de la RFA en % de corde
        = intersection de la droite d'action de la RFA avec la corde."""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return -self.Cmv(alfa)/(alfa*self.Cxv(alfa)+self.Cz(alfa))

    def normalise(self,alfa):
        """alfa ramené dans [0,pi]"""
        return alfa%pi

    def V(self, alfa) :
        """ Vitesse air"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        value = (2*self.mg)/(ρ*self.Cz(alfa)*self.S0)
        if value < 0 :
            value = -sqrt(-value)
#             debug("Vol dos Cz<0, alfa = %.3g°, V = %.3g km/h"%(degrees(alfa),-3.6*value))
            return value
            # exit()
        else :
            return sqrt(value)

    def typeDeVol(self, alfa):
        """le vol avec l'incidence alfa peut être vol normal, dos, arrière ou dos-arrière"""
        alfa = alfa%pi
#         debug(alfa_deg=degrees(alfa))
        if 0<=alfa<=pi/2 :
            return 0,'vol normal'
        elif pi/2<=alfa<=pi :
            return 1,'vol marche arrière'
        elif -pi/2<=alfa<=0 :
            return 2,'vol dos'
        elif -pi<=alfa<=-pi/2 :
            return 3,'vol dos + marche arrière'
        else :
            return 4,'vol ?'

    def Vz(self,alfa) :
        """taux de chute"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return self.V(alfa)*sin(self.γ(alfa))

    def qS(self, alfa) :
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return (self.mg)/self.Cz(alfa)

    def resultante(self, alfa):
        """resultante en Newtons"""
#         Cx, Cz = self.aresultante(alfa)
#         return Cx*self.mg,Cz*self.mg
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        γ = self.γ(alfa)
        qS = self.qS(alfa)
#         qSm = self.qS(alfa)/self.mg
#         mg = self.mg
        return -self.mg*γ+qS*self.CxT(alfa), -self.mg+qS*self.Cz(alfa)
        return -self.mg*sin(γ)+qS*self.CxT(alfa), -self.mg*cos(γ)+qS*self.Cz(alfa)

    def aresultante(self, alfa):
        """resultante a-dimensionnée"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        γ = self.γ(alfa)
        qSm = self.qS(alfa)/self.mg
        return -sin(γ)+qSm*self.CxT(alfa), -cos(γ)+qSm*self.Cz(alfa)
        return -γ+qSm*self.CxT(alfa), -qSm*self.Cz(alfa)

    def finesse(self, alfa) :
        """finesse à incidence alfa donnée"""
        if alfa is not nan and alfa != alfa%pi : debug('alfa != alfa%%pi : %s, %s'%(alfa,alfa%pi))

        return self.Cz(alfa)/self.CxT(alfa)

    def alfa0(self) :
        """le alfa de portance nulle"""
        return newton(lambda alfa:self.Cz(alfa), radians(6.0))

    def objectif(self,alfa):
        """
        La fonction à annuler pour trouver l'équilibre
        Peut être que return self.CmT(alfa) suffit...
        """
        alfa = alfa%pi
        cx, cz = self.resultante(alfa)
        return self.CmT(alfa)**2 + cx**2 + cz**2
        return cx**2 +cz**2

#     def objectif1(self,alfa):
#         """
#         La fonction à annuler pour trouver l'équilibre.
#         On rajoute une force extérieur (acceleration du camion) aux equations
#         """
#         alfa = alfa%pi
#
#         fx, fz = self.resultante(alfa)+self.addForce(gamma)
#         return (self.CmT(alfa)+self.addMoment(alfa))**2 + fx**2 + fz**2
#         return fx**2 +fz**2

    def equilibre(self,**kargs) :
        """
        :param kargs: les variables à modifier, parmi self.VARIABLES
        :return:float, l'incidence alfa° telle que CmT(alfa)=0
        >>> P = Parapente('diamir')
        >>> eq = Equilibriste(P)
        >>> alfa = eq.equilibre(c=29.0, mp=100.0, S=28.0, d=10.0)
        # alfa sera l'incidence d'équilibre pour le parapente P,
        # avec un calage c, un pilote de masse mp, une surface S, et 10% de freins
        """
#         objectif = self.CmT
        if kargs : self.load(kargs)
        #disp=False pour que newton ne lève pas une exception
        alfa, res = newton(self.objectif, radians(6.0), tol=self.EPS, full_output=True, disp=False)
        if res.converged :
            return alfa%pi
        else :
            alfa = res.root
            res = dict(root=fmt(degrees(res.root)),iterations=res.iterations,
                       function_calls=res.function_calls,
                       converged=res.converged, cause=res.flag)
            arn = norm(self.resultante(alfa%pi))
            ar  = norm(self.resultante(alfa))
            cm = self.CmT(alfa)
            debug('Non CV', resultante=(arn, ar), CmT=cm, res=res)
            return nan
#             if abs(ar)<self.EPS :
#                 return alfa%pi
#             debug('Non CV', resultante=(fmt(rn,6), fmt(r,6)))
    alfaEq = equilibre

    def toDump(self):
        keys = sorted(list(self.VARIABLES))
        return dict([(key, getattr(self,key)) for key in keys])

    def load(self, kargs={}):
        """
         :param kargs: les variables à modifier, parmi self.VARIABLES
         >>> P = Parapente('diamir')
         >>> eq = Equilibriste(P)
         >>> eq.load(c=29.0, mp=100.0, S=28.0, d=10.0)
         # sont modifiés calage c, masse pilote mp, surface S, et freins 10%
         """
        assert set(kargs.keys()).issubset(self.VARIABLES), \
            """
    Une des variables parmi %s est inconnue.
        Les variables autorisées sont : %s,
            """%(sorted(kargs.keys()), sorted(list(self.VARIABLES)))
        for key, value in kargs.items() :
            setattr(self, key, value)

        #constantes
        self.xc = self.c*self.l
        self.zp = -self.s
        self.m    = self.mp + self.mv
        self.Cx0T = self.Cx0 + self.Cxp + self.Cxs
        self.Cmx = (self.xc/self.l0)*(self.Cxp + self.Cxs)
        self.Cmz = (self.zs/self.l0)*self.Cxs - (self.s/self.l0)*self.Cxp
        self.mx  = (self.mp/self.m)*(self.xc/self.l0)\
                    + (self.mv/self.m)*(self.xGv/self.l0)
        self.mz  = (self.mp/self.m)*(self.s/self.l0)\
                    - (self.mv/self.m)*(self.zGv/self.l0)
        self.mg  = self.m*g
    def alfaVOld(self, v0):
        """Retourne le alfa tel que V(alfa) = v0. Ne correspond pas à l'équilibre."""
        try :
            return newton(lambda alfa : self.V(alfa)-v0, radians(6.0))
        except RuntimeError as msg : return nan

    def alfaV(self, v0):
        """Retourne le alfa tel que V(alfa) = v0. Ne correspond pas à l'équilibre."""
#         return newton(lambda alfa : self.V(alfa)-v0, radians(6.0), full_output=True, disp=False)
#         except RuntimeError as msg : return nan
        alfa, res = newton(lambda alfa : self.V(alfa)-v0, radians(6.0), tol=self.EPS, full_output=True, disp=False)
        if res.converged :
            return alfa%pi
        else :
            alfa = res.root
            res = dict(root=fmt(degrees(res.root)),iterations=res.iterations,
                       function_calls=res.function_calls,
                       converged=res.converged, cause=res.flag)
            rn = norm(self.resultante(self.normalise(alfa)))
            r  = norm(self.resultante(alfa))
            cm = fmt(self.CmT(alfa))
            debug('alfaV : non CV')
            print("alfa_normalisé = %s,"%fmt(degrees(alfa%pi)), 'resultante = (%s,%s),'%fmt((rn, r)), 'CmT =', cm)
            print('résultats newton = ', fmt(res))
            print()
            return nan
#             if abs(r)<self.EPS :
#                 #NOOOON, si la resultante est nulle et le moment non nul,
#                 #alors c'est un couple, on n'est pas à l'équilibre
#                 #return alfa%pi

    def cV(self, v0):
        """Retourne le calage c pour atteindre la vitesse v0"""
        return newton(lambda alfa : self.V(alfa)-v0, radians(6.0))

    def VEq(self, **kargs):
        """
        :param kargs : dict(VARIABLES)
        :returns la vitesse à l'équilibre"""
        return self.V(self.alfaEq(**kargs))

    def 𝞥max(self):
        """:returns : la finesse max accessible, et le calage correspondant.
            masse pilote fixéé"""
        raise NotImplementedError("TODO")

    def write(self, filename, **kargs):
        """
        :param freins:tuple, les % de freins pour lesquels on calcule l'équilibre, ∈[0.0, 100.0]
        :type filename: fichier ou str
        """

        if isinstance(filename, str) :
            fout = open(filename,'w')
        else : #fichier ouvert. p. ex. sys.stdout
            fout = filename
        dump = fmt(self.toDump())
        print(sapero(), maintenant(), "\n\n", file=fout)
        with fout :
            titre=sursouligne("  Lexique et valeurs des paramètres  ","=")
            print(titre,file=fout)
            for key in sorted(self.VARIABLES) :
                print("%5s = %-6s = %s"%(key, dump[key], self.VARIABLES[key]), file=fout)
            print("\nN.B. la valeur 'nan' (Not A Number), pour signifier que le calcul n\'a pas abouti (plantage, non convergence..)", file=fout)
            print("""\nN.B. suivant la valeur de l'incidence, le vol peut être d'un des types suivants (dans les 2 ou 3 derniers cas... à utiliser avec circonspection !):
            .   > vol normal,
            *   > vol marche arrière,
            **  > vol dos,
            *** > vol dos + marche arrière.""", file=fout)

            self._writeEquilibre(fout)
            for key, values in kargs.items() :
                if key=='mp' :
                    self._writeVariationMassePilote(values, fout)
                if key=='c' :
#                     rdebug(values)
                    self._writeVariationCalage(values, fout)

    def _writeEquilibre(self, fout):
        eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu
        titre = sursouligne('  Equilibre  ',"=")
        print('\n\n',file=fout)
        print(titre,file=fout)
        if 0 : print("Paramètres:\n%s"%fmt(self.toDump()),file=fout)
        CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
        DELTAs, THETAs, ALPHAs = [], [], []
        alfa = self.equilibre()
        ALPHAs.append(alfa*180/pi)
        FINs.append(self.finesse(alfa))
        Vs.append(self.V(alfa)*3.6)
        VZs.append(self.Vz(alfa))
        CPs.append(100*self.CP(alfa))
        THETAs.append(self.θ(alfa)*180/pi)
        r = self.resultante(alfa)
        REs.append(abs(r[0]) + abs(r[1]))# + abs(self.CmT(alfa)))
#         CMTs.append(self.qS(alfa)*self.l0*self.CmT(alfa))
        CMTs.append(self.qS(alfa)*self.l0*self.CmT(alfa))
        CZs.append(self.Cz(alfa))
        CXVs.append(self.Cxv(alfa))
        CMVs.append(self.Cmv(alfa))
        # debug(titre="Je ne sais plus d'ou vient la formule pour simuler les freins !!")
        print(" Frein (%)           |", " ".join(['%-7.1f'%c for c in DELTAs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
        print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
        print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
        print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
        print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
        print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
        print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
        print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
        print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
        self.load(eldump0)

    def _writeVariationMassePilote(self, Ms, fout):
        eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu
        eldump0['mp'] = '***'
        print("\n\n",file=fout)
        titre = sursouligne("  Variation masse pilote  ","=")
        print(titre,file=fout)
        if 0 : print("paramètres :\n", fmt(eldump0), file=fout)
        eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu

        CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
        THETAs, ALPHAs, vols = [], [], []
        for c in Ms :
            try :
                alfa = self.alfaEq(mp=c)
            except RuntimeError as msg :
                debug(msg, calage=c)
                alfa = nan
            idx, _ = self.typeDeVol(alfa)
            vols.append('.' if idx==0 else idx*"*")
            ALPHAs.append(alfa*180/pi)
            FINs.append(self.finesse(alfa))
            Vs.append(self.V(alfa)*3.6)
            VZs.append(self.Vz(alfa))
            CPs.append(100*self.CP(alfa))
            THETAs.append(self.θ(alfa)*180/pi)
            r = self.resultante(alfa)
            REs.append(abs(r[0]) + abs(r[1]))# + abs(self.CmT(alfa)))
            CMTs.append(self.CmT(alfa))
            CZs.append(self.Cz(alfa))
            CXVs.append(self.Cxv(alfa))
            CMVs.append(self.Cmv(alfa))
        rdebug(vols)
        print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
        print(" masse pilote (kg)   |", " ".join(['%-7.2f'%c for c in Ms]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
#         if not any(vols) :
        print(" Type de vol         |", " ".join(['%-7.2s'%c for c in vols]),file=fout)
        print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
        print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
        print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
        print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
        print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
        print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
        print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
        print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
        print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
        self.load(eldump0)

    def _writeVariationCalage(self, Cs, fout):
        eldump0=self.toDump()
        eldump0['c'] = '###'
        titre=sursouligne("  Variation calage  ","=")
        print("\n\n",file=fout)
        print(titre,file=fout)
        if 0 : print("paramètres:\n", fmt(eldump0), file=fout)
        eldump0=self.toDump()
        CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
        THETAs, ALPHAs, typevols = [], [], []
        for c in Cs :
            try :
                alfa = self.alfaEq(c=c)
            except RuntimeError as msg :
                debug(msg, calage=c)
                alfa = nan
                # continue
            if alfa < 0 :
                alfa = nan
                alfa = alfa%pi
                debug ('Fermeture, vol impossible')
            idx, _ = self.typeDeVol(alfa)
            typevols.append('  .' if idx==0 else ' '+idx*"*")
            ALPHAs.append(alfa*180/pi)
            FINs.append(self.finesse(alfa))
            Vs.append(self.V(alfa)*3.6)
            VZs.append(self.Vz(alfa))
            CPs.append(100*self.CP(alfa))
            THETAs.append(self.θ(alfa)*180/pi)
            r = self.resultante(alfa)
            REs.append(abs(r[0]) + abs(r[1]))# + abs(self.CmT(alfa)))
            CMTs.append(self.CmT(alfa))
            CZs.append(self.Cz(alfa))
            CXVs.append(self.Cxv(alfa))
            CMVs.append(self.Cmv(alfa))
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Calage (% de corde) |", " ".join(['%-7.1f'%(100*c) for c in Cs]),file=fout)
        print(" Type de vol         |" + " ".join(['%-7s'%c for c in typevols]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
        print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
        print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
        print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
        print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
        print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
        print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
        print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
        print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
        print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
        print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)

        #On rétablit dans sa config originale
        self.load(eldump0)

    def plot(self, **kargs):
        """TODO : tracés pour variation calage, longueur cône suspentage, CxV, Cxs, Cxp,..."""
#         csfont = {'fontname':'Comic Sans MS'}
        hfont = {'fontname':'Liberation Serif'}

        if not kargs :
            kargs = readPyData(Path(DATA_DIR)/'preferences'/'aero2d.dict')

        if 'mp' in kargs :
            Vs, VZs, CPs, FINs, ALPHAs = [], [], [], [], []
            Mp = kargs['mp']
            for m in Mp :
                alfa = self.alfaEq(mp=m)
                if alfa < 0 :
#                     alfa = self.normalise(alfa)
                    debug ('Fermeture, vol impossible; alfa=%s'%fmt(degrees(alfa)))
                    alfa = nan
                    alfa = alfa%pi
#                 idx, _ = self.typeDeVol(alfa)
#                 typevols.append('  .' if idx==0 else ' '+idx*"*")
                ALPHAs.append(alfa*180/pi)
                FINs.append(self.finesse(alfa))
                Vs.append(self.V(alfa)*3.6)
                VZs.append(self.Vz(alfa))
                CPs.append(100*self.CP(alfa))
                r = self.resultante(alfa)
#                 REs.append(abs(r[0]) + abs(r[1]) + abs(self.CmT(alfa)))
#                 CMTs.append(self.CmT(alfa))
#                 CZs.append(self.Cz(alfa))
#                 CXVs.append(self.Cxv(alfa))
#                 CMVs.append(self.Cmv(alfa))

#         adeg = linspace(min(self.alfadegs), max(self.alfadegs), 100)
#         arad = linspace(min(self.alfas), max(self.alfas), 100)
#         Cxs, Czs = self.Cx(arad), self.Cz(arad)
        debug(finesses=fmt(FINs))
        debug(Vz=fmt(VZs))
        fig = plt.figure()
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222, sharex=ax1)
        ax3 = fig.add_subplot(223, sharex=ax1)
        ax4 = fig.add_subplot(224, sharex=ax1)
#         ax4.set_xlim(left=0.0, right=max(Cxs))
#         ax4.set_ylim(bottom=0.0, top=max(Czs))
        ax1.plot(Mp, FINs, 'b.-')#, label='$ \\phi $')
#         ax1.plot(self.alfadegs, self.Czs, 'b.', label='data')
#         ax1.legend()
        ax1.set_xlabel('$ m_p (kg) $')
        ax1.set_ylabel('$ \\phi $')
        ax1.set_title('Finesses',**hfont)

        ax2.plot(Mp, Vs, 'b.-')#, label='$V$')
        ax2.set_xlabel('$ m_p (kg) $ ')
        ax2.set_ylabel('$ V (km/h) $')
        ax2.set_title('Vitesse sur trajectoire (km/h)',**hfont)
#         ax2.legend()
        ax3.plot(Mp, VZs, 'b.-')#, label='data')
        ax3.set_ylabel('$ V_z (m/s) $')
        ax3.set_title('Taux de chute (m/s)',**hfont)
        if 0 :
            ax4.plot(Mp, ALPHAs, 'b.-')#, label='data')
            ax4.set_ylabel('$ \\alpha (°) $')
            ax4.set_title('Incidences',**hfont)
        if 1 :
            ax4.plot(Mp, CPs, 'b.-')#, label='data')
            ax4.set_ylabel('$ CP $ ')
            ax4.set_title('Centre de poussée % de corde',**hfont)
        fig.suptitle('Session %s : variation de la masse pilote $m_p$\n'%AC.SESSION, **hfont)
        plt.tight_layout()
        plt.show()

# class Equilibriste0Quimarche(object):
#     EPS = 1.0e-5
#     """
#         Manipulation des équations de l'equilibre longitudinale,
#         a partir des donnees Cz,Cx,Cm... de CMARC, ou autre
#         Les attributs modifiables sont :
#         d, b, S, l, s, xGv, zGv, mp, mv, c, zs, Cx0, Cxp, Cxs,
#         les autres sont calculés.
#         Les coordonées sont exprimées dans le repère aérodynamique (0,(i,j,k))
#         dont l'origine est O = le BA du caisson central,
#         et vu par un pilote : i vers l'arrière, j vers la droite et k vers le haut
#     """
#     VARIABLES = set(['d', 'b', 'S', 'l', 's', 'xGv', 'zGv', 'mp', 'mv', 'ma', 'c', 'zs', 'Cx0', 'Cxp', 'Cxs',])
#     VARIABLES = dict(
#         b = 'Demi envergure à plat (m).',
#         S = 'Surface voilure à plat (m2).',
#         l = 'Corde centrale (m).',
#         s = 'Longueur du cône de suspentage (m).',
#         xGv ='Position en z du centre de gravité de la voile, xGv>0. (m)',
#         zGv ='Position en x du centre de gravité de la voile, zGv<0. (m)',
#         zs = '''Position en z du point d'application de la trainée suspente (m). zs<0.''',
#         mp = 'Masse pilote (kg)',
#         mv = 'Masse voile seule (kg). NON compris la masse d\'air emprisonné.',
#         ma = 'Masse masse d\'air emprisonné.',
#         c = '''Calage en %de corde centrale : c = AQ/b, où
#                     Q est la projection orthogonale du pilote sur la corde centrale,
#                     A est le bord d'attaque de la corde centrale.''',
#         d = 'Pourcentage de frein ∈ [0.0, 100.0] => non implémenté',
#         Cx0 = 'Coefficient de traînée voilure seule.',
#         Cxp = 'Coefficient de traînée pilote.',
#         Cxs = 'Coefficient de traînée suspentage.')
#
#     def __init__(self, **kargs):
#         """
#         Calcul approximatif de l'équilibre en vol droit, à vitesse constante,
#         de la voile <session>, à partir des data (globales) calculées par Cmarc
#         (ou un autre simulateur aérodynamique).
#         Les coordonnées sont exprimées dans le repère voile,
#         d'origine le bord d'attaque du caisson central, de vecteurs iv,jv,kv,
#         - iv porté par la corde centrale, vers l'arrière,
#         - jv vers la gauche,
#         - kv vertical, vers le haut.
#         Sauf précision, les unités sont en S.I. (Mètre, Kilogramme, Seconde).
#         Variables globales, dans aperoconfig
#         ------------------------------------
#         :param session : str, nom de la voile
#             (new avril 2019 deux variables globales dans aperoconfig(=ac) :
#                 => session=AC.SESSION, et
#                 => parapente=AC.PARAPENTE).
#         :param parapente : Parapente => self.parapente
#             Si 'parapente' est absent, on utilise AC.PARAPENTE
#         :param session : str, le nom de la session. Inutilisé, on récupère celui de self.parapente
#
#         Variables calculés par self.parapente.data2d()
#         --------------------------------------------
#         :param b : float, demi-envergure vraie.
#         :param S : float, surface vraie.
#         :param l : float, corde centrale vraie.
#         :param s : longueur du cône de suspentage.
#             P et Q étant la position du pilote et la projection orthogonlre du
#             pilote sur la corde centrale, s=dist(P,Q).
#         :param xGv,zGv : float, float, position du centre de gravité de la voile
#             xGv>0, zGv<0.
#
#         Variables par défaut dans <DATA_DIR>/preferences/equilibre2d.dat
#         ----------------------------------------------------------------
#                 puis écrasées <session>.equilibre2d.dat
#                 ---------------------------------------
#         :param mp: float, masse pilote.
#         :param d: pourcentage de frein ∈ [0.0, 100.0] => non implémenté
#         :param Cxp : coefficients de traînée pilote.
#                 [TODO : a faire calculer par PARAPENTE]
#                 ---------------------------------------
#         :param mv:float, masse voile seule, NON compris la masse 'ma' d'air emprisonné,
#             dont le poids est intégralement compensé par la force d'Archimède.
#             Lorsque l'on s'intéressera à la DYNAMIQUE du parapente, il faudra prendre
#             'ma' en compte
#         :param ma:float, masse d'air emprisonné dans la voile [non utilisé actuellement]
#         :param c: float, calage, i.e. position de Q en %de corde centrale
#             i.e. AQ/b, A=bord d'attaque au centre de la voile.
#         :param zs: float, altitude du point d'application de la trainée suspente.
#             zs<0.
#         :param Cxs : coefficients de traînée suspentage.
#         :param Cx0: float, coefficient de traînée voile (frottement).
#
#         Fonctions et param. retournées par Polaires
#         -------------------------------------------
#         :param b0 : float, demi-envergure de référence (utilisée par Cmarc).
#         :param S0 : float, surface de référence (utilisée par Cmarc).
#         :param l0 : float, corde de référence (utilisée par Cmarc).
#         :param Cx : fonction (quadratique ou spline) de l'incidence, le
#             coefficient de trainée voile fourni par la simulation Cmarc.
#         :param Cz : fonction (quadratique, linéaire ou spline) de l'incidence,
#             le coefficient de  portance, fourni par la simulation Cmarc.
#         :param Cm : fonction (quadratique, linéaire ou spline) de l'incidence,
#             le coefficient de  moment par rapport au BA, fourni par la
#             simulation Cmarc.
#         NB : les trois dernières fonctions sont des fonctions f(alfa) qui
#         prennent en argument 'alfa' l'incidence
#         Cmarc fournit les trois fonctions aérodynamiques de 'alfa' pour la voilure
#         seule alfa -> Cz(alfa), Cx(alfa), Cm(alfa)
#         :TODO: l'influence des freins 'd' doit être évaluée par ailleurs...
#         NB (mars 2019): pour obtenir la position de la voile en soufflerie, à vitesse V
#         donnée, il suffit d'utiliser le simulateur itérativement avec une charge
#         (poids pilote) variable, jusqu'a obtenir une vitesse de vol = V
#         La méthode self.alfaV(v0) fait le job
#         [Fixed] => NB (mars 2019) ATTENTION,
#             cmarc donne le moment des forces par rapport au point
#                 BINP9{..., 'RMPX': 0.0, 'RMPY':0.0, 'RMPZ':0.0}
#             valeur par defaut qu'on retrouve dans tous les cmi et cmo
#         """
#
#         #valeurs par defaut
#         dump = readPyData(DATA_DIR/'preferences'/'equilibre2d.dict')
#         dump.update(kargs['userdata'])
#         dump.update(kargs['pdata'])
#
#
#         #Les valeurs de b0, S0,.. et les fonctions Cx, Cz..
#         #dans Polaires
#         pol = kargs['polaires']
#         dump.update(
#             CxvCMARC = lambda alfa:pol.Cx(alfa%pi),
#             Cmv      = lambda alfa:pol.Cm(alfa%pi),
#             Cz       = lambda alfa:pol.Cz(alfa%pi),
#             b0=pol.b0, l0=pol.l0, vinf=pol.vinf, S0=pol.S0
#             )
#
#         #On écrase tout par les paramètres d'entrée
#         dump.update(kargs)
#
#         #on se débarasse du paramètre session qui est en @property
#         if 'session' in dump :
#             dump.pop('session')
#
#         #On a toutes les data dans dump, on peut instancier les attributs
#         for key, value in dump.items():
#             setattr(self, key, value)
#
#         # ne pas oublier, pour calculer les différentes constantes utiles.
#         self.load()
#
#     def Cxv(self, alfa) :
#         """fonction Cx voile"""
#         # alfa = alfa%pi
#         return self.CxvCMARC(alfa) + self.Cx0
#
#     def CxT(self, alfa) :
#         """fonction alfa ⟼ Cx(𝛼) = coef traînée total"""
#         return self.Cxv(alfa) + self.Cxs + self.Cxp
#     # CxT = Cx
#
#     def CmT(self, alfa):
#         """fonction Coef de moment global(BA), aile + pilote + suspentes"""
#         # alfa = alfa%pi
#         Cmv = self.Cmv(alfa)
#         Cz = self.Cz(alfa)
#         θ = self.θ(alfa)
#         return self.mx + θ*self.mz + (Cmv+self.Cmz-alfa*self.Cmx)/Cz
#
#     def γ(self, alfa) :
#         """fonction alfa ⟼ γ(alfa) angle de plané"""
#         # alfa = alfa%pi
#         return self.CxT(alfa)/self.Cz(alfa)
#
#     def θ(self, alfa):
#         """fonction alfa ⟼ θ(alfa) assiette"""
#         # alfa = alfa%pi
#         return self.γ(alfa)-alfa
#
#     def CP(self, alfa):
#         """Centre de poussee de la RFA en % de corde
#         = intersection de la droite d'action de la RFA avec la corde."""
#         # alfa = alfa%pi
#         return -self.Cmv(alfa)/(alfa*self.Cxv(alfa)+self.Cz(alfa))
#
#     def normalise(self,alfa):
#         """alfa ramené dans [0,pi]"""
#         return alfa%pi
#
#     def V(self, alfa) :
#         """ Vitesse air"""
#         value = (2*self.mg)/(ρ*self.Cz(alfa)*self.S0)
#         if value < 0 :
#             value = -sqrt(-value)
# #             debug("Vol dos Cz<0, alfa = %.3g°, V = %.3g km/h"%(degrees(alfa),-3.6*value))
#             return value
#             # exit()
#         else :
#             return sqrt(value)
#
#     def typeDeVol(self, alfa):
#         """le vol avec l'incidence alfa peut être vol normal, dos, arrière ou dos-arrière"""
#         alfa = alfa%pi
# #         debug(alfa_deg=degrees(alfa))
#         if 0<=alfa<=pi/2 :
#             return 0,'vol normal'
#         elif pi/2<=alfa<=pi :
#             return 1,'vol marche arrière'
#         elif -pi/2<=alfa<=0 :
#             return 2,'vol dos'
#         elif -pi<=alfa<=-pi/2 :
#             return 3,'vol dos + marche arrière'
#         else :
#             return 4,'vol ?'
#
#     def Vz(self,alfa) :
#         """taux de chute"""
#         # alfa = alfa%pi
#         return self.V(alfa)*sin(self.γ(alfa))
#
#     def qS(self, alfa) :
#         # alfa = alfa%pi
#         return (self.mg)/self.Cz(alfa)
#
#     def resultante(self, alfa):
#         """resultante en Newtons"""
# #         Cx, Cz = self.aresultante(alfa)
# #         return Cx*self.mg,Cz*self.mg
#         # alfa = alfa%pi
#         γ = self.γ(alfa)
#         qS = self.qS(alfa)
# #         qSm = self.qS(alfa)/self.mg
# #         mg = self.mg
#         return -self.mg*γ+qS*self.CxT(alfa), -self.mg+qS*self.Cz(alfa)
#         return -self.mg*sin(γ)+qS*self.CxT(alfa), -self.mg*cos(γ)+qS*self.Cz(alfa)
#
# #     def aresultante(self, alfa):
# #         """resultante a-dimensionnée"""
# #         # alfa = alfa%pi
# #         γ = self.γ(alfa)
# #         qSm = self.qS(alfa)/self.mg
# #         return -γ+qSm*self.CxT(alfa), -qSm*self.Cz(alfa)
# #         return -sin(γ)+qSm*self.CxT(alfa), -cos(γ)+qSm*self.Cz(alfa)
#
#     def finesse(self, alfa) :
#         """finesse à incidence alfa donnée"""
#         # alfa = alfa%pi
#         return self.Cz(alfa)/self.CxT(alfa)
#
#     def alfa0(self) :
#         """le alfa de portance nulle"""
#         return newton(lambda alfa:self.Cz(alfa), radians(6.0))
#
#     def objectif(self,alfa):
#         """
#         La fonction à annuler pour trouver l'équilibre
#         Peut être que return self.CmT(alfa) suffit...
#         """
#         alfa = alfa%pi
#         cx, cz = self.resultante(alfa)
#         return self.CmT(alfa)**2 + cx**2 + cz**2
#         return cx**2 +cz**2
#
# #     def objectif1(self,alfa):
# #         """
# #         La fonction à annuler pour trouver l'équilibre.
# #         On rajoute une force extérieur (acceleration du camion) aux equations
# #         """
# #         alfa = alfa%pi
# #
# #         fx, fz = self.resultante(alfa)+self.addForce(gamma)
# #         return (self.CmT(alfa)+self.addMoment(alfa))**2 + fx**2 + fz**2
# #         return fx**2 +fz**2
#
#     def equilibre(self,**kargs) :
#         """
#         :param kargs: les variables à modifier, parmi self.VARIABLES
#         :return:float, l'incidence alfa° telle que CmT(alfa)=0
#         >>> P = Parapente('diamir')
#         >>> eq = Equilibriste(P)
#         >>> alfa = eq.equilibre(c=29.0, mp=100.0, S=28.0, d=10.0)
#         # alfa sera l'incidence d'équilibre pour le parapente P,
#         # avec un calage c, un pilote de masse mp, une surface S, et 10% de freins
#         """
# #         objectif = self.CmT
#         if kargs : self.load(kargs)
#         #disp=False pour que newton ne lève pas une exception
#         alfa, res = newton(self.objectif, radians(6.0), tol=self.EPS, full_output=True, disp=False)
#         if res.converged :
#             return alfa%pi
#         else :
#             alfa = res.root
#             res = dict(root=fmt(degrees(res.root)),iterations=res.iterations,
#                        function_calls=res.function_calls,
#                        converged=res.converged, cause=res.flag)
#             arn = norm(self.resultante(alfa%pi))
#             ar  = norm(self.resultante(alfa))
#             cm = self.CmT(alfa)
#             debug('Non CV', resultante=(arn, ar), CmT=cm, res=res)
#             if abs(ar)<self.EPS :
#                 return alfa%pi
# #             debug('Non CV', resultante=(fmt(rn,6), fmt(r,6)))
#             return nan
#     alfaEq = equilibre
#
#     def toDump(self):
#         keys = sorted(list(self.VARIABLES))
#         return dict([(key, getattr(self,key)) for key in keys])
#
#     def load(self, kargs={}):
#         """
#          :param kargs: les variables à modifier, parmi self.VARIABLES
#          >>> P = Parapente('diamir')
#          >>> eq = Equilibriste(P)
#          >>> eq.load(c=29.0, mp=100.0, S=28.0, d=10.0)
#          # sont modifiés calage c, masse pilote mp, surface S, et freins 10%
#          """
#         assert set(kargs.keys()).issubset(self.VARIABLES), \
#             """
#     Une des variables parmi %s est inconnue.
#         Les variables autorisées sont : %s,
#             """%(sorted(kargs.keys()), sorted(list(self.VARIABLES)))
#         for key, value in kargs.items() :
#             setattr(self, key, value)
#
#         #constantes
#         self.xc = self.c*self.l
#         self.zp = -self.s
#         self.m    = self.mp + self.mv
#         self.Cx0T = self.Cx0 + self.Cxp + self.Cxs
#         self.Cmx = (self.xc/self.l0)*(self.Cxp + self.Cxs)
#         self.Cmz = (self.zs/self.l0)*self.Cxs - (self.s/self.l0)*self.Cxp
#         self.mx  = (self.mp/self.m)*(self.xc/self.l0)\
#                     + (self.mv/self.m)*(self.xGv/self.l0)
#         self.mz  = (self.mp/self.m)*(self.s/self.l0)\
#                     - (self.mv/self.m)*(self.zGv/self.l0)
#         self.mg  = self.m*g
#     def alfaVOld(self, v0):
#         """Retourne le alfa tel que V(alfa) = v0. Ne correspond pas à l'équilibre."""
#         try :
#             return newton(lambda alfa : self.V(alfa)-v0, radians(6.0))
#         except RuntimeError as msg : return nan
#
#     def alfaV(self, v0):
#         """Retourne le alfa tel que V(alfa) = v0. Ne correspond pas à l'équilibre."""
# #         return newton(lambda alfa : self.V(alfa)-v0, radians(6.0), full_output=True, disp=False)
# #         except RuntimeError as msg : return nan
#         alfa, res = newton(lambda alfa : self.V(alfa)-v0, radians(6.0), tol=self.EPS, full_output=True, disp=False)
#         if res.converged :
#             return alfa%pi
#         else :
#             alfa = res.root
#             res = dict(root=fmt(degrees(res.root)),iterations=res.iterations,
#                        function_calls=res.function_calls,
#                        converged=res.converged, cause=res.flag)
#             rn = norm(self.resultante(self.normalise(alfa)))
#             r  = norm(self.resultante(alfa))
#             cm = fmt(self.CmT(alfa))
#             debug('alfaV : non CV')
#             print("alfa_normalisé = %s,"%fmt(degrees(alfa%pi)), 'resultante = (%s,%s),'%fmt((rn, r)), 'CmT =', cm)
#             print('résultats newton = ', fmt(res))
#             print()
#             if abs(r)<self.EPS :
#                 #NOOOON, si la resultante est nulle et le moment non nul,
#                 #alors c'est un couple, on n'est pas à l'équilibre
#                 #return alfa%pi
#                 return nan
#
#     def cV(self, v0):
#         """Retourne le calage c pour atteindre la vitesse v0"""
#         return newton(lambda alfa : self.V(alfa)-v0, radians(6.0))
#
#     def VEq(self, **kargs):
#         """
#         :param kargs : dict(VARIABLES)
#         :returns la vitesse à l'équilibre"""
#         return self.V(self.alfaEq(**kargs))
#
#     def 𝞥max(self):
#         """:returns : la finesse max accessible, et le calage correspondant.
#             masse pilote fixéé"""
#         raise NotImplementedError("TODO")
#
#     def write(self, filename, **kargs):
#         """
#         :param freins:tuple, les % de freins pour lesquels on calcule l'équilibre, ∈[0.0, 100.0]
#         :type filename: fichier ou str
#         """
#
#         if isinstance(filename, str) :
#             fout = open(filename,'w')
#         else : #fichier ouvert. p. ex. sys.stdout
#             fout = filename
#         dump = fmt(self.toDump())
#         print(sapero(), maintenant(), "\n\n", file=fout)
#         with fout :
#             titre=sursouligne("  Lexique et valeurs des paramètres  ","=")
#             print(titre,file=fout)
#             for key in sorted(self.VARIABLES) :
#                 print("%5s = %-6s = %s"%(key, dump[key], self.VARIABLES[key]), file=fout)
#             print("\nN.B. la valeur 'nan' (Not A Number), pour signifier que le calcul n\'a pas abouti (plantage, non convergence..)", file=fout)
#             print("""\nN.B. suivant la valeur de l'incidence, le vol peut être d'un des types suivants (dans les 2 ou 3 derniers cas... à utiliser avec circonspection !):
#             .   > vol normal,
#             *   > vol marche arrière,
#             **  > vol dos,
#             *** > vol dos + marche arrière.""", file=fout)
#
#             self._writeEquilibre(fout)
#             for key, values in kargs.items() :
#                 if key=='mp' :
#                     self._writeVariationMassePilote(values, fout)
#                 if key=='c' :
# #                     rdebug(values)
#                     self._writeVariationCalage(values, fout)
#
#     def _writeEquilibre(self, fout):
#         eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu
#         titre = sursouligne('  Equilibre  ',"=")
#         print('\n\n',file=fout)
#         print(titre,file=fout)
#         if 0 : print("Paramètres:\n%s"%fmt(self.toDump()),file=fout)
#         CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
#         DELTAs, THETAs, ALPHAs = [], [], []
#         alfa = self.equilibre()
#         ALPHAs.append(alfa*180/pi)
#         FINs.append(self.finesse(alfa))
#         Vs.append(self.V(alfa)*3.6)
#         VZs.append(self.Vz(alfa))
#         CPs.append(100*self.CP(alfa))
#         THETAs.append(self.θ(alfa)*180/pi)
#         r = self.resultante(alfa)
#         REs.append(abs(r[0]) + abs(r[1]) + abs(self.CmT(alfa)))
#         CMTs.append(self.CmT(alfa))
#         CZs.append(self.Cz(alfa))
#         CXVs.append(self.Cxv(alfa))
#         CMVs.append(self.Cmv(alfa))
#         # debug(titre="Je ne sais plus d'ou vient la formule pour simuler les freins !!")
#         print(" Frein (%)           |", " ".join(['%-7.1f'%c for c in DELTAs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
#         print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
#         print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
#         print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
#         print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
#         print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
#         print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
#         print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
#         print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
#         self.load(eldump0)
#
#     def _writeVariationMassePilote(self, Ms, fout):
#         eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu
#         eldump0['mp'] = '***'
#         print("\n\n",file=fout)
#         titre = sursouligne("  Variation masse pilote  ","=")
#         print(titre,file=fout)
#         if 0 : print("paramètres :\n", fmt(eldump0), file=fout)
#         eldump0 = self.toDump()#On le sauve car je crois qu'on le modifie un peu
#
#         CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
#         THETAs, ALPHAs, vols = [], [], []
#         for c in Ms :
#             try :
#                 alfa = self.alfaEq(mp=c)
#             except RuntimeError as msg :
#                 debug(msg, calage=c)
#                 alfa = nan
#             idx, _ = self.typeDeVol(alfa)
#             vols.append('.' if idx==0 else idx*"*")
#             ALPHAs.append(alfa*180/pi)
#             FINs.append(self.finesse(alfa))
#             Vs.append(self.V(alfa)*3.6)
#             VZs.append(self.Vz(alfa))
#             CPs.append(100*self.CP(alfa))
#             THETAs.append(self.θ(alfa)*180/pi)
#             r = self.resultante(alfa)
#             REs.append(abs(r[0]) + abs(r[1]) + abs(self.CmT(alfa)))
#             CMTs.append(self.CmT(alfa))
#             CZs.append(self.Cz(alfa))
#             CXVs.append(self.Cxv(alfa))
#             CMVs.append(self.Cmv(alfa))
#         rdebug(vols)
#         print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
#         print(" masse pilote (kg)   |", " ".join(['%-7.2f'%c for c in Ms]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
# #         if not any(vols) :
#         print(" Type de vol         |", " ".join(['%-7.2s'%c for c in vols]),file=fout)
#         print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
#         print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
#         print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
#         print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
#         print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
#         print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
#         print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
#         print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
#         print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in Ms]),file=fout)
#         self.load(eldump0)
#
#     def _writeVariationCalage(self, Cs, fout):
#         eldump0=self.toDump()
#         eldump0['c'] = '###'
#         titre=sursouligne("  Variation calage  ","=")
#         print("\n\n",file=fout)
#         print(titre,file=fout)
#         if 0 : print("paramètres:\n", fmt(eldump0), file=fout)
#         eldump0=self.toDump()
#         CZs, CXVs, CMVs, REs, Vs, VZs, CPs, FINs, CMTs = [], [], [], [], [], [], [], [], []
#         THETAs, ALPHAs, typevols = [], [], []
#         for c in Cs :
#             try :
#                 alfa = self.alfaEq(c=c)
#             except RuntimeError as msg :
#                 debug(msg, calage=c)
#                 alfa = nan
#                 # continue
#             if alfa < 0 :
#                 alfa = nan
#                 alfa = alfa%pi
#                 debug ('Fermeture, vol impossible')
#             idx, _ = self.typeDeVol(alfa)
#             typevols.append('  .' if idx==0 else ' '+idx*"*")
#             ALPHAs.append(alfa*180/pi)
#             FINs.append(self.finesse(alfa))
#             Vs.append(self.V(alfa)*3.6)
#             VZs.append(self.Vz(alfa))
#             CPs.append(100*self.CP(alfa))
#             THETAs.append(self.θ(alfa)*180/pi)
#             r = self.resultante(alfa)
#             REs.append(abs(r[0]) + abs(r[1]) + abs(self.CmT(alfa)))
#             CMTs.append(self.CmT(alfa))
#             CZs.append(self.Cz(alfa))
#             CXVs.append(self.Cxv(alfa))
#             CMVs.append(self.Cmv(alfa))
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Calage (% de corde) |", " ".join(['%-7.1f'%(100*c) for c in Cs]),file=fout)
#         print(" Type de vol         |" + " ".join(['%-7s'%c for c in typevols]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Incidence (°)       |", " ".join(['%-7.2f'%c for c in ALPHAs]),file=fout)
#         print(" Vitesse (km/h)      |", " ".join(['%-7.2f'%c for c in Vs]),file=fout)
#         print(" Vz (m/s)            |", " ".join(['%-7.2f'%c for c in VZs]),file=fout)
#         print(" Finesse             |", " ".join(['%-7.2f'%c for c in FINs]),file=fout)
#         print(" Assiette (°)        |", " ".join(['%-7.2f'%c for c in THETAs]),file=fout)
#         print(" Centre poussee (%)  |", " ".join(['%-7.2f'%c for c in CPs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Cz                  |", " ".join(['%-7.2f'%c for c in CZs]),file=fout)
#         print(" Cx voile            |", " ".join(['%-7.2f'%c for c in CXVs]),file=fout)
#         print(" Cm voile            |", " ".join(['%-7.2f'%c for c in CMVs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#         print(" Resultante(=0 ?)    |", " ".join(['%-7.0g'%c for c in REs]),file=fout)
#         print(" Moment(=0 ?)        |", " ".join(['%-7.0g'%c for c in CMTs]),file=fout)
#         print("---------------------|" + "-".join([7*'-' for c in CPs]),file=fout)
#
#         #On rétablit dans sa config originale
#         self.load(eldump0)
#
#     def plot(self, **kargs):
#         """TODO : tracés pour variation calage, longueur cône suspentage, CxV, Cxs, Cxp,..."""
# #         csfont = {'fontname':'Comic Sans MS'}
#         hfont = {'fontname':'Liberation Serif'}
#
#         if not kargs :
#             kargs = readPyData(Path(DATA_DIR)/'preferences'/'aero2d.dict')
#
#         if 'mp' in kargs :
#             Vs, VZs, CPs, FINs, ALPHAs = [], [], [], [], []
#             Mp = kargs['mp']
#             for m in Mp :
#                 alfa = self.alfaEq(mp=m)
#                 if alfa < 0 :
# #                     alfa = self.normalise(alfa)
#                     debug ('Fermeture, vol impossible; alfa=%s'%fmt(degrees(alfa)))
#                     alfa = nan
#                     alfa = alfa%pi
# #                 idx, _ = self.typeDeVol(alfa)
# #                 typevols.append('  .' if idx==0 else ' '+idx*"*")
#                 ALPHAs.append(alfa*180/pi)
#                 FINs.append(self.finesse(alfa))
#                 Vs.append(self.V(alfa)*3.6)
#                 VZs.append(self.Vz(alfa))
#                 CPs.append(100*self.CP(alfa))
#                 r = self.resultante(alfa)
# #                 REs.append(abs(r[0]) + abs(r[1]) + abs(self.CmT(alfa)))
# #                 CMTs.append(self.CmT(alfa))
# #                 CZs.append(self.Cz(alfa))
# #                 CXVs.append(self.Cxv(alfa))
# #                 CMVs.append(self.Cmv(alfa))
#
# #         adeg = linspace(min(self.alfadegs), max(self.alfadegs), 100)
# #         arad = linspace(min(self.alfas), max(self.alfas), 100)
# #         Cxs, Czs = self.Cx(arad), self.Cz(arad)
#         debug(finesses=fmt(FINs))
#         debug(Vz=fmt(VZs))
#         fig = plt.figure()
#         ax1 = fig.add_subplot(221)
#         ax2 = fig.add_subplot(222, sharex=ax1)
#         ax3 = fig.add_subplot(223, sharex=ax1)
#         ax4 = fig.add_subplot(224, sharex=ax1)
# #         ax4.set_xlim(left=0.0, right=max(Cxs))
# #         ax4.set_ylim(bottom=0.0, top=max(Czs))
#         ax1.plot(Mp, FINs, 'b.-')#, label='$ \\phi $')
# #         ax1.plot(self.alfadegs, self.Czs, 'b.', label='data')
# #         ax1.legend()
#         ax1.set_xlabel('$ m_p (kg) $')
#         ax1.set_ylabel('$ \\phi $')
#         ax1.set_title('Finesses',**hfont)
#
#         ax2.plot(Mp, Vs, 'b.-')#, label='$V$')
#         ax2.set_xlabel('$ m_p (kg) $ ')
#         ax2.set_ylabel('$ V (km/h) $')
#         ax2.set_title('Vitesse sur trajectoire (km/h)',**hfont)
# #         ax2.legend()
#         ax3.plot(Mp, VZs, 'b.-')#, label='data')
#         ax3.set_ylabel('$ V_z (m/s) $')
#         ax3.set_title('Taux de chute (m/s)',**hfont)
#         if 0 :
#             ax4.plot(Mp, ALPHAs, 'b.-')#, label='data')
#             ax4.set_ylabel('$ \\alpha (°) $')
#             ax4.set_title('Incidences',**hfont)
#         if 1 :
#             ax4.plot(Mp, CPs, 'b.-')#, label='data')
#             ax4.set_ylabel('$ CP $ ')
#             ax4.set_title('Centre de poussée % de corde',**hfont)
#         fig.suptitle('Session %s : variation de la masse pilote $m_p$\n'%AC.SESSION, **hfont)
#         plt.tight_layout()
#         plt.show()
    
    def tPilote(self,alfa):
        tp = Torseur(A=Point(0,0,0), R=Vecteur())


if __name__ == '__main__':
    import os
#     import aperoconfig as AC
    py3 = '/Users/puiseux/anaconda2/envs/py3/bin/python'
    exe = AC.SOURCES_DIR/'aide'/'simul_FFVL_tests_structure.py'
    debug("Version à utiliser avec %s"%exe)
    os.system('%s %s'%(py3, str(exe)))
#     from tests.aerodynamique.testsaero import Tests
#     session = 'Whizz2'
#     session = 'Whizz2#elasticite#'
#     session = 'elastique'
# #     t = Tests(session)
#     # show = True
#     if 1 : testEquilibre(show=True)

