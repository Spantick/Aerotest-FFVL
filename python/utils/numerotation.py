from utils.debog import debug
# from matplotlib import pyplot as plt
import numpy as np
from random import shuffle, sample
from numpy import vstack, asarray, where, array, empty_like, ndarray, hstack
"""
Une permutation de [0,n] et une application 𝝈:[0,n]→[0,n] qui
représente un changement de numérotation où "𝝈[3] est le NOUVEAU numéro de 3".
La permutation 𝝈 est représentée par sa liste que l'on nomme aussi 𝝈.
Par exemple 𝝈 = [2,0,1,4,3] est la permutation de [0,1,2,3,4] où :
    0→2,
    1→0,
    2→1,
    3→4
    4→3
ce que l'on symbolise par [0,1,2,3,4]→[2,0,1,4,3]
𝝈 est la liste des NOUVELLES positions.
Attention : lorsqu'on dit "les 5 points sont ordonnés ainsi : o = [1, 2, 0, 4, 3]"
on exprime que l'on met en tête le point n°1, puis le 2, le 0, le 4 et le 3.
[1,2,0,4,3] est la liste des ANCIENNES positions
Dans ce cas, o est l'inverse de la permutation 𝝈. En effet il suffit de changer
le sens des flèches ci-dessus :
    0←2,             2→0,
    1←0,             0→1,
    2←1, ou encore   1→2,
    3←4              4→3
    4←3              3→4
soit (en réordonnant la 1ere colone) [0,1,2,3,4]→[1,2,0,4,3]
l'inverse de l'inverse de 𝝈 est bien sûr 𝝈 elle-même
A RETENIR :
    - le tableau des NOUVELLES positions est la permutation 𝝈
    - le tableau des ANCIENNES positions, est 'o' l'inverse de 𝝈,
    - quand on a une liste de positions L à renuméroter
        si on connait 𝝈, on appelle renuméroter(𝝈, L, invert=False)
        si on connait o, on appelle renuméroter(o, L, invert=True)
"""
def inverse(𝝈):
    """retourne l'inverse de la permutation 𝝈"""
#     𝝈1 = len(𝝈)*[None]
    𝝈1 = empty_like(𝝈)
    for i, 𝝈i in enumerate(𝝈) :
        𝝈1[𝝈i] = i
    return 𝝈1

def inverser(𝝈):
    """inversion in-situ"""
    𝝈[:] = inverse(𝝈)

def permuted(𝝈, L, invert=False):
    """
    renumérote L, ne modifie pas L
    :param 𝝈: ndarray(n,int) est la permutation, i.e. la liste des nouveaux numéros :
        𝝈[i] = nouveau numéro du point i
    :param L: ndarray(nl,int), est une liste d'indices, qui vérifie max(L)∈[0,len(𝝈)[
    :returns L1: ndarray(nl,int), L renuméroté
    """
#     L1 = L.copy()
    L1 = empty_like(L)
    if invert :
        for i, 𝝈i in enumerate(𝝈) :
            L1[where(L==𝝈i)] = i
    else :
        for i, 𝝈i in enumerate(𝝈) :
            L1[where(L==i)] = 𝝈i
    return L1

class Permutation(object):
    """Une classe qui modélise une permutation, pour renumérotation dans un maillage"""
    def __init__(self, 𝝈, invert=False):
        self.𝝈 = asarray(𝝈)
        self._invert = invert
        if len(self.𝝈.shape) != 1 :
            raise ValueError("𝝈 doit être de shape (n,)")

    def __len__(self):
        return len(self.𝝈)

    def __call__(self, N):
        return permuted(self.𝝈, N, invert=self._invert)
#         return N

    def __getitem__(self,i):
        if self._invert :
#             return self.𝝈.index(i)
            return where(self.𝝈==i)[0][0]
        else :
            return self.𝝈[i]

    def invert(self):
        self._invert = not self._invert
#         inverse(self.𝝈)
    def __str__(self):
        if self._invert :
            return str(inverse(self.𝝈).tolist())
        else :
            return str(self.𝝈.tolist())

def renumeroted(C, old, new):
    """
    Renumérote les indices de C old -> new
    Ne modifie pas C
    C1 = C dans lequel tous les indices old[i] sont remplacés par new[i]
    :param C : ndarray(nc, int), une liste d'indices à modifier
    :param old : ndarray(n,int), liste des anciens numéros
    :param new : ndarray(n,int), liste des nouveaux numéros, new.shape = old.shape
    :returns C1: C modifié.
    le plus souvent, C contient des indices de old, et peut en contenir d'autres
    """
#     debug()
#     s = 2212
#     old = asarray(old)
#     new = asarray(new)
#     print('C   =',(C-s).tolist())
#     print('(old, new)=',sorted(zip(old-s,new-s)))
    C1 = C.copy()#Obligé de faire une copy, car C est modifié dans la boucle
    for o, n in zip(old,new) :
        C1[where(C==o)] = n
#     print('C1  =',(C1-s).tolist())
    return C1

def testRenumerotation():
    debug(titre="testRenumerotation")
    C = np.asarray([[37, 38], [38, 39], [39, 40], [39, 41], [38, 42], [42, 43], [42, 44], [37, 45], [45, 46], [46, 47], [46, 48], [45, 49], [49, 50], [49, 51], [45, 52], [52, 53], [52, 54], [52, 55], [37, 56], [56, 57], [57, 58], [57, 59], [57, 60], [57, 61], [56, 62], [62, 63], [62, 64], [62, 65], [37, 66], [66, 67], [67, 68], [68, 69], [68, 70], [67, 71], [71, 72], [71, 73], [0, 1], [1, 2], [2, 3], [2, 4], [1, 5], [5, 6], [5, 7], [0, 8], [8, 9], [9, 10], [9, 11], [8, 12], [12, 13], [12, 14], [8, 15], [15, 16], [15, 17], [15, 18], [0, 19], [19, 20], [20, 21], [20, 22], [20, 23], [20, 24], [19, 25], [25, 26], [25, 27], [25, 28], [0, 29], [29, 30], [30, 31], [31, 32], [31, 33], [30, 34], [34, 35], [34, 36]])
    A = [[2, 5], [3, 4], [5, 3], [6, 2], [9, 5], [10, 4], [12, 3], [13, 2], [15, 1], [16, 1], [17, 1], [20, 5], [21, 4], [22, 5], [23, 4], [25, 3], [26, 2], [27, 3], [31, 4], [32, 3], [34, 2], [35, 1], [38, 6], [39, 7], [41, 8], [42, 9], [45, 6], [46, 7], [48, 8], [49, 9], [51, 10], [52, 10], [53, 10], [56, 6], [57, 7], [58, 6], [59, 7], [61, 8], [62, 9], [63, 8], [67, 7], [68, 8], [70, 9], [71, 10]]
    # N = [[-1.041, 0.23, -7.3], [-1.066, 0.372, -6.821], [-1.392, 1.066, -2.489], [-1.606, 0.734, -0.018], [-1.532, 2.148, -0.276], [-1.151, 2.693, -3.084], [-1.351, 3.443, -0.709], [-1.045, 4.495, -1.347], [-1.024, 0.373, -6.821], [-0.965, 1.069, -2.477], [-0.955, 0.733, -0.047], [-0.911, 2.143, -0.303], [-0.804, 2.691, -3.088], [-0.803, 3.434, -0.731], [-0.565, 4.484, -1.361], [-0.373, 4.26, -2.985], [-0.536, 5.105, -2.149], [-0.261, 5.098, -2.154], [0.077, 5.125, -2.134], [-0.978, 0.372, -6.825], [-0.429, 1.062, -2.514], [-0.318, 0.734, -0.018], [-0.304, 2.148, -0.276], [0.347, 0.739, 0.073], [0.276, 2.164, -0.2], [-0.464, 2.672, -3.11], [-0.268, 3.443, -0.709], [0.132, 4.522, -1.311], [0.219, 3.467, -0.648], [-0.938, 0.34, -6.93], [-0.381, 1.662, -4.882], [0.151, 2.146, -3.016], [0.966, 2.18, -0.124], [0.852, 3.494, -0.581], [0.002, 3.049, -3.494], [0.659, 4.559, -1.262], [0.338, 5.149, -2.116], [-1.041, -0.23, -7.3], [-1.066, -0.372, -6.821], [-1.392, -1.066, -2.489], [-1.606, -0.734, -0.018], [-1.532, -2.148, -0.276], [-1.151, -2.693, -3.084], [-1.351, -3.443, -0.709], [-1.045, -4.495, -1.347], [-1.024, -0.373, -6.821], [-0.965, -1.069, -2.477], [-0.955, -0.733, -0.047], [-0.911, -2.143, -0.303], [-0.804, -2.691, -3.088], [-0.803, -3.434, -0.731], [-0.565, -4.484, -1.361], [-0.373, -4.26, -2.985], [-0.536, -5.105, -2.149], [-0.261, -5.098, -2.154], [0.077, -5.125, -2.134], [-0.978, -0.372, -6.825], [-0.429, -1.062, -2.514], [-0.318, -0.734, -0.018], [-0.304, -2.148, -0.276], [0.347, -0.739, 0.073], [0.276, -2.164, -0.2], [-0.464, -2.672, -3.11], [-0.268, -3.443, -0.709], [0.132, -4.522, -1.311], [0.219, -3.467, -0.648], [-0.938, -0.34, -6.93], [-0.381, -1.662, -4.882], [0.151, -2.146, -3.016], [0.966, -2.18, -0.124], [0.852, -3.494, -0.581], [0.002, -3.049, -3.494], [0.659, -4.559, -1.262], [0.338, -5.149, -2.116]]
    agarder = [0, 1, 2, 5, 8, 9, 12, 15, 19, 20, 25, 29, 30, 31, 34, 37, 38, 39, 42, 45, 46, 49, 52, 56, 57, 62, 66, 67, 68, 71]
    ajeter = [3, 4, 6, 7, 10, 11, 13, 14, 16, 17, 18, 21, 22, 23, 24, 26, 27, 28, 32, 33, 35, 36, 40, 41, 43, 44, 47, 48, 50, 51, 53, 54, 55, 58, 59, 60, 61, 63, 64, 65, 69, 70, 72, 73]
    snumancs = [8, 9, 10, 21, 3, 7, 16, 20, 2, 6, 15, 17, 19, 1, 5, 12, 14, 18, 0, 4, 11, 13, 22, 26, 33, 35, 23, 27, 34, 36, 40, 24, 28, 37, 39, 41, 25, 29, 38, 42, 30, 31, 32, 43]
    vnumancs = [44, 55, 69, 0, 281, 291, 305, 237, 518, 527, 535, 543, 474, 755, 764, 772, 780, 711, 992, 1001, 1009, 1018, 1229, 1238, 1246, 1255, 1466, 1475, 1483, 1491, 1422, 1703, 1712, 1720, 1728, 1659, 1940, 1950, 1964, 1896, 2177, 2188, 2202, 2133]
    #Exemple 0
    n = 10
    C00 = list(range(n))
    shuffle(C00)
    print("C00       =", C00)
    print("inv C00   =", inverse(C00).tolist())
    print("C00  ?    =", inverse(inverse(C00)).tolist())
    C0 = asarray([0,7,3,9,4,8,5,8,1,4,5,9,2,3,0,2,7,6,6,1])
    ajeter = (1,2,5,7)
    agarder = sorted(set(C00) - set(ajeter))
    # agarder = (0,3,4,6,8,9)
    𝝈 = hstack((agarder,ajeter))
    print('agarder   =',asarray(agarder).tolist())
    print('ajeter    =',asarray(ajeter).tolist())
    𝝈a = len(agarder)
    print("C00       =", C00)
    print("𝝈         =", 𝝈.tolist(), "𝝈a =", 𝝈a)

    print("C0 avant  =", C0.tolist())
    C0 = permuted(𝝈, C0, True)
    print("C0 après  =",C0.tolist())
    shuffle(C00)
    C01 = (asarray(C00)+100).tolist()
    shuffle(C01)
    old, new = C00, C01
#     shuffle(old)
    
    print('old       =',old)
    print('new       =',new)
    C = C0.tolist() + [1000,2000]
    shuffle(C)
    C = asarray(C)
    print('C         =',C.tolist())
    print("renum(C)  =", renumeroted(C, old, new).tolist())

def testPermutation():
    debug(titre="test de la classe Permutation")
    n = 10
    C00 = list(range(n))
    shuffle(C00)
#     print("C00         = ", C00)
    𝝈 = Permutation(C00)
    print("initial :        𝝈 =", 𝝈)
    𝝈.invert()
    print("après invert() : 𝝈 =", 𝝈)
    print('__getitem__        :',[𝝈[i] for i in range(len(𝝈))])
    𝝈.invert()
    print("après invert() : 𝝈 =", 𝝈)
    print('__getitem__        :',[𝝈[i] for i in range(len(𝝈))])
    debug(paragraphe="renumerotation")
    C0 = array([0,7,3,9,4,8,5,8,1,4,5,9,2,3,0,2,7,6,6,1])
    ajeter = (1,2,5,7)
    agarder = sorted(set(C00) - set(ajeter))
    # agarder = (0,3,4,6,8,9)
    o = hstack((agarder,ajeter))#anciennes positions
    print('o =',o)
    print('___',array(range(len(o))))
    print('𝝈 =',inverse(o))
    𝝈 = Permutation(inverse(o))
    print('   C0 =', C0)
    print('𝝈(C0) =', 𝝈(C0))

if __name__ == '__main__' :
    if 1 : testPermutation()
    if 1 : testRenumerotation()
