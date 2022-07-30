import numpy as np
from eratosthenes import *

def preCalcPhi(n, a,P=None):
    if P is None:
        P = eratosthenesSieve(n).tolist()

    def preCalcRec(x,k):
        if k <= 0:
            return int(x)
        else:
            return preCalcRec(x, k - 1) - preCalcRec(x / P[k-1], k - 1)

    preCalc = np.vstack((np.arange(n+1),np.zeros((a,n+1))))
    for i in range(1,a+1):
        preCalc[i] = preCalc[i-1] - preCalc[i-1][(np.arange(n+1) / P[i-1]).astype("int")]

    return preCalc

epsilon = 10 ** -10
def lehmer(x, alternative=False):
    N = int(x ** (1 / 3) + epsilon)
    P = eratosthenesSieve(N)
    valPhi = phi(x,P)
    valP2 = P2(x)
    piVal = valPhi - valP2 - 1 + len(P)
    return int(piVal)

def phi(x,P=None):
    if P is None:
        N = int(x ** (1 / 3) + epsilon)
        P = eratosthenesSieve(N)
    k = int(np.min((len(P),5)))
    modVal = int(np.prod(P[:k]))
    phiCache = preCalcPhi(modVal, k,P)
    phiVal = phiRec(x, np.max((len(P),0)), P, k, phiCache,modVal)
    return phiVal

def phiRec(x, l, P, k, phiCache, modVal):
    if l <= k:
        q, r = np.divmod(x, modVal)
        return phiCache[l,modVal] * q + phiCache[l,int(r)]
    elif x < P[l-1]:
        return 1
    else:
        return phiRec(x, l - 1, P, k, phiCache, modVal) - phiRec(x / P[l-1], l - 1, P, k, phiCache, modVal)

def P2(x):
    n = x ** (1 / 3) + epsilon
    n12 = x ** (1 / 2) + epsilon
    P = eratosthenesSieve(x/n)
    b = np.searchsorted(P,n12,side="left")
    a = np.searchsorted(P,n,side="right")
    valP2 = 0
    curPi = b
    curP = P[np.searchsorted(P,n12,side="right"):]
    for i in range(b-1,a-1,-1):
        val = x/P[i]
        valIndex = np.searchsorted(curP,val,side="right")
        curPi += valIndex
        curP = curP[valIndex:]
        valP2 += curPi

    return valP2 - (b-a)*(b+a-1)/2
