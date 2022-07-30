import numpy as np
from eratosthenes import *

def intersect2HalfOpen(I,J):
    a,b = I
    c,d = J
    if d <= a or b <= c:
        return -1
    else:
        return (max(a,c),min(b,d))

def intersect1HalfOpen(I,J):
    a,b = I
    c,d = J
    if d <= a or b < c:
        return -1
    else:
        return (max(a,c),min(b,d))

def sieveWithKnownPrimes(P,a,b):
    if b < a:
        return []
    a,b = int(a),int(b)
    length = b-a
    B = np.ones(b-a+1,dtype=bool)

    for p in P:
        start = int(p * (np.ceil(a/p))) - a
        k = start
        while k <= b-a:
            B[k] = 0
            k += p

    return np.arange(a,b+1)[B]

def calcPi(start,B,x):
    return start + len(B[np.argwhere(B <= x)])

def P2(n,P=None):
    n13 = n**(1/3)
    N = int(n13 + epsilon)
    n12 = n**(1/2)

    if P is None:
        P = eratosthenesSieve(N)
    P14 = P[np.argwhere(P <= n**(1/4))].flatten()
    pi13 = len(P)
    pi12 = 0
    pi = pi13
    S = 0
    for j in range(2,int(n**(2/3)/N)+2):
        a,b = (j-1)*N+1, min(j*N,int(n**(2/3)))

        B = sieveWithKnownPrimes(P,a,b)

        if a <= int(n12) <= b:
            pi12 = calcPi(pi,B,n12)

        I = intersect2HalfOpen((n/(b+1),n/a),(n13,n12))
        if I != -1:

            PI = sieveWithKnownPrimes(P14,int(I[0]+1),int(I[1]))
            if len(PI) > 0:
                PI = (n/PI).astype("int")
                L = np.vectorize(lambda x: calcPi(pi,B,x))
                S += np.sum(L(PI))
        pi += len(B)

    valP2 = pi13 * (pi13 -1)/2 - pi12 * (pi12 - 1)/2+ S
    return int(valP2)

def specialSieve(n):
    n = int(np.floor(n))
    N = np.ones(n+1,dtype=bool)
    f = np.zeros(n+1,dtype=int)
    mu = np.ones(n+1,dtype=int)

    #Identify the elements of this array with the numbers where the ith cell matches with the ith number
    N[1] = N[0] = 0
    sqn = int(np.ceil(np.sqrt(n)))
    for p in range(2,sqn+1):
        if N[p] == 0:
            continue
        else:
            mu[p] = -1
            f[p] = p
            k = 2*p
            while k <= n:
                N[k] = 0
                mu[k] *= -1
                if f[k] == 0:
                    f[k] = p
                k += p

    P = np.arange(n+1)[N]
    #primes bigger sqrt n
    Pover12 = P[P > sqn]
    for p in Pover12:
        mu[p] = -1
        f[p] = p
        k = p * 2
        while k <= n:
            mu[k] *= -1
            if f[k] == 0:
                f[k] = p
            k += p

    P12 = P[np.argwhere(P <= n ** (1 / 2))].flatten()
    for p in P12:
        k = p**2
        while k <= n:
            mu[k] = 0
            k += p**2

    return P, f[1:],mu[1:]

epsilon = 10**-10

def getBinaryExponents(x):
    exponents = []
    e = int(np.log2(x))
    while x != 0:
        q,r = np.divmod(x,2**e)
        if q == 1:
            exponents.append(e)
        x = r
        e -= 1
    return np.array(exponents)

def S2(x,F=None,P=None):
    N = int(x ** (1 / 3) + epsilon)
    if F is None or P is None:
        P, fSieve, muSieve = specialSieve(N)
        F = np.array((np.arange(1, N + 1), fSieve, muSieve)).T

    phi = np.zeros(len(P))
    exponentsN = getBinaryExponents(N)
    phiIndexN = []
    for e in exponentsN:
        be = 1 + np.sum(2 ** (exponentsN[exponentsN > e] - e))
        phiIndexN += [(e,be-1)]

    valS2 = 0
    for j in range(1, int(x ** (2 / 3) / N) + 2):
        c, d = (j-1) * N + 1, min(j * N, int(x ** (2 / 3)+epsilon))
        a = 2 ** (np.arange(int(np.log2(N))+1).reshape((-1, 1)) @ np.ones(N+1).reshape((1, -1)))
        for i in range(1,int(np.log2(N))+1):
            jBound = int(N/2**i)
            a[i,jBound:] = -1
            a[i,jBound] = np.remainder(N,2**i)

        for b, p in enumerate(P):
            J = intersect1HalfOpen((x / ((d + 1) * p), x / (c * p)), (1, N))
            if J != -1:
                L, U = int(J[0]+1), int(J[1])
                FLU = F[L-1:U]
                for m, _, mu in FLU[(FLU[:, 1] > p) & (FLU[:, 2] != 0)]:
                    if m * p <= N:
                        continue
                    phiyb = phi[b]
                    y = int(x / (m * p) + epsilon)
                    l = y - (c - 1)
                    exponents = getBinaryExponents(l)
                    for e in exponents:
                        be = 1 + np.sum(2 ** (exponents[exponents > e] - e))
                        phiyb += a[e, be - 1]
                    valS2 -= mu * phiyb

            for index in phiIndexN:
                phi[b] += a[index]

            #k = l -1 to address indices properly
            start = int(p * (np.ceil(c / p))) - c
            k = start
            while k <= d - c:
                if a[0,k] == 1:
                    a[0, k] -= 1
                    for i in range(1, a.shape[0]):
                        lInd = ((k+1) + 2 ** i - 1) / (2 ** i) + epsilon
                        a[i, int(lInd)-1] -= 1
                k += p

    return int(valS2)

def S1(x,mu=None):
    N = int(x ** (1 / 3) + epsilon)
    if mu is None:
        _, _, muSieve = specialSieve(N)
    ks = np.arange(1, N + 1)
    du = (x / ks + epsilon).astype("int")
    valS1 = np.dot(mu, du)
    return int(valS1)

def phi(x):
    N = int(x**(1/3) + epsilon)
    P,fSieve,muSieve = specialSieve(N)
    F = np.array((np.arange(1, N + 1), fSieve, muSieve)).T
    valS1 = S1(x,muSieve)
    valS2 = S2(x,F,P)

    return valS1 + valS2

def lmoMethod(x):
    N = int(x ** (1 / 3) + epsilon)
    P, fSieve, muSieve = specialSieve(N)
    F = np.array((np.arange(1, N + 1), fSieve, muSieve)).T
    valS1 = S1(x, muSieve)
    valS2 = S2(x, F, P)
    valP2 = P2(x,P)
    pi13 = len(P)

    return valS1 + valS2 - valP2 + pi13 -1

