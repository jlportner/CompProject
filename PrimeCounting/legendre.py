import numpy as np
from eratosthenes import eratosthenesSieve

def legendresMethod(n):
    P = eratosthenesSieve(np.sqrt(n)).tolist()
    pi12 = len(P)
    pi = pi12 + int(n) - 1
    k = 1
    while k <= len(P) and np.prod(P[:k]) < n:
        Sk = 0

        def legendreRec(N,prod,length):
            nonlocal Sk
            if length == 1:
                for x in N:
                    if prod * x > n:
                        break
                    Sk += int(n / (prod*x))
            else:
                for i, x in enumerate(N):
                    if prod * x > n:
                        break
                    else:
                        legendreRec(N[(i+1):],prod * x,length-1)

        legendreRec(P,1,k)
        pi += (-1) ** k * Sk
        k += 1

    return pi