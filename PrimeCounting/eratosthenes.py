import numpy as np

def eratosthenesSieve(n):
    if n <= 1:
        return np.array([])
    n = int(np.floor(n))
    N = np.ones(n+1,dtype=bool)
    #Identify the elements of this array with the numbers where the ith cell matches with the ith number
    N[1] = N[0] = 0
    sqn = int(np.ceil(np.sqrt(n)))
    for p in range(2,sqn+1):
        if N[p] == 0:
            continue
        else:
            k = p**2
            while k <= n:
                N[k] = 0
                k += p

    return np.arange(n+1)[N]

