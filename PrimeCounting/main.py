from eratosthenes import eratosthenesSieve
from legendre import legendresMethod
from lmo import lmoMethod
from lehmer import lehmer
from timeit import default_timer as timer

while True:
    print("Enter x to calculate pi(10^x): ")
    a = int(input())
    n = 10**(a)

    print("Calculating number of primes up to " + str(n))

    if a >= 9:
        print("Running Erastothenes Sieve will take longer than 1min, skipping it...")
    else:
        print("Running Erastothenes Sieve ...")
        s = timer()
        print(len(eratosthenesSieve(n)))
        e = timer()
        print("Took: " + str(e-s) + "s")

    if a >= 8:
        print("Running Legendre will take longer than 1min, skipping it...")
    else:
        print("Running Legendres Method ...")
        s = timer()
        print(legendresMethod(n))
        e = timer()
        print("Took: " + str(e-s) + "s")


    print("Running Lehmer-Method ...")
    s = timer()
    print(lehmer(n))
    e = timer()
    print("Took: " + str(e - s) + "s")

    print("Running LMO-Method ...")
    s = timer()
    print(lmoMethod(n))
    e = timer()
    print("Took: " + str(e - s) + "s")
