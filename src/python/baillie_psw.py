from math import sqrt
from Crypto.Util.number import GCD
from miller_rabin import miller_rabin, generate_basis

primes = generate_basis(1000)

def jacobi_symbol(d, n):
    d %= n
    result = 1
    while d != 0:
        while d % 2 == 0:
            d //= 2
            r = n % 8
            if r == 3 or r == 5:
                result = -result

        d, n = n, d
        if d % 4 == 3 and n % 4 == 3:
            result = -result
        d %= n

    if n == 1:
        return result
    return 0

def lucas(k, D, P, Q, n):
    assert D != 0
    #assert Q == (1-D) / 4

    U = 1
    V = P
    k_obj = k
    k = 1
    for i in bin(k_obj)[3:]:
        if i == '0':
            U, V = U * V, (V**2 + D*U**2) * pow(2, -1, n) % n
            k *= 2
        else:
            U, V = U * V, (V**2 + D*U**2) * pow(2, -1, n) % n
            k *= 2
            U, V = (P * U + V) * pow(2, -1, n) % n, (D* U + P * V) * pow(2, -1, n) % n
            k += 1

    assert k_obj == k
    return (U, V)

def strong_lucas_test(n, D, P, Q):
    assert GCD(n, D) == 1
    # Find d, s such as
    # n+1 = d * 2^s        ; d is odd, n+1 even
    d = 0
    s = 1
    while True:
        d = (n+1) // pow(2, s)
        if d % 2 == 1:
            break
        s += 1

    # d is odd
    # using the lucas sequences formulas:
    # U_d mod n = 0 or V_(d*2^r) mod n = O
    #print("d:",d)
    #print("s:", s)
    #print("D:", D)
    #print("P:", P)
    #print("Q:", Q)

    u_d = lucas(d, D, P, Q, n)[0]
    #print("u_d:", u_d)

    if u_d % n == 0:
        return True
    else:
        for r in range(s):
            v_d = lucas(d * pow(2, r), D, P, Q, n)[1]
            #print("v_d:", v_d)
            if v_d % n == 0:
                return True
    return False

"""
    Test N for primality using the strong Lucas test with Selfridge's
    parameters. Returns 1 if N is prime or a strong Lucas-Selfridge
    pseudoprime (in which case N is also a pseudoprime to the standard
    Lucas-Selfridge test). Returns 0 if N is definitely composite.

    The running time of the strong Lucas-Selfridge test is, on average,
    roughly 10 % greater than the running time for the standard
    Lucas-Selfridge test (3 to 7 times that of a single Miller's test).
    However, the frequency of strong Lucas pseudoprimes appears to be
    only (roughly) 30 % that of (standard) Lucas pseudoprimes, and only
    slightly greater than the frequency of base-2 strong pseudoprimes,
    indicating that the strong Lucas-Selfridge test is more computationally
    effective than the standard version.
"""
def strong_lucas_selfridge(n):
    # we also need to filter out all perfect square values of N
    # this is because we will later require an integer D for which Jacobi(D,N) = -1
    # no such integer exists if N is a perfect square.
    root = int(sqrt(n)) # attention les approx de signes en float -> convertir en int
    if root*root == n:
        return False

    # Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
    # such that Jacobi(D,N) = -1 (Selfridge's algorithm). Theory
    # indicates that, if N is not a perfect square, D will "nearly
    # always" be "small."
    d_abs = 5
    sign = 1
    D = d_abs * sign
    while True:
        if GCD(n, d_abs) > 1: # if 1 < GCD < n, then n is composite and d is a factor of n
            return False

        if jacobi_symbol(D, n) == -1:
            break

        d_abs += 2
        sign = -sign
        D = d_abs * sign

    # selfridge reference
    P = 1
    Q = (1 - D) // 4 # is always an integer
    if D < 0:
        #print("ICIIIIIIIII")
        D = pow(D, -1, n)
    return strong_lucas_test(n, D, P, Q)

def baillie_psw(n):
    if n <= 1:
        return False

    if n % 2 == 0:
        return False

    # there are cases where this may save some time
    # test with first primes < 1000 if they divide n
    for p in primes:
        if n % p == 0 and n != p:
            return False


    # Test miller-rabin with base 2 (any base might work in theory but the work to prove its efficiency combined with lucas test has been conducted on base 2)
    if miller_rabin(n, 2) == False:
        return False

    # perform a strong Lucas-Selfridge test on N using Lucas sequences with the paramaters suggested by Selfridge
    return strong_lucas_selfridge(n)

car = 5494192362793378419461850103325205105746618253649659755062958621743505179898846930508801354093259555029072978427447589979886135102666157147211855994196065010559627866228837393549736942388413848102819329
arn = 18554525431820824592280314043934296355402115917852377301311719318802586275716645289777262350595688481395192557164440616580803261204902024969030829323747390730462863817347016368983019917352880671694306089467

#from Crypto.Util.number import getPrime
#prime = getPrime(1024)

assert baillie_psw(19) == True
#assert baillie_psw(car) == False
#assert baillie_psw(arn) == False
#assert baillie_psw(prime) == True

for i in primes:
    if baillie_psw(i) == False:
        print(i)
