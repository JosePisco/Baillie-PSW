#include <stdbool.h>
#include <stdio.h>
#include <math.h>

#include "baillie_psw.h"
#include "crypto_utils.h"
#include "miller_rabin.h"

/* First primes to 1000, can save some time to check if the number is divisible by one of them
might be wise to extend it to the primes until 5000, as it is not time consuming */
static int first_primes_1000[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};

bool strong_lucas_test(int n, int D, int P, int Q)
{
    (void) Q;
    /* Find d, s such as n+1 = d * 2^s
    d is odd, n+1 even */
    int d = 0;
    int s = 1;
    while (true) {
        d = (n+1) / (int) pow(2, s);
        if (d % 2 == 1)
            break;
        s++;
    }
    /*
    d is odd
    using the lucas sequences formulas:
    if U_d mod n = 0 or V_(d*2^r) mod n = O
    then n is a strong lucas pspr */

    //printf("d: %d\n", d);
    //printf("s: %d\n", s);
    //printf("D: %d\n", D);
    //printf("P: %d\n", P);
    //printf("Q: %d\n", Q);

    struct lucas_sequence lucas_res;
    lucas_res = lucas(d, D, P, n);
    unsigned long long u_d = lucas_res.U;
    //printf("u_d: %d\n", u_d);
    if (u_d % n == 0)
        return true;
    else {
        for (int r = 0; r < s; ++r) {
            lucas_res = lucas(d * (unsigned long long) pow(2, r), D, P, n);
            unsigned long long v_d = lucas_res.V;
            //printf("v_d: %d\n", v_d);
            if (v_d % n == 0)
                return true;
        }
    }

    return false;
}

bool strong_lucas_selfridge(int n)
{
    /*
    we also need to filter out all perfect square values of N
    this is because we will later require an integer D for which Jacobi(D,N) = -1
    no such integer exists if N is a perfect square.*/
    int root = sqrt(n);
    if (root * root == n)
        return false;

    /*
    Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
    such that Jacobi(D,N) = -1 (Selfridge's algorithm). Theory indicates
    that, if N is not a perfect square, D will "nearly always" be "small."
    */
    int d_abs = 5;
    int sign = 1;
    int D = d_abs * sign;
    while (true) {
        // TODO: test gcd(D, n) == 1 if not then false

        // if D is negative, take the inv mod n of D, to compute modulos efficiently
        int tmp_D = D;
        if (D < 0)
            tmp_D = modinv(D, n);
        if (jacobi_symbol(tmp_D, n) == -1)
            break;

        d_abs += 2;
        sign = -sign;
        D = d_abs * sign;
    }

    // selfridge reference
    int P = 1;
    int Q = (1 - D) / 4; // is always an integer
    if (D < 0)
        D = modinv(D, n);
    return strong_lucas_test(n, D, P, Q);
}

bool baillie_psw(int n)
{
    if (n <= 1)
        return false;

    if (n == 2)
        return true;

    if (n % 2 == 0)
        return false;

    /* there are cases where this may save some time
    test with first primes < 1000 if they divide n */
    for (int i = 0; i < NUM_FIRST_PRIMES; ++i) {
        int prime = first_primes_1000[i];
        if (n % prime == 0 && n != prime)
            return false;
    }

    if (miller_rabin(n, 2) == false)
        return false;

    return strong_lucas_selfridge(n);
}

int main(void)
{
    /* 2 does not pass but who cares ? */
    for (size_t i = 0; i < NUM_FIRST_PRIMES; ++i) {
        int prime = first_primes_1000[i];
        if (baillie_psw(prime) == false)
            printf("%d\n", prime);
    }

    //printf("%d\n", baillie_psw(2293));
    //printf("%d\n", baillie_psw(2293237));
    printf("%d\n", miller_rabin(2999, 2));

    return 0;
}
