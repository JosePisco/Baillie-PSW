#include <math.h>
#include <stddef.h>

#include "crypto_utils.h"

int modinv(int x, int n)
{
    /* C does not naturally handle negative modulo */
    while (x < 0)
        x += n;
    /* No need to check for GCD being equal to 1 since it will be
    only used to inverse 2 mod n, n being odd so they are coprime
    and the inverse will always exist */

    int b0 = n;
    int t = 0;
    int q = 0;
	int x0 = 0;
    int x1 = 1;
	if (n == 1)
        return 1;

    while (x > 1) {
		q = x / n;
		t = n;
        n = x % n;
        x = t;
		t = x0;
        x0 = x1 - q * x0;
        x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

int modexp(int x, int y, int p)
{
    int result = 1;
    x %= p;

    /* if p divides x */
    if (x == 0)
        return 0;

    while (y > 0) {
        if (y % 2 == 1)
            result = (x * result) % p;

        y /=2;
        x = (x * x) % p;
    }

    return result;
}

int jacobi_symbol(uint64_t d, uint64_t n)
{
    d %= n;
    int result = 1;

    while (d != 0) {
        while (d % 2 == 0) {
            d /= 2;
            uint64_t r = n % 8;
            if (r == 3 || r == 5)
                result = -result;
        }
        uint64_t tmp = d;
        d = n;
        n = tmp;

        if (d % 4 == 3 && n % 4 == 3)
            result = -result;

        d %= n;
    }

    if (n == 1)
        return result;
    return 0;
}

struct lucas_sequence lucas(unsigned long long k, unsigned long long D, unsigned long long P, unsigned long long n)
{
    unsigned long long U = 1;
    unsigned long long V = P;
    unsigned long long k_obj = k;
    k = 1;

    size_t bitlength = log2(k_obj) + 1;
    /* starts at bit 2, bit one is done by setting U=1 and V=P */
    size_t i = 2;
    /* Iterates over the bits of n from left to right */
    while (i != bitlength + 1)
    {
        /* get the ith bit (from left being 0 to right) of a number*/
        size_t bit = (k_obj >> (bitlength - i)) % 2;

        if (bit == 0) {
            unsigned long long tmp_U = U;
            U = U * V % n;
            /* A divison by two is equivalent to a multiplication by the mod inverse of 2 by n*/
            V = (V*V + D * tmp_U*tmp_U) * modinv(2, n) % n;
            //V %=n;
            k *= 2;
        }
        else {
            unsigned long long tmp_U = U;
            U = U * V % n;
            /* A divison by two is equivalent to a multiplication by the mod inverse of 2 by n*/
            V = (V*V + D * tmp_U*tmp_U) * modinv(2, n) % n;
            k *= 2;

            tmp_U = U;
            U = (P * U + V) * modinv(2, n) % n;
            V = (D * tmp_U + P * V) * modinv(2, n) % n;
            k += 1;
        }
        i++;
    }

    struct lucas_sequence result = {0};
    result.U = U;
    result.V = V;

    return result ;
}
