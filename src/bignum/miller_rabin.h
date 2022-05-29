#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <openssl/bn.h>
#include <stdbool.h>

static const int BASE_2 = 2;

/* PROTOTYPES */

/*
    Miller-Rabin primality test, exclusively for base 2
    It has been adapted for base 2 test and shall not be used
    for general cases.
*/
bool BN_miller_rabin_base_2(BIGNUM *n);

#endif /* MILLER_RABIN_H */