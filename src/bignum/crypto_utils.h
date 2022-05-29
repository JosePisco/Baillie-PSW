#ifndef CRYPTO_UTILS_H
#define CRYPTO_UTILS_H

#include <openssl/bn.h>

/* Struct containing the values of the Lucas sequences for the d_th term */
struct lucas_sequence {
    BIGNUM *U;
    BIGNUM *V;
};

/* PROTOTYPES */

/* Computes the Jacobi Symbol for any odd natural n and natural d*/
int jacobi_symbol(BIGNUM *d, BIGNUM *n);

/* Computes the u_k and v_k term of the Lucas sequences defined at
https://en.wikipedia.org/wiki/Lucas_pseudoprime#Strong_Lucas_pseudoprimes */
int lucas(BIGNUM *k, BIGNUM *D, BIGNUM *P, BIGNUM *n, struct lucas_sequence *lucas_);


#endif /* CRYPTO_UTILS_H */