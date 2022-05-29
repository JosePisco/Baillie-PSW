#ifndef CRYPTO_UTILS_H
#define CRYPTO_UTILS_H

#include <stdint.h>

/* Struct containing the values of the Lucas sequences for the d_th term */
struct lucas_sequence {
    unsigned long long U;
    unsigned long long V;
};

/* PROTOTYPES */

/* Computes the Jacobi Symbol for any odd natural n and natural d*/
int jacobi_symbol(uint64_t d, uint64_t n);

/* Computes the u_k and v_k term of the Lucas sequences defined at
https://en.wikipedia.org/wiki/Lucas_pseudoprime#Strong_Lucas_pseudoprimes */
struct lucas_sequence lucas(unsigned long long k, unsigned long long D, unsigned long long P, unsigned long long n);

/* Computes the inverse of x modulo n*/
int modinv(int x, int n);

/* computes the modular exponentiation of x by y mod p i.e. x^y mod p*/
int modexp(int x, int y, int p);

#endif /* CRYPTO_UTILS_H */