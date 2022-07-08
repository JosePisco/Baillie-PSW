#ifndef BAILLIE_PSW_H
#define BAILLIE_PSW_H

#define NUM_FIRST_PRIMES 669

#include <openssl/bn.h>

/* PROTOTYPES */
int bn_is_prime_bpsw(const BIGNUM *n, BN_CTX *in_ctx);

#endif /* BAILLIE_PSW_H */