#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <openssl/bn.h>

/* PROTOTYPES */
int bn_miller_rabin_base_2(const BIGNUM *n, BN_CTX *ctx);

#endif /* MILLER_RABIN_H */