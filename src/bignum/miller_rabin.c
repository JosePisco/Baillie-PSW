#include "crypto_utils.h"
#include "miller_rabin.h"

/* base 2 only for the Baillie-PSW test */
bool BN_miller_rabin_base_2(BIGNUM *n)
{
    bool result = false;

    BIGNUM *bn_two = BN_new();
    BIGNUM *bn_n_minus_one = BN_new();
    BIGNUM *s = BN_new();
    BIGNUM *x = BN_new();

    BN_CTX *ctx = BN_CTX_new();



    if (!BN_set_word(bn_two, 2)) goto done;

    /* specific cases */
    if (BN_is_one(n) || BN_is_word(n, 2) || BN_is_word(n, 3)) {
        result = true;
        goto done;
    }

    if (!BN_is_odd(n)) {
        result = false;
        goto done;
    }

    int r = 0;
    if (!BN_sub(bn_n_minus_one, n, BN_value_one())) goto done;
    if (!BN_copy(s, bn_n_minus_one)) goto done; // s = n - 1 : always even
    while (BN_is_odd(s) == 0) { // s % 2 == 0 -> i.e. while s is even
        r++;
        if (!BN_rshift1(s, s)) goto done; // s /= 2
    }

    /* In general cases, the Miller Rabin test has more bases */
    if (!BN_mod_exp(x, bn_two, s, n, ctx)) goto done; // x = 2^s mod n
    if (BN_is_one(x) || BN_cmp(x, bn_n_minus_one) == 0) { // if x == 1 or x == n - 1
        result = true;
        goto done;
    }

    for (int i = 0; i < r - 1; ++i) {
        if (!BN_mod_exp(x, x, bn_two, n, ctx)) goto done; // x = x^2 mod n
        if (BN_cmp(x, bn_n_minus_one) == 0) {
            result = true;
            goto done;
        }
    }

    // result can only be false here

    done:
        BN_free(bn_two);
        BN_free(bn_n_minus_one);
        BN_free(s);
        BN_free(x);
        BN_CTX_free(ctx);
        return result;
}