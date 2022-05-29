#include "crypto_utils.h"

int jacobi_symbol(BIGNUM *d, BIGNUM *n)
{
    int result = 1;

    BIGNUM *bn_four = BN_new();
    BIGNUM *bn_eight = BN_new();
    BIGNUM *d_mod_four = BN_new();
    BIGNUM *n_mod_four = BN_new();
    BIGNUM *r = BN_new();

    BN_CTX *ctx = BN_CTX_new();

    if (!BN_set_word(bn_four, 4)) goto done;
    if (!BN_set_word(bn_eight, 8)) goto done;

    if (!BN_mod(d, d, n, ctx)) { // d %= n
        result = 0;
        goto done;
    }

    while (!BN_is_zero(d)) { // d != 0
        while (!BN_is_odd(d)) { //d % 2 == 0
            if (!BN_rshift1(d, d)) goto done; // d /= 2
            if (!BN_mod(r, n, bn_eight, ctx)) goto done; // r = n % 8
            if (BN_is_word(r, 3) || BN_is_word(r, 5)) // if r == 3 or r == 5
                result = -result;
        }

        BN_swap(d, n);

        if (!BN_mod(d_mod_four, d, bn_four, ctx)) goto done;
        if (!BN_mod(n_mod_four, n, bn_four, ctx)) goto done;

        if (BN_is_word(d_mod_four, 3) && BN_is_word(n_mod_four, 3)) // if d % 4 == 3 and n % 4 == 3
            result = -result;

        if (!BN_mod(d, d, n, ctx)) { // d %= n
            result = 0;
            goto done;
        }
    }

    if (BN_is_one(n)) // if n == 1
        goto done;

    result = 0;

    done:
        BN_free(bn_four);
        BN_free(bn_eight);
        BN_free(d_mod_four);
        BN_free(n_mod_four);
        BN_free(r);
        BN_CTX_free(ctx);
        return result;
}

int lucas(BIGNUM *k, BIGNUM *D, BIGNUM *P, BIGNUM *n, struct lucas_sequence *lucas_)
{
    int result = 0;

    BIGNUM *bn_two = BN_new();
    BIGNUM *U = BN_new();
    BIGNUM *tmp_U = BN_new();
    BIGNUM *tmp_U_sqr = BN_new();
    BIGNUM *V = BN_dup(P);
    BIGNUM *V_sqr = BN_new();
    BIGNUM *U_V_add = BN_new();
    BIGNUM *P_V = BN_new();
    BIGNUM *D_mul_tmp_U = BN_new();
    BIGNUM *D_mul_tmp_U_sqr = BN_new();
    BIGNUM *numerator = BN_new();
    BIGNUM *k_obj = BN_dup(k); // objective to reach within the end of the computations
    BIGNUM *modinv_two = BN_new();

    BN_CTX *ctx = BN_CTX_new();

    if (!BN_set_word(U, 1)) goto done;
    if (!BN_set_word(k, 1)) goto done; // set k = 1 then increment it to check if we reach k_objective
    if (!BN_set_word(bn_two, 2)) goto done;
    if (!BN_mod_inverse(modinv_two, bn_two, n, ctx)) goto done; // usefull when we will need to divide by 2

    int bitlength = BN_num_bits(k_obj);
    /* starts at bit 2, bit one is done by setting U=1 and V=P */
    int i = 2;
    /* Iterates over the bits of k from left to right */
    while (i != bitlength + 1)
    {
        /* get the ith bit (from left being 0 to right) of a number*/
        int bit = BN_is_bit_set(k_obj, bitlength - i);

        tmp_U = BN_dup(U);
        if (!BN_mod_mul(U, U, V, n, ctx)) goto done; // U = U * V % n

        /* V = (V*V + D * tmp_U*tmp_U) * modinv(2, n) % n */
        if (!BN_sqr(V_sqr, V, ctx)) goto done; // V*V
        if (!BN_sqr(tmp_U_sqr, tmp_U, ctx)) goto done; // tmp_U * tmp_U
        if (!BN_mul(D_mul_tmp_U_sqr, D, tmp_U_sqr, ctx)) goto done; // D * tmp_U^2
        if (!BN_add(numerator, V_sqr, D_mul_tmp_U_sqr)) goto done; // (V*V + D * tmp_U*tmp_U)

        /* A divison by two is equivalent to a multiplication by the mod inverse of 2 by n*/
        if (!BN_mod_mul(V, numerator, modinv_two, n, ctx)) goto done; // V = (V*V + D * tmp_U*tmp_U) * modinv(2, n) % n
        if (!BN_lshift1(k, k)) goto done; // k *= 2

        if (bit == 1) {
            tmp_U = BN_dup(U);
            if (!BN_add(U_V_add, U, V)) goto done; // U + V
            if (!BN_mul(numerator, P, U_V_add, ctx)) goto done; // P * U+V
            if (!BN_mod_mul(U, numerator, modinv_two, n, ctx)) goto done; // U = (P * U + V) * modinv(2, n) % n


            // V = (D * tmp_U + P * V) * modinv(2, n) % n;
            if (!BN_mul(D_mul_tmp_U, D, tmp_U, ctx)) goto done;
            if (!BN_mul(P_V, P, V, ctx)) goto done;
            if (!BN_add(numerator, D_mul_tmp_U, P_V)) goto done;
            if (!BN_mod_mul(V, numerator, modinv_two, n, ctx)) goto done;

            if (!BN_add(k, k, BN_value_one())) goto done; // k += 1
        }
        i++;
    }

    if (BN_cmp(k, k_obj) != 0) {
        result = 0;
        goto done;
    }

    lucas_->U = BN_dup(U);
    lucas_->V = BN_dup(V);
    result = 1;

    done:
        BN_free(bn_two);
        BN_free(U);
        BN_free(tmp_U);
        BN_free(tmp_U_sqr);
        BN_free(V);
        BN_free(V_sqr);
        BN_free(U_V_add);
        BN_free(P_V);
        BN_free(D_mul_tmp_U);
        BN_free(D_mul_tmp_U_sqr);
        BN_free(numerator);
        BN_free(k_obj);
        BN_free(modinv_two);
        BN_CTX_free(ctx);
        return result;
}
