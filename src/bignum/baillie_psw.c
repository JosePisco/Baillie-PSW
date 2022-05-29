#include "baillie_psw.h"
#include "crypto_utils.h"
#include "miller_rabin.h"

/* First primes to 5000, can save some time to check if the number is divisible by one of them */
static const int first_primes_5000[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999};

bool BN_strong_lucas_test(BIGNUM *n, BIGNUM *D, BIGNUM *P)
{
    bool result = false;

    BIGNUM *bn_two = BN_new();
    BIGNUM *d = BN_new();
    BIGNUM *d_shift = BN_new();
    BIGNUM *bn_s = BN_new();
    BIGNUM *bn_r = BN_new();
    BIGNUM *n_plus_one = BN_new();
    BIGNUM *two_pow = BN_new();
    BIGNUM *U_d_mod_n = BN_new();
    BIGNUM *V_d_mod_n = BN_new();

    BN_CTX *ctx = BN_CTX_new();

    /* object in which to store the dth value of the lucas sequence for U and V */
    struct lucas_sequence lucas_ = {NULL, NULL};
    lucas_.U = BN_new();
    lucas_.V = BN_new();

    BN_zero_ex(d); // d = 0;
    if (!BN_set_word(bn_two, 2)) goto done;
    if (!BN_add(n_plus_one, n, BN_value_one())) goto done;

    /* Find d, s such as n+1 = d * 2^s
    d is odd, n+1 even */
    int s = 1;
    while (true) {
        if (!BN_set_word(bn_s, s)) goto done;
        if (!BN_exp(two_pow, bn_two, bn_s, ctx)) goto done;

        // if (!BN_div(d, NULL, n_plus_one, two_pow, ctx)) goto done; // d = (n+1) / 2^s
        if (!BN_rshift(d, n_plus_one, s)) goto done; // d = (n+1) / 2^s -> right shift 's' times
        if (BN_is_odd(d)) // if d % 2 == 1
            break;

        s++;
    }
    /*
    d is odd
    using the lucas sequences formulas:
    if U_d mod n = 0 or V_(d*2^r) mod n = O
    then n is a strong lucas pspr */

    if (!lucas(d, D, P, n, &lucas_)) goto done;

    if (!BN_mod(U_d_mod_n, lucas_.U, n, ctx)) goto done;
    if (BN_is_zero(U_d_mod_n)) { // U_d % n == 0
        result = true;
        goto done;
    }
    else {
        int r;
        for (r = 0; r < s; ++r) {
            if (!BN_set_word(bn_r, r)) goto done;

            if (!BN_exp(two_pow, bn_two, bn_r, ctx)) goto done; // 2^r
            if (!BN_mul(d_shift, d, two_pow, ctx)) goto done; // d * 2^r
            if (!lucas(d_shift, D, P, n, &lucas_)) goto done;

            if (!BN_mod(V_d_mod_n, lucas_.V, n, ctx)) goto done;
            if (BN_is_zero(V_d_mod_n)) { // V_d % n == 0
                result = true;
                goto done;
            }
        }
    }

    result = false;

    done:
        BN_free(bn_two);
        BN_free(d);
        BN_free(d_shift);
        BN_free(bn_s);
        BN_free(bn_r);
        BN_free(n_plus_one);
        BN_free(two_pow);
        BN_free(lucas_.U);
        BN_free(lucas_.V);
        BN_free(U_d_mod_n);
        BN_free(V_d_mod_n);
        BN_CTX_free(ctx);
        return result;
}

bool BN_strong_lucas_selfridge(BIGNUM *n)
{
    bool result = false;
    /*
    we also need to filter out all perfect square values of N
    this is because we will later require an integer D for which Jacobi(D,N) = -1
    no such integer exists if N is a perfect square.*/
    //BIGNUM *root = BN_new();

    BIGNUM *bn_negative_one = BN_new();
    BIGNUM *bn_zero = BN_new();
    BIGNUM *bn_two = BN_new();
    BIGNUM *sign = BN_new();
    BIGNUM *D_abs = BN_new();
    BIGNUM *tmp_D = BN_new();
    BIGNUM *D = BN_new();
    BIGNUM *P = BN_new();

    BN_CTX *ctx = BN_CTX_new();
    // bn_sqrt does not exist
    //if (!BN_sqrt(root, n, ctx)) goto done; // root = sqrt(n)
    //if (!BN_mul(root, root, root)) goto done; // root = root*root -> test if n is a squared number

    //if (BN_cmp(root , n) == 0)
    //    return false; // free everything

    /*
    Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
    such that Jacobi(D,N) = -1 (Selfridge's algorithm). Theory indicates
    that, if N is not a perfect square, D will "nearly always" be "small."
    */

    if (!BN_dec2bn(&bn_negative_one, "-1")) goto done; // use dec2bn and not set_word for negatives
    BN_zero_ex(bn_zero);
    if (!BN_set_word(bn_two, 2)) goto done;
    if (!BN_set_word(D_abs, 5)) goto done; // D_abs = 5
    if (!BN_set_word(sign, 1)) goto done;

    if (!BN_mul(D, D_abs, sign, ctx)) goto done; // D = d_abs * sign
    //int i =0;
    while (true) {

        // if D is negative, take the inv mod n of D, to compute modulos efficiently
        tmp_D = BN_dup(D); // tmp_D = D

        if (BN_cmp(D, bn_zero) < 0)
            if (!BN_mod_inverse(tmp_D, D, n, ctx)) goto done; // if D < 0, take modinv of D

        if (jacobi_symbol(BN_dup(tmp_D), BN_dup(n)) == -1)
            break;


        if (!BN_add(D_abs, D_abs, bn_two)) goto done; // d_abs += 2;
        if (!BN_mul(sign, sign, bn_negative_one, ctx)) goto done; //sign = sign * -1 <=> sign = -sign
        if (!BN_mul(D, D_abs, sign, ctx)) goto done; // D = d_abs * sign for Selfridge sequence
    }

    // Selfridge reference
    if (!BN_set_word(P, 1)) goto done; // P = 1
    if (BN_cmp(D, bn_zero) < 0)
        if (!BN_mod_inverse(D, D, n, ctx)) goto done; // if D < 0, take modinv of D

    result = BN_strong_lucas_test(n, D, P);

    done:
        BN_free(bn_negative_one);
        BN_free(bn_zero);
        BN_free(bn_two);
        BN_free(sign);
        BN_free(D_abs);
        BN_free(tmp_D);
        BN_free(D);
        BN_free(P);
        BN_CTX_free(ctx);
        return result;
}

bool BN_baillie_psw(BIGNUM *n)
{
    bool result = false;

    BIGNUM *bn_two = BN_new();
    BIGNUM *n_mod_prime = BN_new();
    BIGNUM *prime = BN_new();

    BN_CTX *ctx = BN_CTX_new();

    if (!BN_set_word(bn_two, 2)) goto done;

    /* Specific cases */
    if (BN_cmp(n, BN_value_one()) <= 0) { // n <= 1
        result = false;
        goto done;
    }

    /* More specific cases */
    if (BN_is_word(n, 2) || BN_is_word(n, 11)) { // n == 2 or n == 11
        result =  true;
        goto done;
    }

    if (!BN_is_odd(n)) {
        result = false;
        goto done;
    }

    /* there are cases where this may save some time
    test with first primes < 5000 if they divide n */
    int i;
    for (i = 0; i < NUM_FIRST_PRIMES; ++i) {
        if (!BN_set_word(prime, first_primes_5000[i])) goto done;
        if (!BN_mod(n_mod_prime, n, prime, ctx)) goto done;

        if (BN_is_zero(n_mod_prime) && BN_cmp(n, prime) != 0) { // if n % prime == 0 and n != prime
            result = false;
            goto done;
        }
    }

    if (BN_miller_rabin_base_2(n) == false) {
        result = false;
        goto done;
    }

    result = BN_strong_lucas_selfridge(n);

    done:
        BN_free(bn_two);
        BN_free(n_mod_prime);
        BN_free(prime);
        BN_CTX_free(ctx);
        return result;
}
