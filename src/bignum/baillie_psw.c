#include "baillie_psw.h"
#include "crypto_utils.h"
#include "miller_rabin.h"

/*
 * This is a stronger variant of the Lucas test in FIPS 186-4, Appendix C.3.3.
 * Every strong Lucas pseudoprime n is also a Lucas pseudoprime since
 * U_{n+1} == 0 follows from U_k == 0 or V_{k * 2^r} == 0 for 0 <= r < s.
 */
static int bn_strong_lucas_test(const BIGNUM *n, const BIGNUM *D, BN_CTX *ctx)
{
	BIGNUM *k, *U, *V;
	int r, s;
	int ret = -1;

	BN_CTX_start(ctx);

	if ((k = BN_CTX_get(ctx)) == NULL)
		goto done;
	if ((U = BN_CTX_get(ctx)) == NULL)
		goto done;
	if ((V = BN_CTX_get(ctx)) == NULL)
		goto done;

	/*
	 * Factorize n + 1 = k * 2^s with odd k: shift away the s trailing ones
	 * of n and set the lowest bit of the resulting number k.
	 */
	s = 0;
	while (BN_is_bit_set(n, s))
		s++;
	if (!BN_rshift(k, n, s))
		goto done;
	if (!BN_set_bit(k, 0))
		goto done;

	/*
	 * Calculate the Lucas terms U_k and V_k. If either of them is zero,
	 * then n is a strong Lucas pseudoprime.
	 */
	if (!bn_lucas(U, V, k, D, n, ctx))
		goto done;

	if (BN_is_zero(U) || BN_is_zero(V)) {
		ret = 1;
		goto done;
	}

	/*
	 * Check if any V_{k * d^r} is zero for 1 <= r < s.
	 */
	for (r = 1; r < s; r++) {
		if (!bn_lucas_step(U, V, 0, D, n, ctx))
			goto done;

		if (BN_is_zero(V)) {
			ret = 1;
			goto done;
		}
	}

	/* If we got here, n is definitely composite. */
	ret = 0;

 done:
	BN_CTX_end(ctx);

	return ret;
}

/*
 * Test n for primality using the strong Lucas test with Selfridge's
 * parameters. Returns 1 if n is prime or a strong Lucas-Selfridge
 * pseudoprime. Returns 0 if n is definitely composite.
 */
static int bn_strong_lucas_selfridge(const BIGNUM *n, BN_CTX *ctx)
{
	BIGNUM *D, *two;
	int jacobi_symbol, perfect_square, sign;
	int ret = -1;

	BN_CTX_start(ctx);

	/* If n is a perfect square, it is composite. */
	//if (!bn_is_square(&perfect_square, n, ctx))
	//	goto done;
	//if (perfect_square) {
	//	ret = 0;
	//	goto done;
	//}
	(void) perfect_square;

	/*
	 * Find the first element D in the sequence 5, -7, 9, -11, 13, ...
	 * such that Jacobi(D, n) = -1 (Selfridge's algorithm).
	 */
	if ((D = BN_CTX_get(ctx)) == NULL)
		goto done;
	if ((two = BN_CTX_get(ctx)) == NULL)
		goto done;

	sign = 1;
	if (!BN_set_word(D, 5))
		goto done;
	if (!BN_set_word(two, 2))
		goto done;

	while (1) {
		/* For odd n the Kronecker symbol computes the Jacobi symbol. */
		if ((jacobi_symbol = BN_kronecker(D, n, ctx)) == -2)
			goto done;

		/* We found the value for D. */
		if (jacobi_symbol == -1)
			break;

		/* n and D have prime factors in common. */
		if (jacobi_symbol == 0) {
			ret = 0;
			goto done;
		}

		/* Subtract or add 2 to follow the sequence described above. */
		sign = -sign;
		if (!BN_uadd(D, D, two))
			goto done;
		BN_set_negative(D, sign == -1);
	}

	ret = bn_strong_lucas_test(n, D, ctx);

 done:
	BN_CTX_end(ctx);

	return ret;
}

/*
 * The Baillie-Pomerance-Selfridge-Wagstaff algorithm combines a Miller-Rabin
 * test for base 2 with a Strong Lucas pseudoprime test.
 */
int bn_is_prime_bpsw(const BIGNUM *n, BN_CTX *in_ctx)
{
	BN_CTX *ctx = in_ctx;
	BN_ULONG mod;
	int i;
	int ret = -1;

	if (BN_is_word(n, 2)) {
		ret = 1;
		goto done;
	}

	if (BN_cmp(n, BN_value_one()) <= 0 || !BN_is_odd(n)) {
		ret = 0;
		goto done;
	}

	/* Trial divisions with the first 2048 primes. */
	for (i = 0; i < NUMPRIMES; i++) {
		if ((mod = BN_mod_word(n, primes[i])) == (BN_ULONG)-1)
			goto done;
		if (mod == 0) {
			ret = BN_is_word(n, primes[i]);
			goto done;
		}
	}

	if (ctx == NULL)
		ctx = BN_CTX_new();
	if (ctx == NULL)
		goto done;

	if ((ret = bn_miller_rabin_base_2(n, ctx)) <= 0)
		goto done;

	/* XXX - Miller-Rabin for random bases? - see FIPS 186-4, Table C.1. */

	ret = bn_strong_lucas_selfridge(n, ctx);

 done:
	if (ctx != in_ctx)
		BN_CTX_free(ctx);

	return ret;
}
