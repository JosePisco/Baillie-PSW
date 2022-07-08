#include "crypto_utils.h"

/*
 * For an odd n compute a / 2 (mod n). If a is even, we can do a plain
 * division, otherwise calculate (a + n) / 2. Then reduce (mod n).
 */
static int bn_division_by_two_mod_n(BIGNUM *r, BIGNUM *a, const BIGNUM *n,
	BN_CTX *ctx)
{
	if (!BN_is_odd(n))
		return 0;

	if (!BN_mod(r, a, n, ctx))
		return 0;

	if (BN_is_odd(r)) {
		if (!BN_add(r, r, n))
			return 0;
	}

	if (!BN_rshift1(r, r))
		return 0;

	return 1;
}

/*
 * Given the next binary digit of k and the current Lucas terms U and V, this
 * helper computes the next terms in the Lucas sequence defined as follows:
 *
 *   U' = U * V                  (mod n)
 *   V' = (V^2 + D * U^2) / 2    (mod n)
 *
 * If digit == 0, bn_lucas_step() returns U' and V'. If digit == 1, it returns
 *
 *   U'' = (U' + V') / 2         (mod n)
 *   V'' = (V' + D * U') / 2     (mod n)
 *
 * Compare with FIPS 186-4, Appendix C.3.3, step 6.
 */
int bn_lucas_step(BIGNUM *U, BIGNUM *V, int digit, const BIGNUM *D,
    const BIGNUM *n, BN_CTX *ctx)
{
	BIGNUM *tmp;
	int ret = 0;

	BN_CTX_start(ctx);

	if ((tmp = BN_CTX_get(ctx)) == NULL)
		goto done;

	/* Store D * U^2 before computing U'. */
	if (!BN_sqr(tmp, U, ctx))
		goto done;
	if (!BN_mul(tmp, D, tmp, ctx))
		goto done;

	/* U' = U * V (mod n). */
	if (!BN_mod_mul(U, U, V, n, ctx))
		goto done;

	/* V' = (V^2 + D * U^2) / 2 (mod n). */
	if (!BN_sqr(V, V, ctx))
		goto done;
	if (!BN_add(V, V, tmp))
		goto done;
	if (!bn_division_by_two_mod_n(V, V, n, ctx))
		goto done;

	if (digit == 1) {
		/* Store D * U' before computing U''. */
		if (!BN_mul(tmp, D, U, ctx))
			goto done;

		/* U'' = (U' + V') / 2 (mod n). */
		if (!BN_add(U, U, V))
			goto done;
		if (!bn_division_by_two_mod_n(U, U, n, ctx))
			goto done;

		/* V'' = (V' + D * U') / 2 (mod n). */
		if (!BN_add(V, V, tmp))
			goto done;
		if (!bn_division_by_two_mod_n(V, V, n, ctx))
			goto done;
	}

	ret = 1;

 done:
	BN_CTX_end(ctx);

	return ret;
}

/*
 * Compute the Lucas terms U_k, V_k, see FIPS 186-4, Appendix C.3.3, steps 4-6.
 */
int bn_lucas(BIGNUM *U, BIGNUM *V, const BIGNUM *k, const BIGNUM *D,
    const BIGNUM *n, BN_CTX *ctx)
{
	int digit, i;
	int ret = 0;

	if (!BN_one(U))
		goto done;
	if (!BN_one(V))
		goto done;

	/*
	 * Iterate over the digits of k from MSB to LSB. Start at digit 2
	 * since the first digit is dealt with by setting U = 1 and V = 1.
	 */
	for (i = BN_num_bits(k) - 2; i >= 0; i--) {
		digit = BN_is_bit_set(k, i);

		if (!bn_lucas_step(U, V, digit, D, n, ctx))
			goto done;
	}

	ret = 1;

 done:
	return ret;
}
