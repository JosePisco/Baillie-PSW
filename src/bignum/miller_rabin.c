#include "crypto_utils.h"
#include "miller_rabin.h"

int bn_miller_rabin_base_2(const BIGNUM *n, BN_CTX *ctx)
{
	BIGNUM *n_minus_one, *k, *x;
	int i, s;
	int ret = -1;

	BN_CTX_start(ctx);

	if ((n_minus_one = BN_CTX_get(ctx)) == NULL)
		goto done;
	if ((k = BN_CTX_get(ctx)) == NULL)
		goto done;
	if ((x = BN_CTX_get(ctx)) == NULL)
		goto done;

	if (BN_is_word(n, 2) || BN_is_word(n, 3)) {
		ret = 1;
		goto done;
	}

	if (BN_cmp(n, BN_value_one()) == 0 || !BN_is_odd(n)) {
		ret = 0;
		goto done;
	}

	if (!BN_sub(n_minus_one, n, BN_value_one()))
		goto done;

	s = 0;
	while (!BN_is_bit_set(n_minus_one, s))
		s++;
	if (!BN_rshift(k, n_minus_one, s))
		goto done;

	/* If 2^k is 1 or -1 (mod n) then n is a 2-pseudoprime. */
	if (!BN_set_word(x, 2))
		goto done;
	if (!BN_mod_exp(x, x, k, n, ctx))
		goto done;

	if (BN_is_one(x) || BN_cmp(x, n_minus_one) == 0) {
		ret = 1;
		goto done;
	}

	for (i = 1; i < s; i++) {
		if (!BN_mod_sqr(x, x, n, ctx))
			goto done;
		if (BN_cmp(x, n_minus_one) == 0) {
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
