#include <stdio.h>

#include "baillie_psw.h"
#include "tests_bpsw.h"

static const int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999};

void bpsw_tests(void)
{
	/* Complete BPSW tests */
	BN_CTX *ctx;

	if ((ctx = BN_CTX_new()) == NULL)
		goto done;
	BIGNUM *bn_tmp = BN_new();
	for (size_t i = 0; i < NUMBER_PRIMES; ++i)
	{
		if (!BN_set_word(bn_tmp, primes[i]))
			goto done;
		if (bn_is_prime_bpsw(bn_tmp, ctx) == 0)
			printf("WRONG: %d\n", primes[i]);
	}

	if (!BN_set_word(bn_tmp, 2297))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("bpsw : 2297 : ERROR\n");

	if (!BN_set_word(bn_tmp, 983))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("bpsw : 983 : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "64713796208369358001824144539576100181603399893984043952483371071999434497213352403525560173357970261557149261464469317297080737722690627775200014918575930145441973136656843904474027376966799110163773563217"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : arnault : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "5494192362793378419461850103325205105746618253649659755062958621743505179898846930508801354093259555029072978427447589979886135102666157147211855994196065010559627866228837393549736942388413848102819329"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : carmichael : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "16675213704264635359855338582319826599139881187670369606766303648011276513669171437386545724558172168833854673196998468433491303861872481723479169721198573611218928381296412105177684568788527085598972541364147607691655006142778451721220612611799188187750770557721259009955396543567551776802396501325432794976381476775210723134569779554780443295670451640271953769851816031211776828994779252013021094159247355897769604124730071091196654314171749454933025029709544871025024355039907165438184621004838120586980524265072507369438134097822899220153421895988975193331306413676318626385465709286848580363599123478758994747523"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("bpsw : prime 2048 : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "959491208865988256557944966011457536510073817351083162968895417886879744503593692538559369785926840995378444445068108561485601929143510010358069621544003641121478373025752564910898422651413959867616347939213554376090153588726637847022636651202378739187353368559972269295289980182299111431639915236001505106822802768759633474923592003419399920843283963565383193262584302086666533111280974765713358669856702801311198938154325336500546289645515902016579504294225219858609971900580652223730938989668005796665039248963275153507225290607735305406251201088857571815609751681595888705568164910938567846875087593099968524178447965007193035298904760923512762129228271256363692633973298393722599956171716723569356545877488813905083588202963374045233262983749952957463015682959958399696232549511957193070335700470463733015865612049969622208408989913732039620899566196270027025096872535077620197865090356051206814157620286824460589838424452408284975863397145679546375208004975463648911264702462933126328885098164042009448106744928982782501900216013387188560055250161101023705573603132300101630570572572736735233098849674804944795220327468365759516613158470220926527181876078187374281123697114208609629101800048460373829711815457106181195590025033"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("bpsw : prime 4096 : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "190316680009"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : square prime : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "119185643267526116681929"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : square number : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "0f2ad9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : square : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "23511795499245521367381353907826033845413959512754234716909319582718218990844601828440516176385782657177820467113582727425903757512093921812780923509768534574418430699968761836213126769338525755896558788188167211621848801007796153305296777344507158492064242356537308601910305816296760223665287736039173089673224829099442478413679846949501396263075692739679818833943744599023772829448431021942878456136500100143870231595650256912336403410042784597625682408395500502077495815474338012307807527414993090562152792185812293254881449377332845024690551437267399667063015005747883546451693295669354712269552726710363944738481"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : big square prime : ERROR\n");

	if (!BN_dec2bn(&bn_tmp, "222222"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : even number 222222\n");

	if (!BN_dec2bn(&bn_tmp, "178368728938729712"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bpsw : another even number (178368728938729712)\n");

	if (!BN_hex2bn(&bn_tmp, "ff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("small non prime integer\n");

	if (!BN_hex2bn(&bn_tmp, "00"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("small non prime integer\n");

	if (!BN_hex2bn(&bn_tmp, "01"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("small non prime integer\n");

	if (!BN_hex2bn(&bn_tmp, "07ffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Non-prime Mersenne number that is pseudoprime to base 2\n");

	if (!BN_hex2bn(&bn_tmp, "7fffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Non-prime Mersenne number that is pseudoprime to base 2\n");

	if (!BN_hex2bn(&bn_tmp, "0100000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Non-prime Fermat number\n");

	if (!BN_hex2bn(&bn_tmp, "010000000000000000000000000000000000000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Non-prime Fermat number\n");

	if (!BN_hex2bn(&bn_tmp, "0100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Non-prime Fermat number\n");

	if (!BN_hex2bn(&bn_tmp, "123a99"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("pseudoprime square derived from Wiefrich prime\n");

	if (!BN_hex2bn(&bn_tmp, "00bc18d1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("pseudoprime square derived from Wiefrich prime\n");

	if (!BN_hex2bn(&bn_tmp, "04"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("square\n");

	if (!BN_hex2bn(&bn_tmp, "09"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("square\n");

	if (!BN_hex2bn(&bn_tmp, "010201"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("square\n");

	if (!BN_hex2bn(&bn_tmp, "0f2ad9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("square\n");

	if (!BN_hex2bn(&bn_tmp, "01f51f3fee3b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "032907381cdf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "0136a352b2c8c1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "023c3db80e80e53bd1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "0504e8e504fd585e79193ca1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "00b7d84161830e3f6f2231a7a1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("G. Jaeschke: On strong pseudoprimes to several bases, Math o. comp. v.61, p 915-926\n");

	if (!BN_hex2bn(&bn_tmp, "4c6092d9a7a5462b34e5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "22c9a603ee84bb9c4cad"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "437ae92817f9fc85b7e5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "0190e262098f0d746505"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "027a5f7ca7b29ee74d5525"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "008d60a89f3f36cb1fd495"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "02be6951adc5b22410a5fd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "0292a0068ebb0ed3251f55"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "750b703e68cb957ab415"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "02d0facc78aeeb89f5b299"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprimes to 12 or more bases from https://arxiv.org/pdf/1509.00864v1.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "09bdc1c98b9b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Maple\n");

	if (!BN_hex2bn(&bn_tmp, "0ffb48c934842b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Maple\n");

	if (!BN_hex2bn(&bn_tmp, "18444fdb12afb7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Maple\n");

	if (!BN_hex2bn(&bn_tmp, "08e4f37e51"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Mathematica 2.0\n");

	if (!BN_hex2bn(&bn_tmp, "179d55b600e7f1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Mathematica 2.0\n");

	if (!BN_hex2bn(&bn_tmp, "085270bd76a142abc3037d1aab3b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Axioms primality test\n");

	if (!BN_hex2bn(&bn_tmp, "02cb78fe3f36c4f5f05dbe92b82798d5fc18f2bfaaa388ef"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Some primality testing algorithms a counter example for Axioms primality test\n");

	if (!BN_hex2bn(&bn_tmp, "4682f52f0b54308d315b2fbec25065506c77be95912b137bc6eecffad8a299b631c55ce068702b1b3e4ce50958994c289b148fb298a8c603a0959cb0ba5ad4bcba278cf4c87e0ff85a62a25c40849662c53d0f81cf9e4431d8c391586629260e558db473997db20108278b1ae374089140d93bc2c5a808ad3aaf212f60bfc93cc0c788149dcd82f7ab"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A composite q that was acceptied by Gnu Crypto. http://www.iacr.org/archive/pkc2005/33860010/33860010.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "00f67307e54779cfe9120bf862afc5466c5d6d0783d12df5215c0c981c51e4bfc098e9afd574f51b18c820259b692ec0bf7c9d6e56e9bb99fbd3b7ecc4082146a9d7a5b7bc6519d476c4a9975d9c3e3b12bee45b7accb07a6a68ea583ac2523ef32ee6d01bc766b59c43031f9c6980c9b4317da6825be9f7c5db03283d04c13323"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00c1d00b32d63e3ea4fb69ab6b9dee40a17fada46c122e52a53fecd3fe613303f51c07871dc0b5d8d8c1705b484de6bdb7f442efecd7d9f59dc36e495f72905c7619bc4d3706283774e704a3adad7d6c1be42ddeffc2ca5b1c0e31b58ed606f16dc14676e60ecff42ae33e503621e232ba449e91e3a9909e80a8318610aea3b7cf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01c2faadef91d43c9ab1320020e08e2ec3c34012bd0db94a1175170dc5aec26897e867d0b7a7273119fbe1115f02875b522566016f69f319ad5485e7458fcf50205d22ba765cc586a6037be987b6832c46227df19cd8ce0641794b60b73fbdd3c104870ae9bdf0194e772c985536e860b90b7fa3eb205af6b224413f5813836abb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0109fddd44575367466c67aaa921047b367515c9aa579eb60728034ad2d56f10eb01cfadb3ba0abde99f348bc3c70559bc24551b85937ca4c886abc0826cc1c310f14393652c1b4994953881bd2d81de0f2a280839829543f429bc41bf3c6db120bb150173e2707f36d1f76318249851f4fedc39e36aaaca48686de03e6d256973"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00ffd0847cdda5a4fdfd2345bc731f1bc77843478950d33b2830ef0caf8deffdbe6309fe61fb67dded6659e433f30363339dbcc7c0832593f33c24a8b8f0e28038cb6edeed58ae765e6884ac0b66b5218cc758e6247269d24be9f91865d33c105219ffbce00c6c2d6391448643bcf5138268f510258f638b90a6c8b53bfc121759"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0118d077827c6db85bc61d53063edf5676d6ac65b611d836eed07ee7e1d15c02d999a3eb78ce662edaf457f0f7d9c0a0305acc1faec4170400f0610a797de50ebfb08fd0a5da77144a1e0236e2bc6d8d2a6a719e59df071367cd61275f372e23b1c0187d87d15bda5f71f4705b1c3aaaa8ad951d20cee93274b151f3f9a55bd693"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01c09377e15f53b1329b6e8a08bf0f94da27dd29c89be74544d705173a0bdd410935e186dd95ac113732674fe08585690ebe9f749a116a8c64e1b4a281ef0cb28bc70b1639bc1352ff5777783bd72e3b8495c1494ae11fb32bdaba8c80870a3de71c0c27f07983e97500c0ec0321b86c679c53ae7f8c76ddbf6a9cc3ff63e45023"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00f35cac3bb3c7cf5e4e50162f4ca889ac7b875f4aac08c5a2433600e9bc64db6c9895aaccf3ee98783ee2cfd8a5e448b265bbc4cda6cb80d487c7967d5a6724fae1ffd27c70f579e62b49f29819c6221d7659fa9364e8e37795d88611506b552a20533f1f6446a35b41a986d304fdd7a39f484331b4fbf242f95b80788cff39cd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01e9df6f069f5984c080087127f90437f2d38f19385b3592d17a5f23603ec6315c36a88d2012e85eca62a983de7ef27673c605155b5647311840cf8887be8267fbc01cec3f7e0467d5e9a812e5dca577cc8ac93971c84f8cea94637c60c0bfe5d7f4b4f950e60ad077941190afaa905d6d5d570c9b4dab98c32c7abc42346f894d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00c5286502dda772fc22d43b0a2f46823777a91f580f3a1261c47be8e2010a5ad9395e2c036b32813dcdaad33c8f2f4a522593e31ae55ef05c8df8ed58636ac1b9db2b205797d39343e0868ff02bef46d18736bedc6f527730da8594d45d0447e7c7f0e8ca12b285b88aea5e343264874ac22038f5821bd96519d49caf45184f97"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01c29fe8b7e63795218563774685b9fe85eada73691a6420c38f0e9f2f802e89c77ae78716924e4efb5e4c639ca98ddb0c9e35cbc6313196b3327672527404b6da8ff7813915702fb7fa254c1cdc167a34170da57606ccff876ca0ce5e920f443e389fc9d0c071b908c6675b6a9f5903d6d22ad490e6476a7e13adcaf988663b3b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01e8648f8abce82efb0afa9861c96c428f690c5fe33b9c9d47f97198542c982e607fd9700f876159ea404983f4eecbaf2a73b262085da4b7b5de8f6e8ca0b712f5e89c0e8f024033879f858f814275a3ea5543fd539e74f5e099769d0d726ebd8bc74bda6e2f8ffabbb7d043f7818cd8d531180a827731fac59f45b2af35d273f9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00cedc5db464312d6f1ecf53a40bde07ae0d5540ef75a4802ff469142270049dbba2b74e4ece7340d8eb99bac1a3d6f0b52ebb41794d3cd4e4a588431879ff81818abc50bca5e686a06d48461b425be62d3c064321429e346960163f897d21b362dc72f306a6865cfb9c8c5682cc7fcd7dc6ac4202e8d070729ef9e3b526236c71"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0093ec9e6681f1bc1d6534add99d97e0d907828996bb3d7b481f3ceaefbe8f3fdf15698302ce26feb84c08994079c9f368af8171faf76801fe6dfdaecd587fa0edc751d64ff7e9aa73fb7aa51a8469379bac38e9d7941e0bbdcf658633daea40738e81f5605198b04fe8fd49646da4e98c2282a8041c25bb9894252412472294f9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0129fffd0bf1827f2847f45bd490d5423f67d87eb8254535d57078707e19f2ca5ca10602c5eca552fbdc77e30592b7498254f901cad02e0bf59802f5582cbb3059a1979a5e5311855807b1cbeff86a651dbf3818c3b6cf50092c9b744c4831873d1d0d8c23f23b39517ce435a257e5026cfa0be280672e1bba3074b2cdc6474a37"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "017232b942eedc8a0df14f5c1ad4e099f192b242b7d3dff09c50cecfe636c72c6c8ba1c65dde4396282e1a1c823b6d5d9c0c9068b39e202dcba26a9d35a00b7bb6bede272820fbbba503bc1866c6ae183d8b50e28555a921121929862ce87ea4ddde8f9d6ff2e17a8ee7cf9d306faa0815a4d46e8dfd4b7ea538b7399cc1c06c1f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00d3768b43c242fa7ac1de856dc7bd413b79d544bb8d38677bc9f44aa116ac5525c3e7fcf2fb2c1d3de61844931f47646b4c5f7de226031c925acbe57f1cd292fec7e7d4fd25afa128704ffd8da910ef18961e081e88d40bc37582b087f1b1f39fe4d23a03ec6b869c76fa3aed7a3606c469069c4fa1d4ff1c6112da16ba9dcf97"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "011b5119e5c68a710158c36d414597b4e1ccff332d1b437a4d2da2d2269ad2b626fde79e3ba7ed92128e5feaa87556f18ca6937b5a88f4738608d6bb6aacaf4fb719d67561d66dba9690009bcdbea2db4ee48d575722cbafbf1e487bab1c62ba0cde30a34620c7733b3e13d8b27fa035115680fb81016d1ca777b8a2bb7c399a47"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "008e9ee596ea83d06e1a9a4c3b75fc67f3c01de737be4dcdc18f1d10e322df48e455546ac8ac810129dbcb0fbf568987033cadef9d051f6032c8dca2804fc8d8d6e79f5d767963e4b6d72ac29d98d2520c29c8e69ffa59164d6a1e4cb55b7fcc60c7cb274da264203839873ec2f85f4ae377eeb6189e031b17e8603a01ef877b3f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00879d1e0bc0538cd9025110cec61a034305c8fdea2b9709ba80b0c45891e7ffc69c05285f4680b95b5882ad04210342314d3ab465ee1209d0690613a09bf7df0d48de18a7200e09e8b7944e748413ad64057fee2daacd099dcbb19920429cf9776d939c27c74c3adc8c41f1001f98d5293e018b1dde228abc6e79092331804bdb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00a14d02b57eb643499b92b797687a69aa809fc6c5b56be581de2f8668d38936c9921a16c921a18ae91bff15ab595897416ebbbde977244dbab4779d47bccfec14b1bdb255597bb9bb70e9372fc9afe475b2f73754daf575ef2dd565dfb4216208141fa99df428417d84fff2c54b1fba037a4237bb17b07ddac0f39209f83f8541"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00d11e471371b5ce0473a80367ce1b0baeb21d8f8ddfebf1116f3b3721247ec85f6e2786467b63743af0885e69c59d674d2b1a4b655ab15d8003be755fabd56f60ad3a7d2a5edbe942663b882e8c1d9aab7250a45b93feae3f092e8819d5cc2c0eee2cee0c6a098a40331aa12a0efc384e518036d382e4e231de3cf644e8aa8b97"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01fe65939e5a1c520be98522b1ebbd40e4c030adf0677c1878b1b0a58b72873eff6f58712e377457ef467bdbb4666e2f8a4733a13a065aa01e3f5f0cc0fbff0e8a2eb2d8d43b9f2a4931d107315943fa7e1d304f98838903897cd42ab948f7c5ce31a9323a35bdc0cae10eebccb5f318a1239f9b9609d45387805524d67e216477"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00c24248b5f6e52e0ce8c9068ba2b5839489d1a4849feb751b627e12d13722fd5a00cf4597e63c9bfd1a275b68489539f2b0bef36a09504d7539d0e1a346bc0dc5fa2c65c4c23b771a9946ef5bda403dcd27f496dc02233c05d7d7dc73f6438169a0bdc510bad2ca105d84c2c8bbf2a44c4d7d4d0ead980c13bda71a945d1f3f01"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00ab0ef4c1c3be6b7bb39ab0c8a1ffb2c12f8a2fb6c85ea1a8893f452dae161a8decbbc6a84ddc2068bf9df927c0f68a95fff1af8aa9eddd80b0c373b7ea750def2f6df54c0a7e50c16bded071b8d1df6687264e496316be5fcf5f9ab73f5c39b61a876441fb3f467205c92a864d97205032660d6eb2cee3ebfca9649295f6fc95"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01878ec4f236498bbf2320c89679639394b03dda157a9901f2e07486e64f1bb52f6b4823db13786296a71d6e65ad6a17308e46ddbb2608774eab3df41221eec799fc13ec95b567450abfbae8aa04f3c6361df3a1c01028b83560018b729b5924ee5f03f1306267eea55ab65a95591b105810a50111c9041d20b3ddd389e8ded20f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01602a125e7578a82e23051dd12ce12be44f2becfccbd13c2ee18ae1e391356786315832fe9fa6dd5488c83b4f560a5a4b9d9daae4faf0b9b21075fa1b470c7d984b2b43cfca22bc36ec305e52fb4b897445024f2ee536164a5a9a4201db4d9247d4e28e193ad3c62657a91b23727804e8f4bca40691eb41f17c68ab65bb8dd2a5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0109a25eac262687f28e560e73bd95be9894bf2a0006dc217e97547064d29be5fae521312fcbdd2949520961abd90b5a2ebcf55780f0d14ebda3c17825089183fee844a3ba0d132cf3db13ebb8f42905bf24374ac29a7b68f93f76dbce3942d4b1dbd91c611d24251b374bd29ae153cb9e23177115dc7003894269328d960cbbc9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01962b71c4824f2072f59c73cedfa26a49bd976bba7014005b6feecfc61c90caeeaa05ca8954219165f073bcdb73770846c97383ad1d47f0cf656830388fa5847ab9f542e26226d3e9c2a90bdc23819333bd13803f7520272e4cfb80b5c54c92dbc2936ac75f426babec5b49db6a64cd6eee14ecff0402506eabffc8bb11ec6c93"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "0102134c13210c561b22c8f2549e0a1786fb85900e3c69c20905cb46a3f633b7128656ba1644cb6bbfa1b5b0c5a5bed69a7802a543cefceb2132e0db7c596e51b88e62185f3815fdd40e7db9d1aed0b0f135b09c4d90e81fcd4ea7a8e7c150147bb2f0fab2d8a0128f25e1e498813f6dc26722a73a441d6e9ba4f488d96ee6d399"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "01ecee4b07f4311afce14524ab060a72a7198499342f099f681dd6b8a366bc9550a7ddd3288273ef59f62c5daa55c9c4726c78f08c20e0d9a74208db52f732377bbd8ca8f8f1d336bda6bb2defab66506c0db04bf0dd6f7179f52cfe9c5c91179de1c03eab017d7ff867478e45386955c7a5a744e7f8dacf738c80352a99226777"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "019fd1a5266cb6e8dfcff2b755624ec26413d25cf53a9d4341ff5c7b0b4e06e8246e6e1063e185b05d90f38637ca69c298d6a834e9aeb06e02afd001897c1fb097c905445b2e6d27750cef01f40d6030f0328eee55241137afead4f8d358d0be0655782a60265f0b9aa30b275a32b60bdb252c95d8d69b68e8a1e07c2374029bcd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "32fa78d5eb67eb14a53de388e9d03ae6ebeb7ae017dbae8f594b95f82f6ec380d5162f6f498d0cb61bb14d7ae54fa1b427c2a1d819133161576864a86d039200cb22c5d68716fd0e2b8f021cf25e08506d4ce285536bc6a074edb6d9b4a9dc01fd79eda19efd3b168eac045b6a4edc4c880de430dadc5dd3f32886b88d320505f5f0b064e46be0f1e31c57dd160e89738a4f6897975875564f20f82ecd4cc0db"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "5954649e58b4eea73bf1738957727ed4f356fd14891d95b81c7cd40a9ae4b9f1a807fc859d4d419e9a2178a369ae734cebf3b6b9b7069570515a94b5609585625a7aab4e2ff05566be39860b1c2e41910a07b46a555299a573c50b82572a8e40d70cd5949c0c5488582cc2ca544265e1e48ec5501fe611ee65de54946f4543ddd94f5d2c100fad681b6390924e3dbee62bf78133bb2ae6d1592fa5c4b0873635"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "282ca88061946bcd2fa15fecd98e61505b4c98079e5ffd08e9797059673150435ed47f6d94311c9df4ceadce2e13679b4eb1e7120f9f19d7ac393cc090d1885c88136ec24d085ace42e92ab049d8cdf963d8ba7b93b25e3c720367fa9d7d3905eb460c6922f53866fe439bb96f6d5213e66ede6239512bf0c2253ae23c3ff9915dbee4eaa576395e2d6986d40151cd8fe4c9b4d990ba17ec4bcdf6660459858d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "2161895b72aff5d2a865dac7e95cabaf7a28010da0dfb075f9b25c189821c99c1bb599d47d6a688254401511cfad26f1d93f254a3be2752a70f7859acad5e6f741848bfefe449072365616be7251781063e8f8934b59f1826341ebd0839dcf72b1735e21f35301313c683d28fb637f6f93453f575330f74e2a0d661ed5fe54816f8cd38b162d5e769c0bf94dfe83e25b6c05b7705a477ebf52ff4deb6bec6aad"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "71f7dff1a6a0fd66d5228398a7ff1707ed9f83b9b8c660ae57ee4dd40de7493cec1540e50b4586fdda98ee538e6264fb72f51682bb7bb5305285c287f4577023b8350a84fb088005e36121d9d137b16c4528b4a8a3934db88fd27128733b5f9ea78bbaf239c93bd9b6b4b1fb683e2e2ea911eb4da824b5650f186a7304031b62fc145a9a20a269079ba598dbd183f29a2f35a46eb05276b8ac99a8dc72d76151"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "55654725a248e323f3d4050b87acae89736b85dc8dd45a9c143b001685c72a70996f3ce99f40be4cdb83b7b420b520e7fa001eecd49cd43c31500c7c502e8c31e309026c07fcba386f0905da79d34b855861018af444fbd519736483fa79ab2d02182a9f0c0e514528f38cae7ef7668829b25d58b569027e4f286a71c1da3d9257a72a234ccde58d1604954d99115db265ae13c012125b5f317ab3297e5ca3e7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "6af6ed1adb0d772536d2e80f9f048b9a94cea70f6e15f37a6b5cac22794826089a11c8fb421b3bf8c108bd41a3cd7f34d09466aadc8b043a51b0b3e9c18e0c96e4c703343fcf68d45d5f023bf781de530a1d7946f4d2bcde9d7ef44374a2ba94ad56777aa113abb19b57d4802c18bedb58157dcd52eeca7a3837e65aa97d95f3b757e7eec27a5f890f41399aa5c2831f13a724d798aeabfb642a011c52a7c70d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "344b4e93ddadf36e039a4e97783a18c3a84f3d725d5f496f0b3632fd15b1a0c2ddf8f97a0f47401d0bef33c32ef36b2819f5d0f72046ab8bdd68fac28397d1906a1923f5ad96483048254e931a6acb5a3d31d4953212aa58c2f96e94dd5393f1e830e76264af68abfed551f3ff4e8d3bfbc6e6cb296befe2b9d694db4d4dd186cfcd6d697c7aadd92277f9ab85e000dfef3085cd52418d0f9b11605a64719003"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "2947f606c39ded9591b3314918b7fc0586888d42eb0a8d68bfa0890292f83f948280dc92e897c59de2477340c9fb288241737213d63d006a64b5d9c36b010164953fc68b3e4c7d70e4837b707a2b4b3608d878c7e5c122665299c012e2d5b3630b6862b87e4c680cedf13a6fbcc6eea8ce2d1fc394aa2327d6e0f41c4259b00fb8d8922b4a81432a30f7adf6477b5c436102c83bd1896718d8e795cbd5c30b65"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "2d586d8d3e1a38f532ed17011ff9d397084633faf6690129eac51e092c67217fb23e6d08f9cddbc38f7b3fafc308f23375df556f68f8dce22247da756e8aded669cb841b6be2fe5a22da4c0d06dcc6d6fd899d294ad0f62de03a7057e56ea6836ce8967d929f4144c9955460bb924fc32f5210919c79e9566e0552caaa130b6ab2e9be086fc97659bb2097adb0ddf82cca17b472ca511735499c448a8301f379"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "4ae1078c81d196eea211f9c4f762a350b4c060b4d3630bf7fb7dddd2739986b9de2422c9902e5870b3760be7b7926d6aaae633cf0ca9c0e78a2ee03fe193675524e0042073d3be737efe994b7bd93382bf8426f454e4a221fc899764f1059fa30b48ba6db9be33c92e312e449d190b3fa2f1c731277286fa363ac8420668239e0bfc26387ba329720bc4ed0217a772ab214a60d8d2d0889d887960383c420595"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "4acba34e2619592d5cfdbbe195d2aa9eed8762ac0a8336d947c846fc97d1d934c1ff42f1254de674990f76e514be53b2755cfb4ac52edec66a8127685c8e77e84b06bcfeda0684fcbfb20e2ee05c1202f3cb897bfb1c44bcb6301a9843f8e8eed031a1b4eb913bea04f13390ebd2a033ed151ef8b49b511da558e56cf1e3ac89545219ec026b3938ba9732792a1c89ca6d38c3c5e0e400af528ee477ffcf2ad9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "3d809b8c90e877efa20e031ec99d825afc1c1920d8b94e460848b80c3fa0a093ddff5c608963ab74f505a6da96b8068c2c2b3bc1676170dd0c2e65adcaf7cfd0c6b0309634961ad0c9b7f75e2f721f1f57fa9cf5d4f41f60b2ad3fc1d213b8e75fedb69ad157e24ad67f2ecc4099943e19ecfa7e1a34abb9f4bb02cf205906dc159c258973267731ce59d16552d372b9b47f0e630ec677711bc13995e00a41c9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "3de7d0bda6eae8145cc70591c4b78b1dd8d9ecc4a3d7edc1bbb75bf0e98fd3fb8d5cd4e94e4cd3ee246617b22426ceec6981681af9f7e6af08bc02bde7cbfa13301f7b88f607e1751285c4a861af2ac69f20d2d600e27b0de873b9ec7bf2cd0725b31032932f0f817084b347852613af9977931e2b3132a523dcd87f545805730b34db29c8c8dac9df8a50f5aa1e36a056ae41b01d04cd9574acaa98203d84a7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "5e2a15c7d9bee2668dfd689d027bcc37743259309457147ee7785bb3960dae3c8126655cff9e1302086adb3d1c962c3390f50ca3bf5f666e8a004930536c0bedeef4e8bc3f4dedafc3168692109a239a7d4fbd3aef9e6e0c8665c6379caa6ccb05a6f941782379fb13990f2bc104dc7e0007702c7eea3bb7ee42ffb5d570570b2f5409ebe76d7244b1e8392ccabbfda22515beb0bfad6c006c2a02a5e8526763"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "550fda19f97cdfbd13930911ef6e9e1cb2b7b5215a35c215d51ebffeb435642174cbe998f4451bde2d4bd2ce92ab5b9493b657f1d77d9ad4d348550247b903906109c608ecba7f88c239c76f0afc231e7f1ac1cee87b4c34448a16f7979ff4c18e65e05d5a86909615fe56587576962a2cb3ba467d9806445a0f039907601af77ba7d07578eff612364fbcac11d35e243734aa6d9a6cdcf912a2dd0a12ba7e87"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Worst case for Miller-Rabin test\n");

	if (!BN_hex2bn(&bn_tmp, "00907b5573c3d72ca5afda9df723d24066410e3d2b61f89c5c600f90732d0ad7db06a02e209f6792b609fee2ac6f3d73a5805f2b30642d1e2654f7ffd155153e5fbdcb17c76c27fbcc15010ccbfa7a1737cdf032edd5da7edebc9703e51572ce452c2319f1d91bee276d3e1121f9563b1700448ff37346b5a88098c9a682a59ccab86401aeeb74c8ce45dbf8b5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("A strong pseudoprime for the first 46 primes. F. Arnault, Rabin-Miller primality test: composite numbers which pass it, Math. comp. v.64, n.209, p 355-361.\n");

	if (!BN_hex2bn(&bn_tmp, "19bc037ff6b1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "01933ecb87a0c1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "021229a85a2f91"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "032d4a135c4d51"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "07277d9f8417a1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "194f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "0149c3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "1d7503"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "6c7e23"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "00f1f8bf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "0ebbb74637"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "127c6e3a4f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Richard G.E. Pinch, Absolute quadratic pseudorprimes http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.210.6783&rep=rep1&type=pdf\n");

	if (!BN_hex2bn(&bn_tmp, "15179c6582c2a8c42af5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "00800c6ed22988e8353348f28123408551ab4ee482b7961786ea4d90ed7d48bf4cc5bb0d7fbc0346e9ca2dc215540460df3c24bdec561ba766de6d618ce42fedb4fd84a67c5ef94323bfe88d9f55e1b111151edadda5a91cc0056b78c74770ae7f5a1af3741c92af4d87a70f66246fcaac1af0556b0a0bdd511822a01a4b897f0d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "44e282e671aa0c4f85ec68b2447bc29caba0ea0228b2fe7b08cd420955280bcf0ad99a0efbb8688b3b71a90a8f6e4b01911c689db474ff3685813fb2c943ce664f32d2dbc3c07387dec550207461270c323ef25c0992449e142ec3d7c36cb876492ee6a8593c4aa8e992c2f4cb394a88fa7aa9c98dd1c9e18bcf280332fa934b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "00b310aa4e16f59e55df118739db5ac21b65979ff5acd1cd4839716a63eb4ef966afe8a04a877548fa281a252c8a1cd4e62077f2ef5022e855d60d06a24a91cbd042323926aaec1f75fb4cdc4cbaff3a4275903c226d5982c22740e17d3e0bc7bf5bc23e7273b3bf86cad8498e79ffc43054292f38ee035fe9f67d6c542631f833"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "008126e1b6c59a80581221ccb272046804dc8bf7a2893ccbad9e61267f9c56ca5b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "008b741e1c47493e2ac2bd5f69f37c01ff0ec6a28e4ff91fea2ff24e2fad1b3369"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Safety in Numbers: On the Need for Robust Diffie-Hellman Parameter Validation Galbraith, Massimo, Paterson, https://eprint.iacr.org/2019/032.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "351591274f9af9fb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Sorensen and Webster,  Strong Pseudoprimes to Twelve Bases  https://arxiv.org/pdf/1509.00864.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "0331ff3562a8d7ff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Sorensen and Webster,  Strong Pseudoprimes to Twelve Bases  https://arxiv.org/pdf/1509.00864.pdf\n");

	if (!BN_hex2bn(&bn_tmp, "046fe40ff28041a690af557734e885052b879535574af06db2b787f926e85880060199697023504dd9c0d0e23b7e01e922538c586d676c61c972e1356ff053e78fdb481b7e5909c7dcf82155d713e915d8cb694a2f46320cb10868f03b98566022d225a97f1ee3cc26794b1e481abc61458146c48dd452ba81d06fab826c3ea58585500154d36c9076b0e1fd3d47222d2e8ae28fd5586818db16cc2fb9449a399ec9c22551448bde17c1e752506464424123af8de6b690f9407aaf52d8d279d11292fca1c32d0d9c3adb061f530fe10eca96e2bb2e4be1f6df1d7130aa21f78d31a312af5bdf56660247d6651168088ba0f1a7e4ec202f8efe5eade78726abf365c735736f578a57"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0b23c53824cc42b6875b787be423bd8c8aef90a1ccd18f041c8d6164b94e33a5c431217f4572779ef6475407474cb7ee0f49781dda2e903f92f5fe0deb0dabee93d47519b8c2633724e2d2f24062dc79c53add5dcf12a90f389ccd242b82323da265c6db54acbda0105dcce948c5450620166cd27815e22d3c1da9748d4b8640a4a0fc8ba0c11d0ae8965d436539e331bfcb712e4942af901f8e5c5a7d860b92afcb2ac7edd96d715d1d5ebd57232fd74c8bc2e18786aae081704a22efe24b4723b8d7227dc10d5c3e9be23bdd5c646d3f5ca53a3a725bf12009ceb98ed6e83f6ac611a0d582116f4d4caccaeaf150234a88b81b126ec1452dc747f46214d9c01b3005c2bac5fca9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "2085249c586a279f9474255a00d50a111cbe29b474218524fd3cb75b69e0737b9598905d046ff24235075e6df5a07a56e73cbbe0093e19386fe253de96569470a474a843a0211a17013e9132bb8a6f981a18d84b4472985936b72e453401b55c3fe3e7b28398964e2d87788edc03901f95411cb4ab849604caf42a924cf2eb11cc21336efdfce8ec322d27d2744eaff0dddf4ecdf6593485b14d7e7ff50b4d30f4679bbeb9cc0a26cfbaedc0c77c9dde1f54b21b3957c72f396bd7c7e2ed236a3b0dcd763ff85ec0190c7419496d4769a5329a9e8963ad3c9326e46a14b888a18c063e6afe7f350eff3ccea8c630fd4a024c908fa8248fe7cf1c3567f56ee45c1963f4b31225e6c3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "3dd0b362012faecad5221ed17f9dd0a0b1ea8fc23fa1ebaab3177201f76a8121bcd0310c0bf262bfca1b8f17a5eac72f6eac1102e7d68da9e8374e47dfed6619f39a1f51fee008288c72ebf3e0d7f4484d5d5b12a74510793c2200e51f8ec89e45a41b8986aad68ffddf864f912ea12fb889d937c237efb6dddb49ed6ef02e1d1612926c28a2c6f734350d3cfa600f2138dad662f835ecbf166795916c9347a43bac0dc95ebb8b75d9111a1e1efd8f7f6cc8ed276ad027a21090b41699a1b60f5239e7e7e51ccd9f85d10aea334a95fd09b5467c5f6da9bb10e12f22a577b99625be9c7b8046930cfc16ffae77c3733f528d0aee48421fb658d62deee4126d235759f00dfeab84d7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "03263db7214fd0fa5ccb86ed39c03270e0ba52252d80649979ded94b1cd23494533f7d63b47429baaabc9113bf56a785242753301e5a89bd0dc556a173ec596a5f4b93def5f9a1af18bcf228d37b8f615e0feade9b26d498946edad3bbb46183d2e69296a8d96ad6c1397f1e3a64d55c98fe2dc0ce73c3e15672f53e7203d4b658ef17239c4f45b06fc9e30913a8352962e73a47788abc4db223a097ca7f8eb6b404598ca135455758966e6975ef35f077dfb053007a3b63b42f17dc2f4c251aa07ad4f676b2f3c667ff5640470de7fd353e6e62377b0e272f9704f5d4833a9cd6affcd54b0639c594f5f7f1a666c26d6bde51a8590f40201602bb3828225407833a284e618faf89"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0f8451854b84d14366c21be5b7b331d89b1b83c989feaed6430c5e2a85acc3f2b1a09f3c202a99d5b92651d7a38a92059a9fe15ced0358fde59b492266f69dde4f8301d3e7808d3b9d023fcebffbad603908029251ed8a017effe2427527461d7e0d768bc3d726f540eea4cc1de1905301f435bb7ac49175d7bc7a5ed5a8139d5aa7b02d872c982db49b726ae82908ce331dd74c9c8d8056edf8a366e35bb22189d097124588fa9e84f6b8fc2b870851975e280f9b5cdf2f8b7c780454a2129ce315e74ff7e46961404304725303f07c148bbf8eb864ab8f89f6ed75ea2d5766250659f1e5a2c11492869ab3eb8d880f73bee69c7ce27702fedc1f672186df29d6c579fbb7368d6f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "02d248d11dfa04ee4f070722df4c1f03467dd32dff2d18d69ae32e6596567c28a3e21dde873bf6f3410b91a70b8a827bbbc1fa88f3d9c192210c1ae548086023d3ad5a340578af38271ee5bef9e0630b37eb56175cb1bc76cec3cb582bb88fdbe15d5190a5e5ebea44550cb0e2ec9e13098e210910ce2c6372d7a24497e80ebf872e492affade18fc4efc5c2cd34bfed582f06f0da6e969122f22057ce7a9a3474e41ad160db119e82f044319d4aa26419261a1bee786f6003bd6ac854583e7a5489ef1685040162cda798e079a2052fb910f2c36dd9780882738a526a31919420502614542514bf1c4b010ef32cf2e549b0551fb7e0b89cf48cad35ffa29310743d4224fe3ef5b1e5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "04b379213847bad82279fb3dc54d60692e9c128c2e0e5ae46d8388115ee6bf47a922c71e02f2f784e30bd81f56578fe16d901d4ac6060a62932e2dc8d349e1c029c98da5c558ac7da55f07e4422902420fe082018cad6f0d7e024318cb3b8248c87b7baa63d2eb1ecba32bd8051f53c285aad786a8eafc0c05b9d7e365495aa8f1a3afc1301d183be73b689b306c3e1851dfc7c91b88faa3e81b29e23c8c2ae86cfea506168b41eb3ab2a2e19eb4ccf6b1dc73055ce8eae17671110f365e7cf1db7f9a11d66ae816300765868b944d945bedbdd3a275e7faf6ce6b84f2de0a923c7bbec4c6e8f47522eb2fc1bad0f73a96345eb133b9436c505e8c2b8382e067c08f0bf33d1822a7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00ce30bb03f146cde0da64125f5d4df15d9b148b73caea0cd30bd06d6c46db3e86646994b6dbf12fe32eb708862c0e88000b2f44cbeb2244ea4920b15d82852b3b83ca6fd9676615b1e5cd2f4153854e48602684be12254b6eda528539c0eba1304bd37f329568636335db835082095ab4319374fb2aa0d61840ae25dae3d22d5f30a368f9130595c6edd667f0e6051bd0abf7512e973d2a7fc95abe4da8bdfb138740925d2ceaaeaf18fe2244e656d3edf46f6c1c40d7dd44eb116d321a33a48d0641294eeda8759ff5bafd3301b7b916a089b82a725b15dc6634db88dbc092d9dbed575676126f0a60273f24759b24762926a95669148ae8138dee6d84d242a5e9f2b1cb6dfa1633"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "01853dddf18e11020af425c8fb280fb606868aacf59fd8365db931779f858d60fe61fa2441591e24aa4e409dbfce513833619710c68e1da623b9a6e5c594f8cb8fcdab698793529d70c4f0079e1ded6e16aa1b42cd820bd72eb719185c61596db069989b88a8cb496f05e6c8b1917db58f145a679468b6406e15b76b25155402acb4742702e8a5d212e3fbae3d4ff06b91ce6de68e9fda7c5ccf9c591aa0035529fb1c8212a35d74ba5e66cf60ab62c47e7d3a53babac9d4406f3ebab673d2688868b301b7da61e3ab9d8ed91b874a68a3678db9481ee2efb17731c382d232a6303b901054a7b22edc92e31c497034c824b6f065a008670079e0c4564684c986f141d71d0a288a038f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "02556ed72094c997f884b315d0355be63eacb02a918a75907cd65d8b9f105ecc13412a8b4f7a163fd00f62ba434f42d90ff5b30367e9cc21122859ad48b498cf72fb0357672ba88e362a52b211b5b317bb6300f745063bc3685a7f4ffaff32018ecc80f44cec94faa3f35babb50de479433a084662009e70ee4258dd6971aa0973002bd507b4a20e8befde99149b4b9036191149399329e39629b0ccccb5b1760c5ab6f50c32a3b2c1d5f85ca2d33a926ee7c7b35dc363d44d5062edbea7051c4aa38064c196394be4b1b16da35131b02c04bbfec11da64538f3922a582f423071893c129def2be77c738cb37d4ae35623379f6daf129fb44625616ddd886ba1a78c12258f9af7bd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "078af4b6e56f391741bfbc899f3fefd5e20748e7816657e70318f16445f27366f67b931062a8716e3545024edc4d6bdf151f59770772f45fbee812a3056ef42583f37b6f81add2e0522dc11d23f06814f18b379d139cd3773d3c0bf5aef4c82f1dbf69d34180a7720a029f6b283b46cf045c115aae9e5a403b830000c42d592ccc42fb2c6233466e86efb440716fbae0e696114b26f73f8c42f90dae82171ddf96e0755da67c788ef523ca0cce19b432200af05b7314639ac75d26b77d86e08681917ce499f71e8624607217287d0b45898cb69f1323f43abbfbbb758ec3afadf998d27bf30518c613e796bd5f1b7170dac0decd5ac7ea8bc552dc40e2106ce5f793e32bec01a209"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0b9c27e0c46f3de793c85f2840198e51d3c9550751a2dbd855b364bb4da35fa13871bcff3a049631bf9586c5b261626be4e9ed8865b4d9dd435b8b4731c5c9ee7fdf298e2bd6f7a661e360cbe764a7b7a3723fc8df5181b406bfb787dbc1c2e1586f88261af47c8997e71e79a5ebda4e01d5b862f4823e919c4b07a3e1a94acf139aac80d490b8af449d88a9ad1344afb05323d7400a53d17d28e8495ce7b17d182872eef67479f99cf2e8b9abc967618365a4154b4004184db43cfe2476de3f15301708f576712e8bdbf723857eaec4eeccabc8763e5ba2435c184c155909d4ceeb7e34a8fc0acbca6decf8bdd360c63ac4f5bbc307ff2a7ea9901ff48c12cde5b7544ffe9ab55209"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "3c30f69630518ab86d506ccb13f843e64e257e135e68aba14def5c7ca87fb23f606d9a21b43825d46d3249372f6a6734741d9e2a8761c43f151defb35f22a58223a4ea1b512da6741523247dda566b8ebcc070691541e58293b39b3ac06d4055a652d7e599e443ce5c59067700caf6c5c0a9f75af9a1ea7ff95720485fbdd3eb9e3bc28bf26a7bd1f8afa77c99669254e5b88b056af64002bdcd6e1ae8186033c2ba2a92a2894d6a4c1ff15bbf70b8a5773750b8b96976ac93ef39f50b9cd3c54f81c65953629afe6cb0944249d0ab99bed92e57b79244948d03681762661c308ebfd0cb89d6e7925ad2c687b8f30b4536766fc28bfb8486e3791055604b3ee95085cbbd0b328f7f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "41307dac271321e285ffc17e39f2288c668bccb8c967bcbbf86cf833709c6245ee0d7d9c16a65fc414b94ced362790411f378e2221b8522c57da51379f50323f2554ca3ea1b79145fea625ccc2df919cb31a06ebbcd636e00e293da1dd5a6b288fc3d5c5e68491beaa8be6fc815c64dcec6e12963c3458fb57090d6c2c2c26b77606c593d711ec498727cf7fcf362e46f86f24ce85df786ffd302e0d927955e691c5bf2a0ebd9eb8c2742fa8648f82b3ec179b1531749f05cfe67f3559f371bde2627542a7b17262d48fe630fa7c59495cc7edce5489319df977405fd2042ab0a56a62d478115013eab1eac6b37f6de1ae7591d4cd15fa344b05bdfd996c6a200bd2f588daa779"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "131f432c0f4e8f16b2e068bd41500f4ce67285268d21a1888ed5c13225a2890e17be77444f9c9ebafe284c2d36c6f66fe4e3ea5f64092eec66dec3c6d1b80316517fb0908cb67d6d4d783dc98b113f456fd6fa71f066e2e9ef2d5b665600901e6c4f304b2c230eee34c3516bbd547c45d4af2f41dbb6fcb6fe60c76285bdabb82ae6cbac84119d8783a7341fac7872629830a20c17cf5131d2d5d0474a42ac4972d1ba0cc5a18c0af70b6ec820b7d2dc34b94281800112ef1b676cb06ff6be14cd023c3c8e366d04d14118d7299d3aa10986dd1c2df41f19df9cc44fd7c2abf22b59693303555b33210c4ff4d120a4b8f8559e3feaab4ff80b0511f296db95f67ad6a4b0e886f7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0735440b7d3e3aa4be783bb912f644624fb7da694d092d0d3df1d5892dc3d40f96c2e5b5fd3b537a8f6c12b1e1b5931ea92a7ce957f5d08682c6828e3f864e29dac8c2f3f6a4cc1d3c58eb5513c4bccb9ef9da3ff6db38547563d34f94299c73baf7db8bab5a9ff94edfa55d100bc1c1b1a17f75afa619c577019304887914b70fc72c25c7155085dee797fb824b5cc1d4794c26810662d471acdc625949566d06b734408f47a22ee2f9d3566a200df16ed0815ba6965a1ce49b91708c9c53f61db16d102a6fd3d8e1de82425b50d0ba726aadc4013ec0aa8fb0d0a86ae9b025c56d99c9351c58987e89865cf029e4758aee4b03d2e4962ab1e702a46a95986ea380ae3a5e4d3d8f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "027b903357d9edc25b5d45218e0ac3efa2851ea54d84ec92269fac38e533e41ee68a36f86e96d2582d3bcd16afcb7fdedb0a58fb7ff8c94397dd1e1abdec786a4f94fd3acfe15a50f045c2b7bf614612afc4683e0d39f5b100237f52434dbb44eb264da762557cdac6f4aa651f0fea7a9ca7a04952d6f9b0031f2c2f318325b4b84435433578478cfc215506e9a524a8dfd9c7cbd71c81151bc25681261da8fac3220ab32c5c4cc4d94d0febf6353396c6324fc5ed2fffbe6155a63dc74ec3a67f4a38c6f138d91876783d1b9390743eb1503887b041a1f47d1ff564506543ffef691fa56794ffc4258ac0e7aef7e5ecd0749800c68c8835fc8a3e7118166050bde3e9a4e110df5929"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "74e23abe7fcd90a7c0744204a47059f8fe6f4a9d9fdcd011539c97919129f6f46f310ad14f3866d7e82be737cebbcd72b4f1995941e1ab96db7c563444158bef8b60de6b98aa76549fa9eeeee8018485bc55f6f9bb8621321072283d9736acebf0c189453033879fd38f141a316a80f6c2d5d2df7c03165ffae733ac0f060d9d5969446dcb5ab8cad9853486707c1b373f4144a61d1a17a23b3f1171fd06359b98a3b26e4d8f4cb7f83e91bdf9d7a271aec906f596ab47a001c07e78758f7c0ba25857260e3f91cd21462594138e6bf84cf1c0cf60a8ece8cd2e53e6ca73305428af507326babdf37e29483bbcc2b6ee7b058c7d9fe0b407ed9b491ee85e001a4dd9175a5047065b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "010ecf0a6d6fdc6b91c46ed7fba35496708c494b8772bc007bba48958a55e584a40c9a34598c31ace09afc982944c860f6794d5f91f5d07ac3f79758aba4739b592146dfc4aad9bad99aeabe97960b7245b3e62e04f49cea54b33ab2caaefd84fcc202902da5e35ea446c0057d6015833f4e63d793cfd6192cea8736c0ca4a6c4a7a9cf0d3c8a5820384ff1728ea09900c0b2c3eba1fe588719e7d1ddd750508b28b4c5fed49a03b250a424260ac27ad46df6b08554c09b75f80505c1f31021fffc5118e40f523fe2ea437025acf3a8e6a23ea6a2863b460ffab45e47a00c5f8a01427e3986cb7520b549db4aafdc9277fd122787808b519d4d7cad2d225b5c85f253e9b5aa19d7625"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "01a097c0a285f02b54ba79a7ddb126b709f41b47bb9d8913fcffab4db0f3ba01766b502e9b3350bbfec70addd52ace387622abc21095eb7a019acf6873b9d2074edcca2c5eff5998a3e984dbee023a71d62c0bf9c771d84eb16dee06fbeb7babed577e77cab8785951af65086fdbcbb15f2e1c018192a20d6add44db22cbb14edf2a140cc04f4dadf8284fa77fdc780d9ea34eccba9480288b6b776f09f7e7f4b9ea702359c5fd3cbeff5469530413c891d8df081f9a25d65173b14b313a8c3b75f97b56f053b879f7e31b6d5cd093a47227b2a16afa4af36fd2e91ff1827be9b5d59f537082e535d59788eedcf07c7a61431bc30c55bc9cb93db60fbc7c747badd057908fc0fa0be9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "05dd50ba142f4e7f831ff9ac6f4dc4a4582d245cea39319308c2dac8ff1636314b0cda157d0d5ce8dca3dfab608922d9a7071b8478e4e5311f7283d30469be556a5924bcc85338f1f06a4da7e13fa7dfde9fa6db76ddc5e558619d5933a2022d633f1a9ddbb2047c8ae585723f04a69e8c2e01f09e9ee53d3bab6d7902893d9ad725e08b0ed4a25b778addcd20a9439da8bede1a96cb7fa1efa149d047bb08771b59e22763ebd098ae394ece2912d5b2af85f2499b44bd4ea2878021a33a3f305ddf9e1860bf670fbf72d1f09ccbf87b00cf996a719d5b5c2728ed3963e13682784f00ca7b6f96eae3879923108e6432fefd20481f72dc6c6d10ca95db052e0c54e294283e68248ce5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "01dab8d833a71ed5abfea060c1a9f2ae09ee0931516fc1f38c14b959cddd92c3aef0574b9f9c9f2991c38fe43a7536c081e14e2b87b7b0d495e834650ea9e466783d4cf3068382cf8ad00651669959feafbf336f4be62bb4dfc891794e097ea53cd8f800d79818127258f89a7e7d6b3d05e1e3f0da7ca6d1e343d24f82ffa9d96fe2db279d2809f0b6482262d53d32677f57500aa703e5ba9df500367ca255d051d7ff7018fe687c907a2520a2b992cd4f7a10b70b1ba3f1b2e5ab07de0e06de76affc27e6b29aa2730198454a8fe529963c27260729c8fafc6e14594604d8e1046ecb7f88d8ed100280f42feb39058f17a5c239848d08b85afa976efadd0711c3253410eb8d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "08f969a48c87087c160937ab35e3a80a04c58eb3620eab93184c7e1e2cff1d958e92faa1c3ff3bc17360c15f221aeeb6af889a95df029dc5c4f9974fb77c86601a4b13f872e57482dfed06c4d0055bd478408c40472d599bdce63c79f91240d448560a554673841ac071518c627fb0f22ea0c56b88a1ba5fdf427d5dc68e8d25d944e0ddb61827aaa1e224f0812acfb1158e37805d84e0957c6895b07913141db56d4b41996e3043977259ab2aae56409146421da6f89efbcec0c2cd6c173949cece2e402139e9d5c8cc1a0a1832926985811dd052cad509454c51ce4c2ef5b08cb04d6c497431adc86d43a27bd4c7647a12208ca663f5fce246f4045fcfb9ee8ae5d48a4838a9f797"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "3dcf125a77c4a3797dd742d5f4647f22107fcc23597d261d42126ecfd63c979f0069c0a44c9d3bc3eb47f94e83041cf43c90c9d685f61d1784af6560826a7858807fdcc4a62e8fa2e24c060bed22bf5164ea8b193248698a59df6a8ebcf2831a746be18e5a7fff4ecf202ef6f872a773463acc99233dae6731e20db3b6a1b7d71171754866a9ec6c20fe99a06589b2f940a076068d3cc2e5e199a48804b6361548620877d3b65f2b652ab5029b7e964b465bafd5725add9461c399db82688b0f2ab510384fe387e8f289c7982d3952bdb61944c37fa1474a67a07008f3cc7115ac907ed22448808842c247d554a3f3e36e6665ba30d489723a08a8342e59dd2f5942a54f0302f3cfd5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00e10f450df7bd6e44aaa5f66994d5e11a57da6947969b8ffa84fc942725d0bdd57fea4cea3907cb5c8200f432453ca855e77b4e89766b1f339f1e79f1b3bd5a477d7af9c21d97ea2f025ec6810101b103d496141715d61764193ba8bd63b00162161b213da888df612e610c8b3ce100b57ac59d0eaae65f6ef136d5c42c82104dd37b483d68345216689fca1122fe3e2957d357df3e1bc0a7a23b3f3789103fc8c8bcd6a6a966e2661652e892c059612770425b251b8bd02a0955fb5d895ca90a447e560d13b5d065a241777320c3dd839212a9be7ec0dcf792e5d0383ddcc98cb3cdfb85b05d3cfe3c6c8117c76411d76e5de85b1b117b22521d01728da606bb28491e2dc93b917f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "034b9c645b94c535846280e92897209efe58970459e1557d61a3a0178a8e6f2e522bb4629291cd32b6357ab7b0fc121a7c62fb7e3de939391847383b1be7d27ee8394d561d11532159cd3e3ba5e93d49466d1ab5f0196fcb3ca72c4fab0fba4abe918cf22972af7c34168e49a5ffebfb893dd0badba1355ab22daf54422271333b2565d31298f87eb0c9ddb32afa15155c611249f3500045e17aa830dfceee724215a633559f9e65d9603b3b8a848025fd6ec8eed39f9e4d095b08221edad29372c97df63d151f68c3b5b502a12423bb961e51a9626ae8ce0f08b7ac969e1d0ef1e5a04fea3302868c28e02e85eb79ef16c1c7e45d6f68cc32c292205b74ab40bc02cbd5990fc7b92b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00c1221e5f3877bd767a7c56286fcf77a3a3c96e1c81e15a59933d90cef676f95fa6cde09404c8988f5094edbf5589a01abe9d612b858068ec2c1620b0b8d49cb3e431982ef99104dbfa95f8008bf5915cc42354a1ee2d8888bfa0d2b964d9f664503be6a1c6a99a121853651a063c33bc96ba1021bac44151fcf92c8fda6107bcdb4ab61bab8588e94ff38adc65da325b42b1525c635cf096da2da789bb9d97edf07a1d292d9b8dd7169f6292182dd89e2d9cd7169e20b6cce19f951c08d48b3466e134664a6a45ab508e502e3a17271bbd44293b871ae3a61c5168608545be5ffc889ec8f2357b21a628c9af10c1edb37c8442f8f7676663ad9fae6ba6567115a89f90d1cf05e999"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "1306e7699113ad21f3d91f0b6444a2c65b3261a2a5ff51ef0362866d81305e8f13bc3e112f62005b7853974f9021b4a30b559e699282673c893a7bb91c07969b572b98c460b483ccae8acf42f713da00eb6d65c7123212cfdb538e98b5865f5d9b20ad1f7f9b64887f9efbcb598d7c864a6812bb2f7d9d2cf8ec3d3bbc7004d5316556dc8b663bcd285741ea061dd735b31316160d869b097e44c3042546befbe43e63f24bfc870dbe0f7a20b887c384eae7eea4cac974d3ba610ca6392b75a6fa4646b111a43a6a729835edee935f7019f3cb0929c8858b390d4097d9f6b4cf4665f925bbc8e85da11b996980556b3e230eb6d59ed8ba337018745f16d7c6f7310db87a615257"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "4eb3ea092cab164f3fd7a73136de87896de4479e92ba918fe1a29797902db20c2ded396a4351b61aaa66a0d142bd2d4f4b44d39ddad927fd38c1c8e773f993f9d49d6aed9af93191408711a0774de82da243279435594c49950929f074b4b95f2e50f7d57a9c523bcf30b8c627dd142529e9679bf4a4ccdf76b2d0077b40a6006ba8721703378b8538064afebfb97c1fa8c49bc704b99675db97de4eb52e9cb78a907909221d492165f074421033428baac5c23c508c959d43276ba840d8be98baa38f89dd30f2c67d27dcf60e69af725541538cbbffc2ea804a34f861fd06ed03c682c0bce11cc0c16ad164d846c478a55787f162d2943b577b2cf4483eae13b2fce80f436a15e3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "020f1c50b3640950fbf9008742b8d993dbd9657026d755691cfb088f0c6f9c4b98bc0be2e0c8e47881e7c9d6ce6f35c08fed549bb40b0f4fd0f79432df4b5a5f78d2b54df8cd958da785d9ab1c727c1efb0e667bfe216e7d2955dd490e868e783f1409d0ebf2e079f1303f57b50ecc3987a53afea5d824f8dc9a89438fd32f1f4b3a729c5482a3f66cd69e712b1fdd3ed25836dc8157079053bed47f5e500ba698ffa7b6d02100f70993e43bda086dac726e72f9eacc01a1d623edceded81e0a446c0713b06d9224488df1f4239a7d99daf16d5273e0bbbb11360dfef18ce33613441ab6947a0daa61ecf0cf732c4b8141e951b232934b61073455e454e131e442edf39cf62e30b101"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00811bb71ab010d948e4ab07149e752d4c0b9bb6aa11842a1146011c0da66cd8d597f7dbc48f26b9445a92374684098b2c87db94276481d79958a8425263ed6ebfd97d4e42a5239b36e079cc14fe923d3312f62800b153c0bb4e4e97396182b6f1ca5eb6f33ec61d4e7c2d822964b679ca01314712c931a8f430011644ab9d47ff485ba18041a564464c806c0b445a69d4fd4469939dfa304d8aa11fc2e9c98b450441b5658be9ce498f638aaa6076baee06c31f66751b440f977543ae6b268da016aadee31ae4866ecbc9f57d077a0cd23f802d27875b524898ad2dcd19e91334b88ff23a7532323984b040c3d50e6b37044b89d6471f92d03ddb3862530b8a95ec1e10e40b768ebb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "1a110258b0f5aadf8b223c58bae17d256aa1b66ab381dcb1ec128c4cf8d1d3bde3dbdafb45808865c919babcb5822f5121d6efd80e32496a66acb9c642fa93b9dce7181295085009f2427e1e0dc6bf322b8f6b45219b37640119bf01f468a16def4fdee8ae8bd10829481a918069de36d161dd5a00f426fab3267ad043c601a2109e4e40568e76bf97b8c64dbf55b442ad484ac3faba1d654c1e27ff6cc5a215ea6a695b55ad0cd71a14e3288b1c2221c387b8667e8a37eeafb5703b4f64b13444330cd9a292395f64e26ff8e27ffdb041a7b3d559b187a39df9be773916a4b7ef968892ea923fe79291138a8de437e4617b9e43e3f4d0a0c3933b6a0babb54cf69756e6457025bf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0315ac7638fe2dde5f466b9ae990a6b4ca6e4529c86812f1148e65c2268f24569aa0d2fb1a9b4ef4059cb4c93b2c537a63cadbc5de9ed118d4ffe2fe845a57b9fdb1dfaf19e50fc09469beafa470a45baadf99d46d86b23b0f9a7e2211fd5119db7fa220a819bf270ed8cc37df6cdd39413f566158375a8c6ba19d33e59b517f23bfbb7ab72a4253d3b2f25450cbd4dd2795bedbd6267d4bab9c58cf7accf9090e44e932886546d30865fa3675dc31d88c16e223553f4c50e4407ef44c1937b2da3447bc9a9db838e8cb709194b84d155d7dbde917c485a6b95a884dc1776e96c51641445c015bd709b6d1f0b64349092dc3675b51d15b86b6d73de9e08d61bf3da3e7d3be9d1be689"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0090fc577b63378b614de9a87438496201917d1d98ee56b34d1220aea2608c9296cc10686b0ae9d554447ed47c5bec9b489f0d4456eb42cd7559bde32a3556a2e7b6b61c868ee49d85b8aec0d5d17993c7165ca2c0accf59499c743cc4c6a50836b0363284b0c7552d8435f2a25257bb6f82d484b1233ddcbb7c8a3f54027d0bbbf935f067dc3579973d1d819d90f4311fd9ff2ae23b3e8d5e049da85d70281cda755de9c57ab09eba0961ff025f3bce5bd3974883836d8d3b9d37af73cde87700a46cb49f424c2264cffeeb0941cd7ffeba9202b6f789d2749860e46b27209e9eba449cf1794944470ece94b47c092572616ad2f4aa3adf17099dd1dcf434aaadd457f62c18c05d51"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "016ee5c4a1146a72841a52e1add93bfe32ae14ad0ff82b2879ae691f5347cb3866daff3b94cf40fcaf2efbf8be197a1aad8408493becee6f4fcebe53a43eeab0c444b6ff50cc9661d34b3671effb555ba5f7425d3c99520b29c5fcb937de5c45f0a80f7089fcf6a5e212cae6b68c6811ec22e71706d86dcc2636cca099bf8066336b9da793f86b4780c838145a5f4d079257fd383116cd00b878dd617a984e3694f4ec7d134653946b81b12308457dd4027116526964099f52f2220778cd954515a705080994d4bfa1327168121ed942f69712d1d8a21cee6a510d38421472179d085908e9993749c2973774b9020cfde097dddc28d7694d65d28a04684640ff90d5a86a0c037f2bcb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "09b96de09b3c269edab8cd497efcc8cb84cc2693656dd8b8454423efd8ac844ddad1ffbd5de39f4fbe47db100ee56131691b80c974019ebf319900068a646a6ac837d69a0d3470f4fdb309481bf1b1df29aa70b1c793094c0a78645216279a4f592ccaf49a39740ec82f2656fc8e343fe58eb4f205afb197d488843fa3054f9023064cd534823b87f69f808c24690ba57f2307c47c6261d1f240aa35c2c47bb0d89b18f071e7f96359fd91f8a5adf68bf86d49aff7030c5a106a39ba388b471bad93b49c69aee8d8e2aa12c6ab8ef318507b24603665ced96f8451c5cef5a3340bb4bac1f577cf0be337c1ee8764ea2b48348089ce0d070a0d7e1a5bd735f636baa88f9d282efdfb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "25ca8e8844298f700c87fdb4156abdfbf2540b4eb16ffdec9d6742a43514e48346040b4ecab2aadb7ac43b59fd113ae2c5636459c964306150e880e2c688272ab74a9e0fbeeffc29c60df8d8d7e696396ec21e80c2529e12bda83a1e8dcb9858e568afe89a79fdb00f766e5979a0c7b48168ef845ae674ca5bafed340cc93d51ca130e72dad8497b2ad8e321e498e169898e6c1491a12f05dbffc31a81c859c27657b510a37914676fdf828c43d4f308e6ec42de80c44cd49b835f6efddeb89df5fe10026c3eb0c6f580bf1a2322468b56ea60e9adef61f06b211b8c072f9a52593ab333dcde7c4109d6c628e44b20fc0e19476a72956f53fb0c03cafa56d0a3ec0e07fec558fc31"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "4a6e7bc26df7ed0d7c09e82af5b1905f9836705ac7cb854f4be1316dd2ab97505cbf70d090c1feb774e3e11bb4c0bbfa92074b2a59e49d2d6afe3e0de31824d407735b9b7e3b5ac9dc2bfaee3548d8d3ce02e83a275af3933803e301e23d4244a543fd80ff79e1fe751f9540ae7ddd23da5930f01e041a095bd5b505ca33868000588a2000938245e75135744dd8a4da04a0288e78fd73ac0160f3cb108f212576418482a581bcc71902f598d9844676dde99fff86be9c10e85036a60925703b80831dbf6bdd75c61b24bfe1ea22b48d5502e5a52036f59ce0332c71836623c22e2dcb9f2958cd4067041d4c4596ff98a88ef53cbb82f011f4346debe204f5389863a0637379888b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "049a4a9991305451a4c682030ebbffff6a0101f04b9ce286965fa4afe83376fd028134a9e45b3d024bb331e6c80365398649f591ca0f32517171ec860bb9d9f7b415ce4f4a702aa3cee416a82b51182ce23088beb94d5afcb7d1b0c8b2a2e47e7ff63afaac28aaffe7b2459628d1979a1ccacd028909db31641a40f3a3f742fc993aa36de8543c19cc05fa3bc6031db33c56a5810c279a0f872ead931e85c5b55f71ba7232f6f0d50e2c7a614f9cd87938b6df53df6a68e492a0715aee49c235b954aa2fb6ac13c9d64daeafa16ff4addf7605400538ffe6cfb17bd8d694b3a28eb18c2dded3be5167b357a124bc38376c74c970f394e4acae0b0bdebb5d4479073ebed1829abe37a3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "03a48ea2550158e6910298c3a4b6162e9849bb91378d93672c95aaa20c8470ae964d4a11d3edb400dc032f3eabd44d0016255c57379e2765db31bc00b83a7914b048a28799aef1a74a35abba31755fbcef113c96deb380c86b404e961c28a3fc4bd1beb71f788e98141b1b7ba70365e3063ced78b814e543405ae80f6135c9f4a9c129bdb8f29a25889a07767339a1de2d5720f491a8394651d6d34fafbd6a63724028809acc69b9c542f107b2368a74db0cab8f00b4f7006dc619ad1a0b2d10c38cd7d05407b117a6bebd54cefbb552af1b0b81ff21c7bf542140f43cc2e10f270180bfb7b1665f09d36ca5cb86aab4ba9015c9fb6d47b954decdddedb1b81c7faa84671bbd71e9b5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "05725333dee4707b0e174b4d985516515e17ca12661d4714f75f33757ef58771d6979d6eabca6cedccbcc4c793afc4d7d6a429d17b7c6bb5f1472227e67f69d6a32cae34b8eb51474308a3ec274dc5c9c8e4ab1fc330957a6d0a8c10abb4565977fd780f74905f7597d01485010360ebfdf7d22ce9fdd09b6234e4a3c9e7f27aea55622672f89e4d545d9b7cd73bb0d4312ee9bb614592864df9c3eb88d50f2b445c64df29e7c49b7f394c5d5f20cadaa88cc9f647e89cf1c66827ee51f47c4d2ac989ea5a87062cd8f0f083de4aa30de074d5f2efa3ccc4a931fc861e8fa64fa4d1db13d86046d608cd14eed45ac217ac67e9b70f566688d9c1f6e74aeb54bfa008d5b12206025d6f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "121df27a2bede84a4d79f74687ced480d4e12a1330a2e078aff07474fa91ac35a47ec70ba79249ceee6c55e6ba1fe8262ba677e794a51bb03dd6ffc81577b5a40332bb43d7e621f8eaa91ea32fd481275583e20aa046cc5fc6f0890d5bb68f59672d1adae312b2f03c070b36bcd1645569e421497c7cdd0f7dcab105d3b349ab0c6014d95dc666e35ef99854b2a7d75c533cc77b9925d92fc9278effa8d7f3b667796dd689499b90b324f6d8b770d250b4275ef62541b58efcda44834c934caa0d9e648a72ad7c61d20ba1b457ea953c968bf4eca2dc1cf46f0e33663b07479ce3849ab14c1240c177cda2234fab6971694465debd591a1fa8cddcbf8b2d52cc711fa9caf58fa39a65"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "013af9b522dc6f0cf460e489cdb857f58cfa82ff572ce560f8ca917a25a3f840db925c836eeeb9c9dd1cf790624902eb513e957073e315b050493f65a4ae753f122f49619bbe3c13d458d4875aec14b0adfcecdd0e8928a2e76d2658788b21110e48a87d877f7fb1de6530adf1613dd3a719bd89bd5b9811c052bdd481510853c3cd9dcbd2237174b601ea589c6bfbf8113067c5f17b8c6ffc2f761ea06184319fb92048fcae4618093716c6f746bdba565b6cc5e234f7e4b6ddf82f4937f6adb12aff8eb6d0ee83c2482d488669fe63938a41426ec09165168fd4bd284294f6e8975ab4a523801aa3ace206b5abfacc1fe4a33e30a2fb1b2a0c36db02f702158d156fec37a73e06dd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00d5761fb975294618b26559075aed6606c562d0c39fb1605a20522cef749f42dc0a2a6bdf4ac277b5033b6e3a3392fc52abd1b705df5e0b9099df20b8603cea9b76519153f43508248cd474ce8d85b657c440ee715149ad70eecf6f4a362730d006d09ec09f5c5eadf3738d8e254f208c80f1f2610c94381031c20ea82caf5d5ea2ea35ec51e4f98d352809058a41f6433d7efec539e461695ecc39131443e0bae23ca985bd0fa133b8945124da374fa465ca3b18fd260197a21ab19c38a7964c47b42bf3afb6ca7acf5af0f2741fbc02d3b894b8a09168139f5024c74c709648935c06f91918fee75987979b8e045987451dd887d502db27aaa8171f50442b6b014eb219495babb3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "768acf7297a99d332c6b2e86ecd5d546a82aaa2236996da7ecd23f52e4b49e350714d5c213193bd9f29dc215cc513972cff3f6d5ca930588f81f306392ab12a155ab7b567fef2ea299c6e1a20293894f6282799364817fd6d74f5a2d8d0e1a096846a87f9976bdf4095b64470df394da237b1d6c5cbd0840959ee1bc563dac61a3abe087e5786332bc05456ed1a3034faad3ad4488ef90576115c5422ca993e154843856c7c49dd82c7f24eddc0841e9a483ad2e3cd5d9ad52d465807c0d84f61f2c70bfd372b2d57e6a2dc973e9345cad0b5a4a3965564f21054153cbc037ab1a9236de907c26d0d959a08427d48663bbb35a5b3a62071951f05e139b0dc9ab6c9600c62e91b6dbdb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "135de6a816fe05856578378de617bc97e7b98164510c2e3322605861048f937f20b92544cf34bcd6ba51047129bdbd476e1d8f94485793e9b0a4114c544b9d56f0fd8cffac0cb8041a19bc4c9b6a19f0918f5872db68d393fdcb04ba921dd977ed4e867859efe3e1d7ea6dacf381a4177363d35a011454bc07c3e619f8eb6c6d507b2e261270bbe379b4d83e3ab6066d567505fd1d2c68d36fb8379ae88f145bf2e732906660ea5f6fe6936980b9ed54c4b4330f0910e026bc637219a7f43a9683c7433592c4ddb94545e0957f8355a5c3a4e2819a8c40e56456eccd1c5aacf44d35eb55c32133e970fdeeaabd4c5da00ffce060caf40fdc6e171e5f7522364a06411b7756c5d4b799"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "18140470bce8e40cc53e33926c4a2850b61a15fa099da53bfed9a0393b50844c42520b0d259378f9b8246eaebab0f18c19c24171e051612e8708298820773c86235b0268544b04f737877fb261f40e867d9cc54b5048063f643d1117da90565a366b9e99f754d75b12ec8a3001f9adf20021e7a128e9ec9e8d83871f63f8090c7874f58a4fbb87bd4fcda7ad61c8a722527bb7ebcf4685ea61782d2aac6421ef1259afc80820f216f2ef143a9c3d6f79958b6eedf7b84fae64e497c3072558855e38274117b676a6a85f444abf3ddb66397d4381f34e048a6d9646c188fee062d67c314164e11a30acc9ba5dd2c189beedac5052f02a231ba42ba32e4336f7d3f64c2ca8117d19e491"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0f32a306229ec4471f1bbf598c61ca486db9378e07563992c45cdbd0bb0ddec9638ea89dcd7bdd153140df285c600a2716b5d2411af9a5e4f935dc3d9b85f8a85d5c24685d81ce9a5e3bdd3a498c97bb2d2842674dd3401b299c753fe8cb88c25fb0883d47976594bb22c1e34abf074b2417d3abcc787ad77a95b5793abc5182cb3c1eff8179409c19ab5c4162a84b9a68377bd1e5258d6c97afdbda4baf3818d949862df6a7005b343ead167b5f324b137a5697969fbbbac6f26c7f54288d93965255f7d82e1b42ab078d5234346d21bcfeb7aa07078624f14a5ff86fa5063a81d730d0598a0fc4f2b3ab05a6526cdad2c5dd2f52eb1c0055085b5f1844f0554e4af09ecb1321a625"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "31bdfcaaf28d96bfd3c847ea06a65ca47f051b0b37e0ec4f06b30a4d0096a08da9861ad5d48b45cd2a79998df5732295115ebf69a7bd3a44271bba2ae68ff0f1d6f473c573602ecb5fad2f6f634dd3f1719fd609c036188eebb93be61c7b740c78db83b043cefe1cb0fdd7d533bb4310ee81fd0665529b1405f9e2515e092d052192f1017d1389b385b1871bebd23469c96bf0e744852a925a3549cba9e3a96f43dd5477a533e1458c861f170c90ddb7e8a3e13bb10a138dcea6e0abba9aac91835c0e2396d5d2ee4b4190c62067f1c481bbd38efb85fa66532060747b67f9e6502fb26593c1649c14a8527a0a6bd7015b5b9eeda2d7b4fce79a1f0e18b457a21641303fd56064d9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "16f7c5076d402edb9f2c0d0a73f1045851a0908d692f2a2d3a5fd95da4cf27f3e9549ddc7b0e19a001401b6eecefc272aa1a611595300e9b5ae60b9f747c8cce149cfbe7f4c402ff5dc1eb7a5e90f8b7f9240110b746abb111359c5b78623d08ef7e3ad122c180c9e34ca16cf9b36d7eeb493c9219b6978dabc3678f4b5a9ed0606bafef78db1294839f928e5ddba036013e78e16b03a4f49d8828f36831b2a772e93f93b75484b777fce3b421f55039a1f1f69244e8dc2d6b379b796af19571969984bb07ce0f55583a5c318677ba84a64de56653a6053cdc11211e7e6c73ca3f88814b57f8b0c1096c4fe83615640622220e386f85a5927d3c5309002c69717464108f2b5801"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0139618d4a8e3cafaf6f87b21d88d8295013e8508ba87f9f02ffdb35daa80b5ad571acf56c5350315d14009533938c6c59ef8d26cee37dfe26dd3ba75a4c12434fc38f4ebbe03a43b8f08181b4d06277028bb73352ceb9239c6467c4a1178b0a6724c8d73609fe05a45f0000177fc5a377b7784a88484f0e057b2bee1a6dec07a7746f9e8ee6a5bcc0fffd3824edf510c656abb19a089c5566ac87c4e4d33d9d8b4e0544a12bd795b2dbc2801b39b89f501bd8416a38d3750f7917947778b779b8c5e7923dbb8eec24c4b38a57423706ce4836518892ef45b2f8a12e0debb5f9586e90501365c333bc07eba6836920c5c5a7390a6f07b153d28a47e8d8a334ffc8460134fa1b26977b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0323369d1f4a9668958bad24e36c2be20e47b26a89d4d3037d68789e825c71846f54258dc5a4d2d1d894c8496fdb71acb2d63afa97ce26a3c0f90fe58fdc0f38f099affc92dfb268e93fdd8962a8afdc4a0a118773482e3899034ea3f98dc91d041573a29acc74637cc3a5fadaef4de2d60b67b20ba5b4abc6458e592c954db61212561ea5c2bb3bb681d6e6b9119211aa1428e8c22b820ccd1d7031fb97ca1716fb475aa596f14c251ac3e5865359146ec20a2f4ef84b988c362e5c73c1616e0067842ef836761fdd23da0a00ad76fbc9c1b12086e1323d969446cb2beed23246ec38f16674992468dfc1cedc799ab73aef9f1d2819a441c8afe663e657eb09ebf186e16debd98eb7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0541d97b2b84f153c25561e677fa48c26218d47d2a69038b0575e0b2bb507e81e3d125d335a91ac635aacc011892ffaa69b9b11982d802d15dbe02954f7399406191f56335935bfb6e70f90d7a9419f4084b6ed9d730e27bbe900a6a58a8c20da5215e0eaa2c4a1a0f61429d1bfbf2a323d57a5f09bd4c360ff8c473a0b1a2e4ba41eba51a407dd57db67a8daf1fd0a13345177bad5bca13ecd8a6ae693234f7c7655c248e221dac222897271f089c545b9dca45625de984d497c9ad8b6a8c2778531dd262bbb7e729f5def7de782ee8d66f6ecc6ba745f5b16bfe67f47c158d8661be8472125c48da201a2b7808300f7c76fcc7b4c5f574ffd79edeb2197a660a1e7673a3862f26cf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "09f1f1e62ace5bc98c1b41854570c54c30cf0cd96a7bf7356a158c0515923e80d1cce2ba3131fbdeb8d247188b26230a4a02877a729fdce4d8a9972d1fd0af4740c34f23bb4941058e45cf249d8f9218c24754e3e2917da994abadf571b789cad51d41074fee1c714d012034b20bbd258d357b1d06c3376fe32c68093580fdcd99406bb4284153652c63e43115d6729812be9f6428318f7c23eb9cc3ef556918696aa535b002b725e3a1e122795e23d4f28fe8c1745d3d829882c7bc7db7a1a74df1d1d9f9a20ec26181c7f57eee9e176dd8049111fb363ffb493aa94e84e1ad575139decc3d9dec5ec59309f6d05219f0236b66c725d3b53a636db3b566fe936d50b59a73f4dad939"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0082bc995b99c1af472ed2c81440bd573726479b6d63e6baddb97fbce259f3873ffcec2ad96bc500356829e77abc7600a838ac9b387662f2b06c7653c7ca3fa50bad6a1d6518eb1af34781828a1a3a13992ddeb38f73e070d9ea0d2a342bc8c32efd98e5ee15721e9cd1c50770243d8ced5fe4c162f2cb0fc68e851e243dced19ec6d8ac03be899cddfc14871c18937265ae37b36f3dc374c794699b3347e9925c407b694c2541c8315126e69b0ea010c72426e83cb83726dd23fc1199b279fe200a785e256437524d20a07c8fd81f17abdb3ece6bdd67378c0e1186ac58a56cbc1aac6c07152a931d55462670df79b5c1502f92ec824abdd6e5faebc68b8420a03a19c54a00db6225"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "134b2a8e130d49d3cbc0d60e29a75bdf380e066e64902d955d950af0cf776234e0ff14b52e8b0a02134b7a1666db814213251285ae98373e905c587cdc7d98523bd54d347144ba070c6334c10b2f92b5d7d244df99f2009e30df035bcd41ef6b17ded1bc71d82a49b9451534713601cc646d4fb8771aae2cb9045acadf0c75f0be89e933b3eca66f955572d266b14c86dbc825ecda7d2f3078ef5d1a9a6676e6f9fc81d1e9bf50391917c778ae647886fde1b1445bd15ed8eb7138150e63f7fc6282b3aa9e407d921220ce87cefe410f6388b6fe1e26429b0929af911c630084a9f6e29d3754a79443d668b5e6a53e40104ccb5aaf6e3e9e72fd17df02285700dcdf8f127443869747"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "4f99ada92e63b19c09ebe9b49b1447a93a6c86032bec51d0dd8cebe8f3153d9e3487e7d29c47a41ee94feff54115d128b14fca7309fec87f3e54a6baa05e21a2c320dcc4c0baea03cce223a08ac6ddc876a7d66db61a59d25d7e38cdb37c3e59d1a2e34837d3b4457116c84a7e7b30daf8e9950b758d7ffa3432bcf337463d124004353628172698c41e9135dd77c23c9dfccce30652845d2accbaad8cdaf2857c6a086ba5ce494b659f4756c0e4b56793d2cdc9d12937545dc966398f9a3f4524e3c3f4801b19f1e42cc0bb8fe16d8605155f2178dd3d9008d9e9bdfd6ab8b71ab9bbcd2254650f4a5f0330b2ceb7441bc0a2f7609cf838eb39b8b1b109fa32ee04c083bc8c646b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "1bb439731b7ae12b6874783162b16dca098c378752625fc220fe53c2b66debc9ff7c2834af26a220bd15956b0f7b07b8d49926abb58ea5cd0442bbd9b4614ae7aef2a2d7e7223603bd736f43513524fabf0f958bd997e17cd6fc7c1cd6fc14577fb3691ff7e4273f266e7cbafd8e5494d96a057b4003b5330e7a52819ced9eaefc1b6eece210fbe7cb389ad87df6ad493d8fb7fe327e409721fbe1ecc99718896fda3e06845a466db1fbcc0a669cf42a7c29688f4da5d99f4e77cc6c6734f9ece0305ea20a219b28aad2a0bb0bd24a23b37f5b411cb92d38264c05ded2d96c0bae6905c864f38ef9d452c40fba1a3ee1591564280ad6e0561942355e341f1508a66f31c83297bf8137"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "5620bbda712c8e2664a235bd5a11e233509970d682650dd60fe1fd98c5c9b6cdd16aa79ea00fb7ee221730ad0a0df6709f6ffa699470f1e8270b346075dee913345f708fcdda9b5011327756d34d064c6e595794cb0215dd3d3121a515706bc3a44318b56587a02f9cf8770078f963d243d9f1e12228cef2ed22e961dcfe00cc552f9786fbd21d128f7eac7afa06a3abf6027161e88596699865b7623cf3107dcd4ac8a4bcf61a35583575dc3acd519c547ca8d15a43cb6a01a68a8614570ad497cc5f167a517b436a06241b68202255482116f4c1ffa9a0594191d367a5d22a2039f1ecd6320b8ec640c66905108fb0ec6b5d0d23ae386c1e883a7547aa65f76e22d5f182f17b6437"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0490f2692213ad6c3f56744b1b3be3d419f8dd020994282968ec8ab382cb608b62387f4452e778e2437f5e84a5277506aeb2e41b4dde83c1bea26761d61a19a5c0e460ce69697235fddb2b0cc4343b4b3ecfe4e802074484abfa5461383a85c9eb4daff7c84d45ef0475f9e6746c6b2c2925deb068bbbcc11c62f36b825e790b8224c3dd9325cd43a8e9138985e8e7446e914b23950567549a234e7f59ac45eebb5b5968d0bae4900b9291408fc654c9cea7f2b31d648b62e3f080d136cffb31ab1d4c2c134741e5ab77097e1c137641ddef70bd7a9fb520cc930f1335436d2cdd096b6b7ea816e1de56c45506e1579bfa4027ae066c27dbef4004ac5af682804322016b3b5cb33159"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00b4a85de3e2010095865ce98a144e4b88d42a1a6729155557b5d90d3ba80315a031b4e7480d6d4970c7764de7dc46dc1ee4c7c433826d7002fadf05c7f46c0bbc3298c66421742db268e77fbc022d42327969c34bab396590482f6cbc97dad8ddf778112a617629b559f9f4f6ee394b08080472d7dca1f35901b72b64262ad19150ab6fafc7327ef4be191712f475dcdf126db77422140f4e145d531e800230cf9fb4d481e152bd41d89b835366566e309602c48635e75f7ff14764900289959c5b0b42c0339e8fe4fae02c618fd07456f4210b0282bf93764b8d84ac93a7f056000a3a89b508fb99a0b6dcf610a31fb3d801d5bb3dafa52552aafb8012b8f8c651d2254c44c2753fb9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "222fe3024b8e5a4c11e233d1e00a75139968fcbaf8cfc9c5f0a8e2337a14c8052f126cb0dd7c784db408d3d75103854c7d275bd583855051e1bc1da13d39b2f34988d09d04455f077008e49088b6a2905ae270c7920f70110fcc6a4289a660722cecac7d32380f7c80e14062bcd7ae6ff0999f653d48ec6894ca28822c7328fa6d3cf06e98c9bd1b180a413134d217142b99e1ebff444406c96dfaac8c1cac8e5d5d432db84b0ca787a52320180e6704e9f367103c0f440b2740347dda3ebeacf7ee02e09cf93f916ba91e37e2b19a360b779944321083dd1d2e1185de2f7e12e275648cd5ea7a5ba44ebc8b2bb84693a6978632a65d08dcc89df3ef740809e981c41df9a9d1d50b07"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00fe578feddb333e11b3ee0a8870fd085e7b5739c5ac7d2f8d58614ea72ec6f6df92b7857dcda930ce6dbe69df4974a4eeb077d3ba18b2705e50a9412eac82dea651d26873584d829300f578f8f5ca0bc0631aa8a9448208f0586bbebd57b01cc72b4de6f6a22fcfdff3f0149069059bff027751028102b9c466da961217fd8ccaab0cf499bc9637f7b0768ec86e911e7907bce37723cc81678ed143ab3058b179ac68d07352b8a6a4a601510eaaee5ec3455d8b3ab8b4a12d18d49e431cefe69013419bc372f29106a01fc78c089461594c04b23f6493b3815ba91f19f99d01b903c93ec1565e21b0ec880b546330c9ed4b8e973ed7decbca936ac5d19ca26e051ea124267d0223013d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "462eee2ddc554144bd3bc464dc0fb2854597196d9d4c319a8c3ad7946256e7b2bfcce9455ab4c43d81f3ba2560626afb8e4ce2e3ce8e88dd7affef7aae67574956106e50f98c57dec0c375fe18b4a91ff6f451cde2fd5a4cde201466c393f01f7852bc89705db09e64f392ab438549d66b42caca88a65fe846970a76d590ac682a65d9863411b10f1a1e23b3a78b36dfaa63ecfe8818f993ccc3a04089cba1c778228c8200768a75b5c11a7c6c5d17aaaad01be3ede5bab8f5393cc18eafa3d00f46e326ec83c18436e210d86049ce5a0ffacbbd301b0e7b22d4d1b67b6f06e0cd9901642430201981e1fd0073695a771d5a648cca247a9d5a7df699afa74430ddf5fc9614fe59d7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0843859178a5aef0cf76fbf7a6fdbb261fb21f24aff24a97e37b3f87569cbbbd9eda5c7b5a56327e7b98892cd4aab05b53508eca2aea6e03d43b16bc2a5a857a93e9936e7cd3a128c3c12c644a5f24fd688e9d0daea4bc4a68dbbd5056ae5763dcde9ef60a8cb36a5705724f94572b65f0d210d0a8825d27b0fa2f1e3b6696e6ded4ef8adbeaf14dd34b3133469ac0895f5106c69caacc697dd18b2db57455087105f2fbde982501bb19a80d9bd2dc7083ae4249f4945961c65333e772726d08b4ed030a8942feb48ef390dd1b560c1270196c276ac7037863ae3d87b31c034e66561a620c271581a4f3054ee38ce8ba19094c49190ddbadab842adc5b62b4efd19de0f2f53c8ee9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "11193e4b78a93ecfad44551c9913729b5e5d57e903635ec76022ea9fb939f14ae7ab1c895c3f81c6d5a5f03abc7c778d9c6ba0aaba159af380a809963706568b051e8b0667b96bf839a8958cdcdb49b31119998145fb6ca0c7a2ff0e2c7c0d90d24c618257b9246b7c322fac9f92ca13d85df961c40d0afac9fb404497635bee784f297226e47f2ab0f6606263319a4ebacfa0abffbcb9781fa868d7038c66bde1afa539efdbaf175510f48c9d9c760fb64b4d3d34f10fbbfbc00e6d67b7697d68c751bae07f3970d7bf6440814b6ea570b6adddee2248b73d572256627197a907e285a6301fbad76013d41598786fb9d33e580d63950f2c904caf351d342d36d460763c7a1f1f5957"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "012f617c5a054349a3e93ed6854e6a540038bf60044262f151b9c825b36f564bf4079f6f9e5a818fa2a639ad7e6f638d3d9e201f3187150a66c0295088736ac020958974e0296ebd595ba4a0b63d6e11961736cf56e9270dadb49cd77bcb1a9d89ebfa9e2298aa7b3cdaffa3c675c5dccd222f0b38204287d1fc269c6440b4d39dc7b2b193b455747f75cc24c360da99df57c8fe2c3e92609e7507cfa3784a8c464c1a83e7b91e9ef6336576a56bb6636cbfb8ae1f6e9724e9c393576e5cfedb29e221550dd2f39e8511b0fc9b606225a49d5e5086502229b61989e3b6feebfa090e4474a325071498de0ecb789c0291a8ae6c04c57516d00e487f8a60417af209a60f0c34180337621b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "5553413dd9466d66d0ba3eac38d68084784bc3c8f04c414bf46c42971e3c09ead8da7e8290a5ec2a8d0c55b9a27af71bf44edc1d0f9c175cfdb92651c2fef12355de421ae9d4e463d97e1b7d5b75de136697f55729151fc55a4cf365f5f80cc0c65659bc29e0ac360c837ba8288bafa80b46b20f3598ab782f0d4dda7dcdec673d2ff7ab05868466e5fef586e63586b22100e92ef9bf3151ac53e59805adcbb011162a1963fa9baead1517eb8731ef4c84350e7d73ca47dfc5fadc1ff7404de4ffb2d3203b1e61f525c9eee55672ac1526730bcb4b4d4b54c2707b876710d0dc0f2ae38b813f4b738e381bc4e05c4cf0b0f4db1dfce5107b6d4c1199dbe5ab89bb063a591fa0c879"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "410a34783e03400a52d847bea34213583afff783d80f3aedbeea0807d347b680ca6a72c463ab1b2779308732b486309864311e4a542b46d2123033fabd0bdb958890cb4b9ce864e091ad7edf2f406f6a542f2d34eeb21b8155667ed89a7fa9434fe479a79356ab10155e097a329c3b11a06baa8dddeb4fce45d01254c6bd3ad03881f2344c3378dd8fbb69feadf8f2bb0103d164a9ea52214f92e0fec3377530610f16247693245db19e2c6fb9ceae1e45be90a2028bb6fda1038483099635d09e2bb794fb5fc9140c1f775ac10d094e6400f9ea767f5cc2b33f47fff7d67fa00e00dadf8fb4f93f1d38d695d9dafaf80cda2176b2585e4be4983025beb099143dc541c7a9c268b1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0101d935080df531c252bf2d2116719dce5eb155862ad5984fb325c75e00b9825ce9fdcdb39db9a50a4890920d69815849e62d478584417161ae96a98c306d080efa88356c23ce7d7b1e17b9c99987391a7adc0121a9bbf2df21dbfc7dcf4cb529d9c3c42975624cb34750051a8eb6b0c2f4e89c904ec614958a071a5cb1acb2a5f9dd6249231fc7a167d84bf2f187642eb1df50db6824b54c50336d1a3b189c42aa74c613667620849b5fc38f9bdc8b8cd07e0b2d4597ad6ceef65d6bdb03bea4fc62135d28bdd5f47f7403fb17924adec5e138d4cbbbca0efc141dd7cac26259ebfc9ecfe1473c6c9408309767308a83afd1a3294a00f258989e1049a9a984549b6d98e4e1248acd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "1a5fe69fe95b9137923fc26b908e86bbba2de8456cce10b64574d78ccaa5534f19facabae2320f307bee35de24d929ab0fe6fc94bae4567787ceb74263f656c850a4395188a4a84690aafb5de68543c63bd0ee3cb16669aa37178d8b77dc019a5180233eb13245464babbf180a1254fe2e5d430dacf7dec7548d2d75ce795fd1979b0951413807227f8104a01a617fdb8a0b6605ae4ba2d57a55c4121999526d8e79585313199b321ee6e1956573853c2a80df7111f2b28070cd0646e66952c7344845d218924347568ab3f5270503ac9e85cb8a20944ad8b6811f90cca6c8a9d9dd29fe747296a818fcaf50813044cf5970293e51df42ffc5fcca42987c0e2dc6e0c459251166ed67"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "028a4f88abd93fb89bc21645e872a7887f86e6e484f935cca6c8b4c0d5d40612b3de89a931d97b437737f6ff6a5f617c95061103a13ce3b4f145c13e9712014b9c2130d87e2eb4136fa2e53c280327b5bb6ae4c4bb0be5a3c265f987cd2b4c1232b1eb93deee4bc8298236454394fd25e1f8a54b778997b5fae1809382c308738b8f011c812ff155e9a52dbea46ce0e1682d78fadeff30f24bce8c998bf2bcfdf9487e7e938e746ad1f73178b5803964451234d06a8a84afa9fe0b7e3507892fe07fc1c36667761386d289a053b7b2df38fb4081e5ba407a9227233135919f25a30ec20d1494fb7783b75375aa17263e9c45d4667dc34bbb40116fcf67e089f72f46b8e567a9aeebb8c7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "008eb22ca7dd25fa7c230083277360e3e29e3ddb30206f48037c4b31ce2ee37bcad6fc314e4a27192d945efc1ef02bb5a514b82baf98b442b3155577a689c195fe3a59db2c21204e00ff09d40b34fa7bb6cfb61a727f61b024fe5fe182e3c051f07776b3856b9af0c24bb887aa5fa1a02143ffe2349d7555bb4ba1e3b597e1cd082c276d7abcaf8c4ff951a594e3a0ae3460f37ef0b651598cbf3e991de6e44e76b4ecb1d76e87c2f6011a9262d3f1a2790a534f8bc3f49119d7afb8130aceba030397304b1582fdab9b00a76ce97acdf7055faa2d85a037f9601ea697fd2bf7f6fc600aecb5a43f18c8fce7eb1221a8d80e0e8c8abab3ba3979d98654c9db24114d0f2204b5914af377"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "02581eae696be714b5d7ebbd232660202fbd411dc9ddde7d73a48bd7b692d0b0c216b3c39e51a707ac47b6a5dc39de6b249cd0b0b27e67c525af2308d5e9763a572f5b2d3fad53a8567ea0cdac5dd3f90a41c05cd79c89e9c94d592b8ffbc7c1f07a2e7bd820dd5ce55c38e21c7890855bb3871e276be9305ae02867697d32eb3fdef35e3e9cc6f8a1ea3970f0b43f520849577143631f3414914a7af4ed82faa3f3540db1d6c249f5b642c83301055326a4ada3579eb05994b2d8c89139aed660d577337ac3a1d6a1321e58a3808aaa9f78954880d11b2b1b6b3f1b65a36ffb25914ec28d3a7bf707ba663f49e067ccc2097125a12e993f13899ce0fef029d99a6a3faf8a54aa1d7f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00a09c7c19c421b546cc94c42e8b8662b6ba80d5ec245157378b366ff039cb4d29a037993976dadded5cbb9baeae4443ae6fbd8e0dbec9da4237e68f83dfd6fdf2dd1349fba3a6cbc92832177837eea57e0271c6f5a1bfcbc8082eb14a0ec5b230854dd52fb0bf1d56a3575b2bc518e28853d87997c8a36ff1572a128e4fedc64ca012294751e0bfa9cf97c878e049144ed49d562cb70051c5a91fef5b4c4a55d3284a4093fadefe283d668e62bede1b84dbb1fb55137a34979a0d71029a71e3a6f6cb24db24cf950aeb182ea9acbca214170f3e3749677ec688cf47b35b9820e6e1e1731503c3de31256e13aa641f98fc48ca39249d1b8332d596566116d59488a87e5f1ef8adde9163"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "0a98538007d0939057ed39cf18a18e8f1fca061c849099a771aa5dbe42bd08643c42433d2f24a9c598e537a7b5448cb92af1b605ab97b4c817056e7cbad446120a5b91ce77db2244ebb821ef727977d731e8b080bb71db4111deb82815621bd378696ebc2679e46f764b7946c91c37ce4fe6ee08e8d6db3d4e9223bfc48485af6bedfc0e06c546ac208ccb9153dbbc47d1dc56974aa6faf1235cd691059aaadd4fc6bc6c5bdbb2c1691561f84f0436c7f3df2153ab01aaa33b214779b152e6d3422b5a98931a5a08bd6db1a9c23ca1eacae57d0cb0d66b9e405fe71dcc1ab5d74e3f17800eac2ed9624a0a66d71a08087611db077a97a3131ea6e9f4e0ea21b27812654984e8220807"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "023481defdb2d2224a9b315baa087dffdbdc511c8c581e43d7c6209d8b7e0eba75cbbc38b6b676dae33c46abb511064ed248cb2296b0246d1805f8c141e77d0bce4608fd9c9519bc0a827aa9ef5f877a658ff79efef187beb25adec65c2d8b8928549f534f25e90def3c4442e58cb5683c532d1f6ebc9d63e1da8d4fe012880ad1b357834d426ce5a0c6ccf09de4a917a59441e7977d46d3ad83e96af125d02cb7dfdd78047ffea7157763f669aa399937b05b46362e3b6db6d08b6e03c6b334ae6d7fc11bff839e155625e7e6b5eaf9372b086bcea47caded8e4eff4756df5e3904bf9facc73d9cec4a4a8f0ee1665752f5743710cbed7d0ac01e9a8253d3310f37c1e72bdf2785e685"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "4f41312cc1981a2007eb29a87fa2f3ada1915937e4f0807454ba5f18a03a32f954f45ae8c2df13ec3dc2045d92cca2d3135d25100d7b1c4f5c334513e9ecbf02f27ada2f6d6023ad737148840609befcf9edf2bb574782209f5a712d5b35f6cc3373580e89c6148bf2cd435e4f401d678940cb9226e8da87a33d992da172da41f710410cd450a1dd2a8608033efd87f5a23a5af1a1158feeb422e7391b960139f0cd7401a52c747803500dc96f130e2ce224d6b60169b210401b35fd3eab56bbf0887c6af3acc54473f928adb538ca44702789ac934151bb0013368883c9fa9a44405039a37a837d9ddc35d2d019e102ec3e1ca5c9e9cc0b317a3b6b107aaa657b8bdab892dd19e7a9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "029116984040a488254941b331e73395ab58a597fea5e15bd92042b3e2b51e883789ec6ae538aebe3d84c928952ef9f239bccf3fc267a85d81c2ddd209e28f56e62dd9996c6bf56f781b02ee3a2c5c05e711ac28a48098cc1c4b9019206510271388db9550b46dcea014e90b8e2a37b5f311af9167343e3dbaad5717bacd0a3e8e2be61bae423f9ac9c9d479ace3bb5e65ccb4560b98325d6fd467e533b2a39d2982a22dd237f7ce38b4c1d6540636a1e39931a9f441ad5ee3c1d06e5c4eb42db1a965e66b59b7642fba3d49eca9c2287f33e6815954c4bda231fe5c60e46d003642fafc3079a2e79da673833b1b6895f795ecc1209f5d73b7131f0c68514115908353a7930b301184d5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "09c877d3a8b3cfaf8b5192de0ebc2ac29c9306f49eb83e6cdf1839f56e46529a72570354ace0e538d630d10cc715c03e1db78a70a0b5b0871a41408dfb16c2b207d1cd59ff28f2f6bf835a1d0ca72ba144630f71e521f1671dfa16bafa9aea3ce69b9d5c53e4c1c546d50b4f1ef847d7b9d256af9da829b5529d249dd4a62784a2b0a01411a42d30358161acff7ae2ba69cefe7b4efc7370f49ca482f1832c578450718d014e6d08de7b1412bbfae9e5123ae2218e37063b480945325005671f67f3d0305ac02b1813b02342fe17113972c796a3d14fd173a169cfe7a05f389a90f78f5bafd13550bc3778d12bbc00195a626550fe40b763c3fd713fcbdbe786cc0af7a2424d9226c9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "1508d017f18635efa2cc25c56f40992ee2aed4868ffe810d1973741a4d57390a4af51e8431d90f65b442ff504300bcc54fa0bcc5208a8dc59ca4fd2d227482ea4aca9fe2ca79d1eaff70054c0b0abb794dad9bc3900779446283c95aff750770f9ddedb1a66b0deb264c39a8460329585776546d934094465517a802adbd15aed60aec187f52a20dcd5b0f023480cb784b12248374271f4ca43f276d7c911f2bc5ed25eed2aa2c9a408284652fca768eccff61e7eff25830e660c9110f78c325d2374fcd727ac5739886adc1d7908a07fd803a08b7bf3ff7d8f55fd76668886f85999eff5a0eb704f746ab1357bc942d57a37fc558abef850984ae5bc46d771781117533c7166b54f7c1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "014160a4d755f7f9ea44b90f75a2c370694d3f5e0715a1146bb309a04aea9169bb6117dbc43c90381de4d7fead0177179241be84f468ad38a629068108288e4c99e2c9ec86ccfa11dd62a22aacf486cf5c23f03162dc0d981e705837b8bec7d1b13123af2b331bd3bc9b2054f59f317853ed63cc0af0dc4ef15fc751a43e83f731cbec2191ff6f3990e8787ca0e4e47793ff7bbec63c012e07f24647b484fe8d0da12215c0d5cdb0ed2d44a253d5c825de3c42bdab327260c300d0d806695a717b59fd352e68f9e0828d7d546a57333578d1dc0e2c48791dc2a659fbdfaada59c7071b1b440eca073697bb7ddfa3f98131e23430c17ab6d4e34b44b11ed1b5a1cb7918b99e7bd713ab"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00d050af27a2af0aab019083053bca9a2318f1d3a322852073c21fa9109c7fe2fe5a16126ea0064c7655ebd9b1a67c9a61a028117fb9df03cbf774aa4a341f278f1570bbc0b3efafaa855d3878ab6039b2ffbd8c5f1fb9f04cc282d072eaf7904f5feb42b950b4236da9e67b7c5f4889533ba66bece01c0c35cafdd5b2b310d91173ddbbce856f5c4671c0f61b28defde2eadc7d6b96229e4dd12baecae8680aa038c104df148191a40e5f0cca2b25b456957bd8f2145529e71d25762fedbb3b6cf3023dfacf47200b91b6a4bacbfaa92ffaf4a760ec132868b9e7e3f3d0f7cf77a1426645ad54a2e057fc01e223682e7c56afeec356d4f53a08528e5d2684b8be5eab78a3d9b46cf331"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("Carmichael number with 3 prime factors\n");

	if (!BN_hex2bn(&bn_tmp, "00f4e8aaa62114c404219ed23f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("S. Mueller, Strong Dickson Pseudoprimes which are not Fermat Carmichael Numbers\n");

	if (!BN_hex2bn(&bn_tmp, "07ff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "05361b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "14f5d5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "008a8d7f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "018271b1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "3e9de64d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "00bfa17dc7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "011baa74c5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "518dafbfd1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "01053cb094c1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "323ee0e55e6b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "1c6b470864f683"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "081f23f390affe89"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "00ffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("bound for deterministic tests\n");

	if (!BN_hex2bn(&bn_tmp, "02"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x02\n");

	if (!BN_hex2bn(&bn_tmp, "03"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x03\n");

	if (!BN_hex2bn(&bn_tmp, "05"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x05\n");

	if (!BN_hex2bn(&bn_tmp, "61"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x61\n");

	if (!BN_hex2bn(&bn_tmp, "65"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x65\n");

	if (!BN_hex2bn(&bn_tmp, "00fb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x00fb\n");

	if (!BN_hex2bn(&bn_tmp, "0101"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("small prime 0x0101\n");

	if (!BN_hex2bn(&bn_tmp, "7fffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Mersenne prime\n");

	if (!BN_hex2bn(&bn_tmp, "01ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Mersenne prime\n");

	if (!BN_hex2bn(&bn_tmp, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Mersenne prime\n");

	if (!BN_hex2bn(&bn_tmp, "7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Mersenne prime\n");

	if (!BN_hex2bn(&bn_tmp, "07ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Mersenne prime\n");

	if (!BN_hex2bn(&bn_tmp, "02611501"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "00f17a60a5d627ded85b6a9a397c2ba63bb27910ccf7e3135d4d1ae8c9f5cc1e4bf01ea704abb2000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "01e764f3171d1e44a5f0c50c6537730168041cd93fa34898140da93d3df2939adecf61802daa63eaf08428d72148d63f267f22bd24cd411b7f25984b057bda5c11510000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "3a7c596683f12898e64bf1355bb9bc85f01d91307e568d01afdc9cf0b3fa9e464b140d899d9bf62a0c61c2bf0a8bca1de36f6d36a5be4aa212681896def96f583c8a7cfe362b4e823bd244f813e575391a029df7012e738d3e2e8e0181ea40000000000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "72b20ce22e5616f923901a946b02b2ad0417882d9172d88c1940fec763b0cdf02ca5862cfa70e47fb8fd10615bf61187cd564a017355802212a526453e1fb9791014f070d77f8ff4dd54a6d1d58969293734e0b6bc22f3ceea788aa33be35eed4bdc1c8ceb94084399d98e13e69a2b9fa6c5583836a15798ba1a10edd81160a15662cdf587df6b816c570f9b11a466d1b4c328180f614e964f3a5ec61c3f2b759b21687a122f9faefc86fe69a3efd14829639596eb7f2de6eab6b444d06233d34d0651e6fed17db4d0025e58db7cad8824c3e93ed24df588a0a4530be2676e995f870172b9e765ec2886bce140000000000000000000000000000000000000000000000000000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "2c47a7947e4ef970e990c8b4a793b5f7d49b9af95a12b9f08475e1cf58f31046fd224c3ef20a736d7cae39a2f989d934c2aa644483aa6e348bd41c34a6819d7c08fdbd93a7f7c24a4756bb7dd97516287e161af87e56735c06d61918cb2fd4ae9dc1c7f2cbb5749934626af5f4db5bde6b748072c004110d45f6db0fe51c4889ff053bb2a24f83bbb80798b94e5d7a189599d85792807626de78a61a7468eab70a2c4dda6200e0c8328408e0327897220bbe009cf8bbbb23fa1cb5fbd3713f7172f8186d059d0b97c2ef5b096c558ec61f66e81116be44f2940f4c93b67d7cd3564c266540fbf0bb95cc3c52c9dbc71aa6a424457131aec3285e6ba46e828e635f3455e30b6db3e4680ba04c580fb569145f6371a0d352f40321751cd26623e92a6c5c9e83eb655338c9077826148e23c3705b8f11b15a00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "1774015499125eee9c3c5e4275fe37ffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "22d4fb39eb23880b4674bcffd06a18547ee73e7e77f1fb29c0dbfa66ed52cb8b22bbe0ed9b2a2b779c9037d7b412a389bec5ffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "371196ced90a51b120fd9171fa388fe8c2e634f9ee10f4bcddddfd269ebda2f3eed661eaab3cfbe6914395a73735701d7d65e278f76842b02d1f8f5d941d652067ead60bf9bb537ae7e13404711ee80b35bbf5936641be34d53d4b3bbd025bed4be7fa44113cfea3ffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "0120dd73742e20e30f56d82ace2d9ff917e66b2c92024a1444490511d41a39685a9901187f206b5a248b9e52d82f15820801be21beb73ff9e0c8150c69334f02fe9593493b55d48229601857a3ca4449a444d2c0566936deadacd46310d04480265834fe9b5e733357b0c73a0d1e23d85e401e8c3b60571045a6bfb1a19f4940140736098dac2d705dc1339370f1ac19252b931c450bb260800bb40aa404dc54199b7251abcb50d26fc9de82de037c3b9926a2958bd6a1d8690805c0681f5cb5d90b1447cb7e5d81c436b913d743372be382e3bb2d1cd7185948136957af2496888060c7b7ea519b173d5f190c27c70f3dffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factorial prime\n");

	if (!BN_hex2bn(&bn_tmp, "07c97d9108c2ad4329db02eb8f166349"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "010001"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "663d81"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "00b161194487"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "08112264cd9bb77f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "0b73493decfd9b68318ef9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "3d30f19cd101"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "126cf51772d253cba3f5a7cf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "00d3eafc3af14601"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "013540775b48cc32ba01"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "3a294c585a8f5c7073e36ee3637cab2586d049baa0ba2c911801"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "03f1cb0fdf0fbef0f3747f239f5a8983e72b455488b792c8e29308f8c78e7f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "37a5f7f30fd2d1f46cd794e8337106ccebced1189c1f5b6b3c525b64b6c36768785f7912013f"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("Factor of Mersenne number\n");

	if (!BN_hex2bn(&bn_tmp, "00c4ec4ec5"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00c18f9c19"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00a08ad8f3"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00fcfcfcfd"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00c71c71c7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "3d70a3d7"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "3ef368eb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "69d0369d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "51b3bea3677d46cf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "7e3f1f8fc7e3f1f9"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "43fa36f5e02e4851"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "3454dca410f8ed9d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00c5b3f5dc83cd4e93"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "593f69b02593f69b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "008f6ec07432d63dbb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "101767dce434a9b1"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00fafafafafafafafafafafafafafafafb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "0c934ff1a0c934ff1a0c934ff1a0c935"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00d2f87ebfcaa1c5a0f02806abc74be1fb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "7880d53da3d15a842a343316c494d305"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "6a850096a850096a850096a850096a85"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "0098dbdea62334302c77d10fbfc4b593eb"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00df0041ff7c0107fdf0041ff7c0107fdf"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8af8b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "7f26fe4dfc9bf937f26fe4dfc9bf937f26fe4dfc9bf937f26fe4dfc9bf937f27"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "009b8f4f9e02732385830fec66e3d3e7809cc8e160c3fb19b8f4f9e02732385831"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "64a9a50bc0a383524478973fdf4c22bf1b14f339bd92a6942f028e0d4911e25d"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "7f95438b41e0500d578e97c3f5fe550e2d078140355e3a5f0fd7f95438b41e05"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "5f1bbd6c9500cae5d85f1bbd6c9500cae5d85f1bbd6c9500cae5d85f1bbd6c95"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00967300c9a633fcd967300c9a633fcd967300c9a633fcd967300c9a633fcd9673"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "00a305942530f7f11f9cd2c027abb32354eb8b77a1c8368c165094c3dfc47e734b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("edge case for Montgomery reduction\n");

	if (!BN_hex2bn(&bn_tmp, "feff"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 1\n");

	if (!BN_hex2bn(&bn_tmp, "ff3b13b13b"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 2\n");

	if (!BN_hex2bn(&bn_tmp, "ff38e38e39"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 1)
		printf("negative of a prime 3\n");

	if (!BN_hex2bn(&bn_tmp, "ae4c415c9882b931"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 4\n");

	if (!BN_hex2bn(&bn_tmp, "a6c0964fda6c0965"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 5\n");

	if (!BN_hex2bn(&bn_tmp, "ff05050505050505050505050505050505"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 5\n");

	if (!BN_hex2bn(&bn_tmp, "ff20ffbe0083fef8020ffbe0083fef8021"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 6\n");

	if (!BN_hex2bn(&bn_tmp, "ff5075075075075075075075075075075075075075075075075075075075075075"))
		goto done;
	if (bn_is_prime_bpsw(bn_tmp, ctx) != 0)
		printf("negative of a prime 7\n");

	printf("tout est bien :)\n");

done:
	BN_free(bn_tmp);
	BN_CTX_free(ctx);
}
