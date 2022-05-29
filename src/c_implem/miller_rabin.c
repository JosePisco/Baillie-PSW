#include <math.h>
#include <stddef.h>

#include "crypto_utils.h"
#include "miller_rabin.h"

bool miller_rabin(int n, int base)
{
    /* base 2 only for the Baillie-PSW test ! */
    if (base != 2)
        return false;

    if (n ==1 || n == 2 || n == 3)
        return true;

    if (n % 2 == 0)
        return false;

    int r = 0;
    int s = n - 1; // always even
    while (s % 2 == 0) {
        r++;
        s /= 2;
    }

    /* In general cases, the Miller_Rabin test has more bases included in a loop*/
    int x = modexp(base, s, n);
    if (x == 1 || x == n - 1)
        return true; // TODO: not sure

    for (int i = 0; i < r - 1; ++i) {
        x = modexp(x, 2, n);
        if (x == n - 1)
            return true;
    }

    return false;
}