#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <stdbool.h>

/* PROTOTYPES */

/*
    Miller-Rabin primality test, exclusively for base 2
    It has been adapted for base 2 test and shall not be used
    for general cases.
*/
bool miller_rabin(int n, int base);

#endif /* MILLER_RABIN_H */