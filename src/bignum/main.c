#include "baillie_psw.h"
#include "crypto_utils.h"
#include "miller_rabin.h"
#include "tests_bpsw.h"

int main(void)
{
    /* test suite */
    bpsw_tests();

    return 0;
}
