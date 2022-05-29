#include "baillie_psw.h"
#include "crypto_utils.h"
#include "miller_rabin.h"
#include "tests_bpsw.h"

int main(void)
{
    /* unitary tests */
    miller_rabin_tests();
    jacobi_tests();

    /* test suite */
    bpsw_tests();

    return 0;
}
