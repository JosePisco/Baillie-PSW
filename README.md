# Baillie-PSW
The Baillie-PSW primality test is an arithmetic test aiming to determine wether a given number is prime or not. It lies in the category of probabilistic in contrast to deterministic ones. It relies on a base 2 Miller-Rabin test combined with a strong lucas test.

Note that its reliability has been proven with base 2 Miller-Rabin but in theory, "any" base could work, they just haven't been formally proven. No composite number under 10k digits can pass this test.

## Implementations
Here you may find three different implementations of the BPSW test: one in `python`, one in `C` and the other also in `C` but using bignums for real cryptographic purposes (the python one being way too slow). It uses the openssl Bignum library.

### Python


### C


### C with bignums


# What more ?

For more informations, check my conference on the subject: https://www.youtube.com/watch?v=MMX7LeDe9GQ (exclusively in french) and my slides at: https://www.lse.epita.fr/data/lt/2022-04-05/lt-2022-04-05-martin-grenouilloux-carmichael.pdf