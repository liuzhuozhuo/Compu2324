#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, const char *argv[])
{
    int i;
    // GSL's Taus generator:
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    // Initialize the GSL generator with time:
    gsl_rng_set(rng, time(NULL)); // Seed with time

    // Get uniform numbers
    printf("Uniform random numbers:\n   ");
    for (i = 0; i < 10; ++i)
    {
        printf("%.4f ", gsl_rng_uniform(rng));
    }
    printf("\n");

    // Get binomial:
    printf("Binomial with n = 100 and p = 0.1:\n   ");
    for (i = 0; i < 10; ++i)
    {
        printf("%u ", gsl_ran_binomial(rng, 0.1, 100));
    }
    printf("\n");

    gsl_rng_free(rng);

    return EXIT_SUCCESS;
}