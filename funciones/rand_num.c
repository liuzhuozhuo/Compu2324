#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, const char *argv[])
{
    int i, n;
    // GSL's Taus generator:
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    // Initialize the GSL generator with time:
    gsl_rng_set(rng, time(NULL)); // Seed with time

    printf("How many random numbers:\n   ");
    scanf("%i", &n);
    // Get uniform numbers
    printf("Uniform random numbers:\n   ");
    for (i = 0; i < n; ++i)
    {
        printf("%.4f ", gsl_rng_uniform(rng));
    }
    printf("\n");

    gsl_rng_free(rng);

    return EXIT_SUCCESS;
}