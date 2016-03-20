#ifndef GAUSSIAN_RANDOM_NUMBERS_CPP
#define GAUSSIAN_RANDOM_NUMBERS_CPP

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gaussian_random_numbers.h"

double some_randomness = 1;
double initial = 1;

double* gaussian_random_numbers(int no_of_iterations){

    double *uniform_random_numbers_1 = new double[no_of_iterations];
    double *uniform_random_numbers_2 = new double[no_of_iterations];
    double *gauss_random_numbers     = new double[no_of_iterations];

    //some_randomness = sqrt(some_randomness)*log(some_randomness)*pow(2,some_randomness/1000);
    some_randomness *= 2;

    if(some_randomness > pow(2.0,32.0)){
        some_randomness = initial;
        initial += 1;
    }

    srand(some_randomness);

    for(int i=0; i<no_of_iterations; i++){
        uniform_random_numbers_1[i] = (float)rand()/RAND_MAX;
    }


    //some_randomness = sqrt(some_randomness*1.5)*log(some_randomness*2.3)*pow(2.5,some_randomness/1000);
    some_randomness *= 2;

    srand(some_randomness);
    for(int i=0; i<no_of_iterations; i++){
        uniform_random_numbers_2[i] = (float)rand()/RAND_MAX;
    }

    //Using Box Muller Method to generate Gaussian Random Numbers
    float U1, U2;
    for(int i=0; i<no_of_iterations; i++){
        U1 = uniform_random_numbers_1[i];
        U2 = uniform_random_numbers_2[i];
        gauss_random_numbers[i] = sqrt(-2*log(U1))*cos(2*3.14*U2);
    }

    delete uniform_random_numbers_1;
    delete uniform_random_numbers_2;

    return gauss_random_numbers;
    
}

float* gaussian_random_numbers1(int no_of_iterations){

	float *uniform_random_numbers_1 = new float[no_of_iterations];
	float *uniform_random_numbers_2 = new float[no_of_iterations];
	float *gauss_random_numbers = new float[no_of_iterations];

	//some_randomness = sqrt(some_randomness)*log(some_randomness)*pow(2,some_randomness/1000);
	some_randomness *= 2;

	if (some_randomness > pow((float)2, (int)32)){
		some_randomness = initial;
		initial += 1;
	}

	srand(some_randomness);

	for (int i = 0; i<no_of_iterations; i++){
		uniform_random_numbers_1[i] = (float)rand() / RAND_MAX;
	}


	//some_randomness = sqrt(some_randomness*1.5)*log(some_randomness*2.3)*pow(2.5,some_randomness/1000);
	some_randomness *= 2;

	srand(some_randomness);
	for (int i = 0; i<no_of_iterations; i++){
		uniform_random_numbers_2[i] = (float)rand() / RAND_MAX;
	}

	//Using Box Muller Method to generate Gaussian Random Numbers
	float U1, U2;
	for (int i = 0; i<no_of_iterations; i++){
		U1 = uniform_random_numbers_1[i];
		U2 = uniform_random_numbers_2[i];
		gauss_random_numbers[i] = sqrt(-2 * log(U1))*cos(2 * 3.14*U2);
	}

	delete uniform_random_numbers_1;
	delete uniform_random_numbers_2;

	return gauss_random_numbers;

}

#endif
