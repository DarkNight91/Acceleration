#ifndef MATRIX_MULTIPLICATION_CPP
#define MATRIX_MULTIPLICATION_CPP
//This is a s file which multiplies the random numbers and the upper triangular matrix

#include <stdio.h>
#include <stdlib.h>
#include "matrix_multiplication.h"



void matrix_multiplication(double** random_numbers, double** correlation_matrix_UT, int size_m, int size_r){
    //One thing i will do is to replace the independent random numbers with the newly generated random numbers
    //The other assumption is that the matrxi which i am receiving have proper structure.. I am not checking here
    //whether are can be mulitplied or not!!

    double* temp_memory;
    temp_memory = (double*)malloc(size_m*sizeof(double));


    for(int i=0; i < size_r; i++){
        //intialize the temp memory to all zeros;
        for(int allocate=0; allocate < size_m; allocate++){
            temp_memory[allocate] = 0.0;
        }
    
        for(int p=0; p < size_m; p++){
            for(int q=0; q < size_m; q++){
                temp_memory[p] += random_numbers[q][i] * correlation_matrix_UT[q][p];
            }//end q loop
        }//end p loop
        for(int rep=0; rep < size_m; rep++){
            random_numbers[rep][i] = temp_memory[rep];
        }
    }//end i loop
    free(temp_memory);
}


void matrix_multiplication(float** random_numbers, float** correlation_matrix_UT, int size_m, int size_r){
	//One thing i will do is to replace the independent random numbers with the newly generated random numbers
	//The other assumption is that the matrxi which i am receiving have proper structure.. I am not checking here
	//whether are can be mulitplied or not!!

	float* temp_memory;
	temp_memory = (float*)malloc(size_m*sizeof(float));


	for (int i = 0; i < size_r; i++){
		//intialize the temp memory to all zeros;
		for (int allocate = 0; allocate < size_m; allocate++){
			temp_memory[allocate] = 0.0;
		}

		for (int p = 0; p < size_m; p++){
			for (int q = 0; q < size_m; q++){
				temp_memory[p] += random_numbers[q][i] * correlation_matrix_UT[q][p];
			}//end q loop
		}//end p loop
		for (int rep = 0; rep < size_m; rep++){
			random_numbers[rep][i] = temp_memory[rep];
		}
	}//end i loop
	free(temp_memory);
}
#endif
