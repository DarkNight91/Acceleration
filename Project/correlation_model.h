#ifndef CORRELATION_MODEL_H
#define CORRELATION_MODEL_H

#include <stdio.h>
#include <math.h>
#include "iostream"
#include "string"

#include "Eigen/Dense"
//#include "Eigen/Eigenvalues"

//#include "circuit_class.h"

using namespace std;
using namespace Eigen;

#define DEFAULT_GATES_IN_A_COLUMN 2

class correlation_model{
    
    friend class circuit_class;
    friend class gate_class;

    public:
        int rows;
        int columns;
        int min_x;
        int min_y;
        int x_grids;
        int y_grids;

        double **corr_matrix;      // This one for length
        MatrixXd corr_matrix1;
        double **corr_matrix_width;
        //MatrixXd corr_matrix_width;
        double **corr_matrix_oxide;
        //MatrixXd **corr_matrix_oxide;
        
        // I will use only three different parameters for the variations which are length,
        // width and oxide thickness. I can use a three dimensional pointer for various components
        // but for the sake of clarity i will use three different two dimensional pointes for the 
        // three components
        double **random_numbers; // This is for lenth: Not changing now so as not to affect the previous coe 
        double *eigen_values;
        double **eigen_vectors;

        double **random_numbers_width; // This is for width 
        double *eigen_values_width;
        double **eigen_vectors_width;
        
        double **random_numbers_oxide; // This is for gate oxide thicknedd
        double *eigen_values_oxide;
        double **eigen_vectors_oxide;
    
    public:
        //int*  parse_placement_file(const char* placement_file_name, circuit_class *ckt);
        void generate_correlation_matrix_quadtree(int x_grids, int y_grids, const float *corr_levels, int no_of_levels, string param);
        void generate_correlation_matrix_grid(int x_grids, int y_grids, const float *corr_levels, int no_of_levels, string param);
        void print_matrix();
        void get_eigen_values(string param);

};

#endif
