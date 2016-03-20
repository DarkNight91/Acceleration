#ifndef TIMING_ANALYSIS_H
#define TIMING_ANALYSIS_H

#include "define.h"
#include "parser.h"
#include "correlation_model.h"
#include "circuit_class.h"
#include "iostream"
#include "stdlib.h"
#include "string.h"
#include <stdio.h>
#include <math.h>
#include "map"
#include "bitset"

#include "cholesky.h"
#include "matrix_multiplication.h"
#include "gaussian_random_numbers.h"

#include "Eigen/Eigenvalues"

class Timing_Analysis{

    private:
        circuit_class *ckt;
        correlation_model *corr;
        float size_of_die[4];
        int **adaptive_blocks;
				bool only_FBB;

    public:
        void Init_Timing_Analysis( string benchmark, string placement_file, map<string, vector<float> >  gate_parameters, double required_output_time, int x_adapt_blocks, int y_adapt_blocks, int *actual_blocks, bool only_FBB1);
        // TODO: The following function is commented as of now as I have to discuss the 
        // adaptive blocks with the other two guys
        // void Create_Adaptive_Blocks(std::map<string, int> adaptive_blocks);
        //The monte carlo method is special and it is used for comparison with the SSTA
        void Monte_Carlo(int no_of_iterations);
        double DSTA();
        void SSTA(map<string, pair<double, double> > &mu_sigma, map<string, pair<double, double> > &arrival_time_slack, map<string, vector<float> >  gate_parameters, double required_output_time, map<int, bool> adapt_block_exist);
        void get_pca_params_plac_info(map<string, vector<double> > &pca_params, map<string, vector<float> > &plac_info );
        void get_probability(map<int, vector<float> > &probability);
		pair<double, double> get_power(map<string, vector<float> >gate_parameters, map<int, bool>adapt_block_exist);
        void get_timing_yield();
        double get_timing_yield_temp(map<int, bool> adapt_block_exist, double threshold, map<string, vector<float> > gate_parameters );


        void calculate_all_possible_cases(map<int, vector<bitset<64> > >&cases, int no_of_levels, int no_of_blocks, vector<int>&act_blocks, double threshold);
        void third_method(vector<int> &blocks, int no_of_levels, double threshold);
};

#endif
