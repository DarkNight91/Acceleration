#ifndef TA_H
#define TA_H

#include "ocl_opt.h"
#include "CLHeader.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <Eigen/Dense>
#include <stdio.h>
#include "gaussian_random_numbers.h"
#include "cholesky.h"
#include "matrix_multiplication.h"

using namespace Eigen;


#define GRID_SIZE 900

struct mu_sigma_struct1{
	double mu;
	double sigma;
	vector<float> k_parameters;
};

class Block_inst{
public:
	Block_inst() {};
	int id;
	vector<int> gates_idx;
	float w_slack;
	float prob; //for High vdd
	bool adaptive;
};


class covmat_class{

	friend class TimingAnalysis;

public:
	int x_grids;
	int y_grids;
	int no_of_pc;

	float** MC_matrix;
	MatrixXd cov_matrix; // I will use this matrix for both L & W LOL!
	float *eigen_values;
	float **eigen_vectors;

	void get_eigen_values();

	float** random_numbers;
	float** random_numbers_width;

};



class TimingAnalysis{

	//friend class Block_inst;

public:
	int no_of_pc;
	covmat_class* mat_class;
	void Init_Timing_Analysis(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, igraph_vector_t tmp_nodes, string place_bench, float time_constraint, map<string, Gates* > Gates_Map, \
		map<int, vector<float> >gate_sizes_input, Gates* gates, int no_of_gates, int* edges, int no_of_edges, igraph_t circuitGraph, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool>& adpt_block_exist, int& no_of_adpblk, bool cl = true /*use CL*/);
	void update_gates_delay(Gates* gates, int no_of_gates, int* edges, std::map<std::string, std::vector<float> > gate_sizes);
	void get_req_slack(Gates* gates, int no_of_gates, int* edges, double time_constraint);
	void SSTA(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool> adpt_block_exist);
	void cl_SSTA(Gates2* gates, int no_of_gates, int* edges, int no_of_edges, covmat_class* cov);
	void sum_function(mu_sigma_struct1 max_strc, Gates* currentG);
	mu_sigma_struct1 max_function(Gates* currentG); //return the index of the max delay gate
	float integrate(float);
	float integrate1(float);
	double integrate2(double);
	void update_gates_delay2(Gates2* gates, int no_of_gates, int* edges);
	void update_req_arr_time_and_slack(map<string, pair<double, double> >& arrival_time_slack_local);
	void find_blk_PO();
	void calculate_prob(map<int, vector<float> > &probability);
	void adaptive_SSTA();
	void calculate_yield1(Gates*);
	void calculate_yield1(Gates_cu*);
	void get_max_delay_prob();
	void cu_SSTA();
	void cu_MC();
	//scheduling
	void naive_schedule();
	void new_schedule();

	//
	void backup();
	void backup2();
	void prede_test();

	//Monte Carlo
	void Monte_Carlo_Gen();
	void Monte_Carlo(covmat_class *mat_class, bool adaptive, int blocks);

private:
	//0: max X, 1: max Y, 2: min X, 3: min Y
	float size_of_die[4];
	map<int, vector<float> >gate_sizes;
	Gates* gates_t;
	float req_arr_time;
	int* edges_t;
	int* edges_cu;
	int no_of_gates;
	int no_of_edges;
	int no_of_adpt_blk;
	int count;
	float** cu_k_param;
	float* k_para_matrix;
	float* mu_sigma_matrix;
	vector<Block_inst* > Blocks;
	bool* adpt_blk;
	igraph_vector_t tmp_nodes;  //topological sort
	Gates_cu* gates_cu_p;
	Gates_cu* gates_cu;
	Gates_cu2* gates_cu2;
	Gates2* gates_cl;

	int* tmp_nodes2; //level scheduled
	int max_level;
	int* l_start;
	int* l_count;

	igraph_t circuitGraph;
	int* index_tmp;
	int no_of_iterations;

	

};




























#endif