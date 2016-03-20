#ifndef GH
#define GH

#include "Celllib.h"
#include "parser_helper.h"


// Gates struct contains all the computation parameters
typedef struct {
	int id;
	COMPONENTTYPE type;  // save the "INV"-like info
	string name;
	/*PI,PO flags*/
	bool is_input;
	bool is_output;
	/*Arguments to find neighbors in the adjlist*/
	int start_in;
	int no_of_in;
	int start_out;
	int no_of_out;
	/*level*/
	int level;


	//timing analysis parameters
	float gate_resistance;
	float gate_capacitance;
	float load_capacitance;

	//position
	int x;
	int y;

	float output_time_of_gate;
	float delay;
	float slack;
	float req_arrive_time;
	float gate_mu;
	float gate_sigma;

	float delay_mu;
	float delay_sigma;
	//pair<double, double> gate_mu_sigma;
	//vector<double> k_paras;		// use seprate data to store K values
	float* k_param;
	int no_of_k;
	int adapt_blocks_id;
	int nn;
	int in_degree_v;
	float k_v;
	int ready_cnt;
	int label;
	bool finish;

	//test
	int selfV;
	int getV;

}Gates;

typedef struct {
	int id;
	/*PI,PO flags*/
	bool is_input;
	bool is_output;
	/*Arguments to find neighbors in the adjlist*/
	int start_in;
	int no_of_in;
	int start_out;
	int no_of_out;
	/*level*/
	int level;

	int x;
	int y;
	float delay;
	float gate_mu;
	float gate_sigma;

	float output_time_of_gate;
	float delay_mu;
	float delay_sigma;
	//float* k_param;
	int adapt_blocks_id;
	

}Gates_cu;

typedef struct {

	bool is_input;
	/*Arguments to find neighbors in the adjlist*/
	int start_in;
	int no_of_in;
	int x;
	int y;
	float delay;
	float output_time_of_gate;

}Gates_cu2;
#endif