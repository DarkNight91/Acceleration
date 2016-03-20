/* project body*/
#ifndef OCL_OPT_H
#define OCL_OPT_H

#include <igraph.h>
#include "Gates.h"
#include <string>
#include <fstream>
#include <cstdio>
#include <vector>
#include <map>

class TimingAnalysis;

//#include "CLHelper.h"
using namespace std;

/*-----Reduced size of Gates for GPU computing, the parameter to find its fanin and fanout may no longer use in the future, cuz the ssta computing takes majority of the time-----*/
typedef struct{

	int id;
	COMPONENTTYPE type;  // save the "INV"-like info
	const char* name;
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

	//int adapt_blocks_id;

	//timing analysis parameters
	float gate_resistance;
	float gate_capacitance;
	float load_capacitance;

	//position
	int x;
	int y;

	float delay;

	float req_arrive_time;
	float gate_mu;
	float gate_sigma;
	float delay_mu;
	float delay_sigma;

	//float* k_param;
	int no_of_k;
	int adapt_blocks_id;
	int nn;
	
	}Gates2;



class ocl_opt{

	public:
		Gates* gates;
		
		map<int, string > Gate_Id_Name_Map;

		//map name and the gate, value is a pointer
		map<string, Gates*> Gates_Map;

		int* adj_edges;
		igraph_vector_t edges;
		igraph_t circuitGraph;
		typedef std::map<std::string, COMPONENTTYPE>  GMAP_INDEX;
		GMAP_INDEX gGate_Index;




		bool init_lracc_database();
		void print_edges();
		int find_wire(string& name, vector<wire_par* >& wirePar);
		bool createEdges (map<string, wire_par* >& wirePar, Gates* gates);
		bool  load_graph_and_tiles(string  &filename);
		void adjlist_generate(Gates* gates, int no_of_gates);
		void adjlist_regenerate(Gates* gates, int no_of_gates, int* index);
		int gates_levelize(Gates* gates,int no_of_gates/*,int* adj_edges*/);
		int find_level (Gates gate, Gates* gates/*, int* adj_edges*/);
		void QuickSortbyLevel(Gates* gates, int no_of_gates);
		
		void load_at_raq();
		void update_gates_info();

		//SSTA
		void init_timing(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, string place_bench, double time_constraint, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool>& adpt_block_exist, int& no_of_adpblk, bool use_coff);
		bool load_celllib(char *filename);
		void calculate_prob(map<int, vector<float> >&   block_prob_local);
		void SSTA(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool> adpt_block_exist);
		void pass_pca_plac_info(map<string, vector<double> >& pca_param_local, map<string, vector<float> >& plac_info_local);
		

private:
	int no_of_gates;
	int no_of_edges;
	CELLLIB cell_lib;
	map<int, vector<float> > gate_sizes;
	igraph_vector_t tmp_nodes;
	TimingAnalysis* t;
};

#endif