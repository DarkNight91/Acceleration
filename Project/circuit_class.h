#ifndef CIRCUIT_CLASS_H
#define CIRCUIT_CLASS_H

#include <string>
#include <map>

#include "define.h"
#include "gate_class.h"
#include "correlation_model.h"
#include "data_structures.h"

using namespace std;


struct string_compare{
    public:
        bool operator()(const string str1, const string str2){
            return str1.compare(str2) < 0;
        }
};

class circuit_class{
   
    friend class correlation_model;

	public:
		//map<string, gate_class*, string_compare > Gates_Map;
		map<std::string, gate_class*> Gates_Map;

		map<int, vector<gate_class*> > levels;

    // This will make the circuit work on either low or high depending upon the value
    // true = high
    // false = low
    bool *high_low;
		int no_of_adaptive_blocks;

    double **ran;

        //This maps contains all the adaptive boundary crossing nets in which the first key is
        //an 'int' which tells the adaptive block number and the 'string' is net or output of the 
        //gate which lies in this adaptive block but have an output in the other net
        map<int, vector<gate_class*> > adp_bnd_crossing_nets;

        //the following map contains the probability of the the adaptive block going high or low!!
        map<int, double> prob_high;
        map<int, double> prob_low;

        map<int, double> high;
        map<int, double> low;
        
        // These vector are helpful when we want to Levelize the circuit and during the pert algorithm
        vector<gate_class*> Primary_Outputs;
        vector<gate_class*> Primary_Inputs;

        map<pair<int,int>, int > granularity;    // 0 = Low Vdd, Low Body Bias.   1 = High Vdd, Low Body Bias.
                                                     // 2 = Low Vdd, High Body Bias.  3 = High Vdd, High Body Bias

				bool *adapt_block_exist;

	public:
        double **random_numbers;
		void add_gate_to_map(string gate_nm, string input_output_or_gate);
		void add_gates_input_output(string gate_nm, string input_output_name);
        void Levelize_Circuit(circuit_class *ckt);
        void update_gate_delays(map<string, vector<float> >  gate_parameters, double required_output_time);
        void gen_gaussian_random_numbers(int no_of_iterations);
        void pert_for_timing(int iteration_number, correlation_model *corr, int **adaptive_blocks);
        double  print_time_at_outputs();
        void determine_boundary_crossing_nets();
        void calculate_probability_low_high(int no_of_adaptive_blocks, bool only_FBB);
        gate_class* trace_input_of_adaptive_block(gate_class* gate);
        void calculate_probability_RBB(int no_of_adaptive_blocks, bool only_FBB);

        void calculate_required_arrival_time(double required_time_at_output);
        void update_required_arrival_time(gate_class* gate);
        double get_max_delay_prob(double threshold);
				double Critical_Path_Delay(correlation_model *corr, int x);

        void Monte_Carlo1(correlation_model *corr, int no_of_iterations, bool adaptive, int blocks, double threshold);
        double Run_timing(gate_class* gate, correlation_model *corr, int iteration_number, bool adaptive, int x );

        pair<double, double> Calculate_Power_Monte_Carlo(int no_of_iterations, correlation_model *corr, map<string, vector<float> > gate_parameters);
        double get_sigma(double *power_c, double mean, int no_of_iterations);

        //Functions for PCA
        void   SSTA(correlation_model *corr, int no_of_prin_comp);
        double calculate_sigma_of_gate(gate_class* gate, correlation_model *corr, int no_of_prin_comp, double adapt);
        void   calculate_mu_and_sigma(gate_class* gate, int no_of_prin_comp);
        void   SSTA_Addition_Function(gate_class* gate, mu_sigma_struct s1, int no_of_prin_comp);
        mu_sigma_struct SSTA_Max_Function(mu_sigma_struct gate_1, mu_sigma_struct gate_2, int no_of_prin_comp);

        //Functions to calculate the Required Arrival Time
        void SSTA_Minus_Function(gate_class* gate, mu_sigma_struct s1, int no_of_prin_comp);
        mu_sigma_struct SSTA_Minimum_Function(mu_sigma_struct gate_1, mu_sigma_struct gate_2, int no_of_prin_comp);
        void Calculate_RAT(int no_of_prin_comp, double required_output_time);
        void calculate_rat_recursively(gate_class *gate, int no_of_prin_comp, double required_output_time);
        void Prune_nodes();
        void prune_nodes_recursively(gate_class* gate);
        void prunes_gate_in_line(gate_class* gate);

        void Monte_Carlo_RAT(correlation_model *corr, int no_of_iterations, bool adaptive, double required_output_time);
        double Run_timing_rat(gate_class* gate, correlation_model *corr, int iteration_number, double required_output_time);

        void compress_gates();
        void compress_gates_recursively(gate_class* gate, vector<gate_class*> &del_gates);
        
};

#endif
