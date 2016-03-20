#ifndef GATE_CLASS_H
#define GATE_CLASS_H

#include "define.h"
#include "string"
#include "vector"
#include "correlation_model.h"

using namespace std;

class gate_class{

    friend class circuit_class; 
    friend class correlation_model;   

	public:
		vector<gate_class*> inputs;
        vector<gate_class*> outputs;

        float* Gate_Gaussian_Random_Numbers;
		    string gate_name;	
		    string gate_type;
		    string gate_type_add;
		    double delay_of_gate;
        //double arrival_time;
        double required_arrival_time;

        float  gate_size;          //This parameter is an input to the timing analysis , currently set to 1
        float  gate_resistance;
        float  gate_capacitance;
        float  load_capacitance;
        float  output_time_of_gate;
        int    gate_level;

        double leakage_power;

        float Vth;

        pair<double, double> mu_sigma;
        pair<double, double> output_mu_sigma;

        pair<double, double> rat_mu_sigma;
        vector<double> rat_k_parameters;
        pair<double, double> rat_input_mu_sigma;
        vector<double> rat_input_k_parameters;

        double rat;
        double rat_input;

        vector<double> k_parameters;
        vector<double> output_k_parameters;
       
        float  x_grid;
        float  y_grid;

        int adaptive_block;

	public:
		gate_class() {};
		gate_class(string gate_nm, string input_output_or_gate);
		void   update_gate_type(gate_class *gate_pointer, string gate_type);
		void   add_inouts(gate_class *gate_pointer);
        void   update_gate_level(int level);
        int    update_gate_level_recursively(gate_class *gate, int n); 
        int    get_gate_level();
        float  get_output_time();
        string get_gate_type();
        void   update_timings_of_gate_outputs(float delay, int iteration_number, correlation_model *corr, int **adaptive_blocks);
        void   assign_gaussian_random_numbers(int no_of_iterations);
        void   print_gates();
        void   update_delay();
};

#endif
