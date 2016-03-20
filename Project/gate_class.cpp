#ifndef GATE_CLASS_CPP
#define GATE_CLASS_CPP

#include "string"
#include "gate_class.h"
#include "algorithm"
#include "iostream"
#include <math.h>

//#include "gaussian_random_numbers.h"

using namespace std;

gate_class::gate_class(string gate_nm, string input_output_or_gate){
	gate_name  = gate_nm;
	gate_type  = input_output_or_gate;

	gate_type_add  = input_output_or_gate;
    gate_level = 0;
    gate_size  = 1;

    gate_resistance  = 1/gate_size; //Assuming the gate resistance of unit size gate is some 1000 units.
    gate_capacitance = 1*gate_size; //Assuming the gate capacitance of unit size gate is some 1000 units.
}

void gate_class::update_gate_type(gate_class *gate_pointer, string gate_type){
	if(gate_pointer->gate_type.compare("dummy") == 0){
		gate_pointer->gate_type = gate_type;
		gate_pointer->gate_type_add = gate_type;
	}
	if(gate_pointer->gate_type == "PO"){
		gate_pointer->gate_type_add = gate_type;
	}
}

//This Function add the inputs to the gate.
void gate_class::add_inouts(gate_class *gate_pointer){
        vector<gate_class*>::iterator gate_iterator;
    
        //Add the inputs of a gate
        gate_iterator = find(this->inputs.begin(), this->inputs.end(), gate_pointer);
        if(gate_iterator == this->inputs.end()){
		    this->inputs.push_back(gate_pointer);	
        }

        //Add the outputs of a gate
        gate_iterator = find(gate_pointer->outputs.begin(), gate_pointer->outputs.end(), this);
        if(gate_iterator == gate_pointer->outputs.end()){
            gate_pointer->outputs.push_back(this);
        }
}


//This Function updates the Level of the gate in the circuit
void gate_class::update_gate_level(int level){
    this->gate_level = level;
}

//This function recusively assigne the level of the gate
int gate_class::update_gate_level_recursively(gate_class *gate, int n){
	
		static long long t = 0;
		t++;
		//cout << "t = " << t << endl;
		//cout << "Gate Name = " << gate->gate_name << endl;
		//cout << "Gate Type = " << gate->gate_type << endl;
		if(gate->gate_type == "PI"){
			//cout << "A primay input "<< endl;
			//int wait;
			//cin >> wait;
		}
    if("PO" == gate->gate_type && gate->gate_type_add != "DFF"){
        //gate->update_gate_level(n+1);
				//cout << "PO " << gate->gate_name << "   " << n+1 << endl;
        ///return 0;
    } else {
        for(std::vector<gate_class*>::iterator it  = gate->outputs.begin();
                                               it != gate->outputs.end()  ;  ++it){
						if( (*it)->gate_type_add == "DFF" ){
							(*it)->gate_level = 0;
							n = 0;
							//cout << "dff " << (*it)->gate_name << "    " << n << endl;
            	this->update_gate_level_recursively((*it), n);
						} else if((*it)->gate_level < gate->gate_level + 1){
							(*it)->gate_level = gate->gate_level + 1;
							n = (*it)->gate_level;
							//cout << "normal " << (*it)->gate_name << "    " << n << endl;
							//cout << "N = " << n << endl;
            	this->update_gate_level_recursively((*it), n);
						}
          /*  if((*it)->gate_level < n){
                (*it)->update_gate_level(n);
            } */
        }
    return n;
    }
}

void gate_class::update_timings_of_gate_outputs(float delay, int iteration_number, correlation_model *corr, int **adaptive_blocks){
    //Deciding the grid of the current gate and depending upon the grid i am using the random number
    float voltage_level;
    float new_delay;

		float tox = 1.5e-9;
		float mu =  0.01186;
		float L =  45e-9;
		float W =  135e-9;
		float E =  3.9*8.8542e-12;
		float kdL;
		kdL =  2*tox*1e-9/(mu*W*E);
		float kdW;
		kdW =  -2*tox*L*1e-9/(mu*E*W*W);

    float DdL;
    float DdW;
    float p;
    float Vth;
    float Vdd;
    float p_cal;

    Vth = 0.5;

    Vdd= 0.85; //1.2

    p = Vth/Vdd;
    p_cal = (p-0.1)/pow(1-p,2) + log(19-20*p)/(2*(1-p));
    //p_cal = 1;


    DdL = (kdL/Vdd)*p_cal*this->load_capacitance;
    DdW = (kdW/Vdd)*p_cal*this->load_capacitance;

    //cout << "DdL = " << DdL << " p = " << p << endl;
    //cout << "DdW = " << DdW << " p = " << p << endl;
    
    int g;
    if(this->gate_type == "PI"){
        g = 1;
    } else {
        g = 1;
    }



    int   row, column;
    int   random_number;
    row    = floor((this->y_grid + corr->min_y)/GRID_SIZE); //corr->rows);
    column = floor((this->x_grid + corr->min_x)/GRID_SIZE);//corr->columns);
    
    random_number = row*(corr->columns) + column;
    //cout << gate_name << "  randm number " << random_number << endl;


    if(iteration_number == -1){
        voltage_level = 0;
        new_delay = delay + this->delay_of_gate;
    } else {
        voltage_level = adaptive_blocks[this->adaptive_block][iteration_number]%2 == 1 ? 0.5 : 0;
        //voltage_level = 0;
        voltage_level = gate_type == "PI" ? 0 : voltage_level;
    
        //this model is the cell delay model
        // I will assume that the variation in length is 15% and for width it is 10%
        //new_delay = delay + this->delay_of_gate + (corr->random_numbers[random_number][iteration_number])*(0.10)*this->gate_resistance*this->load_capacitance - (corr->random_numbers_width[random_number][iteration_number])*(0.08)*this->gate_resistance*this->load_capacitance;
        //new_delay = delay + this->delay_of_gate + (corr->random_numbers[random_number][iteration_number])*(0.05*45/3)*DdL/100 - (corr->random_numbers_width[random_number][iteration_number])*(0.01*135/3)*DdW/100;
        new_delay = delay + this->delay_of_gate + (corr->random_numbers[random_number][iteration_number])*1*g + (corr->random_numbers_width[random_number][iteration_number])*0*g;
        //new_delay = delay + this->delay_of_gate - 0.2 + (corr->random_numbers[random_number][iteration_number])*(0.17)*this->gate_resistance*this->load_capacitance + (corr->random_numbers_width[random_number][iteration_number])*0.1*this->gate_resistance*this->load_capacitance;
    }

    if(this->output_time_of_gate < new_delay){
        this->output_time_of_gate = new_delay;
    } else if(gate_type != "PI") {
        return;
    }

    if(this->gate_type == "PO"){
        return;
    }

    for(std::vector<gate_class*>::iterator it  = outputs.begin();
                                           it != outputs.end()  ; ++it){
        (*it)->update_timings_of_gate_outputs(this->output_time_of_gate, iteration_number, corr, adaptive_blocks);
    }
}


//This function will assign Gaussian Random Numbers to each gate
void gate_class::assign_gaussian_random_numbers(int no_of_iterations){
    
    //this->Gate_Gaussian_Random_Numbers = gaussian_random_numbers(no_of_iterations); 
}


//This Function returns the value at the output node of the gate
float gate_class::get_output_time(){

    //cout << this->gate_name << endl;   
    return this->output_time_of_gate;
}


void gate_class::print_gates(){
    int i = 0;
    cout << "Printing the Inputs of the gate " << this->gate_name << " " << endl;
    for(std::vector<gate_class*>::iterator it  = this->inputs.begin();
                                           it != this->inputs.end()  ; ++it ){
        cout << this->inputs[i]->gate_name << endl;
        i++;
    }
    
    cout <<"\n\n" ;
    i = 0;

    cout << "Printing the Outputs of the gate " << this->gate_name << " " << endl;
    for(std::vector<gate_class*>::iterator it  = this->outputs.begin();
                                           it != this->outputs.end()  ; ++it ){
        cout << this->outputs[i]->gate_name << endl;
        i++;
    }
    cout <<"\n" ;
    //int wait;
    //cout  << "Waiting... enter an integer " << endl;
    //cin >> wait;
}

void gate_class::update_delay(){

    float load_capacitance = 0;

    if(gate_type.compare("PI") == 0){
        this->delay_of_gate = 0.0;
				this->gate_resistance = 0.0;
        //return;
    }
    
    //This models the elmore delay of the gate
    //Note we are not modeling the wire delay as of now!! Don't know what will be required later
    for(std::vector<gate_class*>::iterator it  = this->outputs.begin();
                                           it != this->outputs.end()  ; ++it){
        load_capacitance += (*it)->gate_capacitance;
    }

    if(this->gate_type == "PO"){
	load_capacitance = 1.0;
    }
    
    this->delay_of_gate = load_capacitance*gate_resistance; 
    //this->delay_of_gate = load_capacitance; 
    this->load_capacitance = load_capacitance;

		//cout << "The ouput delay of the gate  " << this->gate_name << "  is  " << delay_of_gate << "  " << load_capacitance<< endl;
		//int wait;
		//cin >> wait;
    
}

#endif
