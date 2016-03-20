#ifndef CIRCUIT_CLASS_CPP
#define CIRCUIT_CLASS_CPP

#include <string>
#include <map>
#include <iostream>
#include <stdio.h>
#include "string"
#include <math.h>
#include "integrate.h"
#include "gaussian_random_numbers.h"

#include "circuit_class.h"

double* ran_1;
double* ran_2;

using namespace std;

//This Function add the gates instances to the map
void circuit_class::add_gate_to_map(string gate_nm, string input_output_or_gate){
	if(this->Gates_Map.find(gate_nm) == this->Gates_Map.end()){
		gate_class *gate;

		gate = new gate_class(gate_nm, input_output_or_gate);
		this->Gates_Map[gate_nm] = gate;

        if(input_output_or_gate == "PO"){ // I will also push the DFF gate as a primary output too !!
            this->Primary_Outputs.push_back(gate);
        }

  			if(input_output_or_gate == "DFF_I"){
						if(find(Primary_Outputs.begin(), Primary_Outputs.end(), Gates_Map[gate_nm]) == Primary_Outputs.end()){
    						this->Primary_Outputs.push_back(this->Gates_Map.find(gate_nm)->second);
						}
				}

        if(input_output_or_gate == "PI" || input_output_or_gate == "DFF"){ // Here we have marked the DFF also as a primary input too !! this is done from the parser file !!
                                          // ^ what i have written in the comments needs to be checked properly !! -> TODO
            this->Primary_Inputs.push_back(gate);
        }

	} else if(input_output_or_gate.compare("dummy") != 0 && input_output_or_gate != "DFF_I" ){
		Gates_Map[gate_nm]->update_gate_type(Gates_Map[gate_nm], input_output_or_gate);
  	if(input_output_or_gate == "DFF"){
    	this->Primary_Inputs.push_back(this->Gates_Map.find(gate_nm)->second);
    }
	} else {
  	if(input_output_or_gate == "DFF_I"){
			if(find(Primary_Outputs.begin(), Primary_Outputs.end(), Gates_Map[gate_nm]) == Primary_Outputs.end()){
    		this->Primary_Outputs.push_back(this->Gates_Map.find(gate_nm)->second);
			}
    }
	}
}

//This Function add the inputs/outputs to a vector List of a Gate
void circuit_class::add_gates_input_output(string gate_nm, string input_output_name){
	this->Gates_Map[gate_nm]->add_inouts(Gates_Map[input_output_name]);
}


// This function Levelizes the circuit
void circuit_class::Levelize_Circuit(circuit_class *ckt){
    int n;

		for(map<string,gate_class*>::iterator it = Gates_Map.begin(); it != Gates_Map.end(); ++it){
			//levels[it->second->gate_level].push_back(it->second);
			it->second->gate_level = 0;
		}

    for(std::vector<gate_class*>::iterator it  = ckt->Primary_Inputs.begin();
                                           it != ckt->Primary_Inputs.end()  ;  ++it){
        n = (*it)->update_gate_level_recursively(*it, 0);
        //(*it)->update_gate_level(n);
    }

		for(map<string,gate_class*>::iterator it = Gates_Map.begin(); it != Gates_Map.end(); ++it){
			//cout << it->second->gate_name << "   " << it->second->gate_level << endl;
			levels[it->second->gate_level].push_back(it->second);
		}
}


//This Circuit performs the Normal Timing Analysis
void circuit_class::pert_for_timing(int iteration_number, correlation_model *corr, int **adaptive_blocks){

    //this->random_numbers = random_numbers;

    for(std::map<string, gate_class*,string_compare>::iterator it  = this->Gates_Map.begin();
                                                               it != this->Gates_Map.end()  ;  ++it){
        (it->second)->output_time_of_gate = 0.0;
    }

    for(std::vector<gate_class*>::iterator it  =  this->Primary_Inputs.begin();
                                           it !=  this->Primary_Inputs.end()  ; ++it){
        /*cout << "The inputs of the circuit are " << (*it)->gate_name  << endl;
        int wait;
        cin >> wait;*/
        (*it)->update_timings_of_gate_outputs(0.0, iteration_number, corr, adaptive_blocks);
    }

}

//This function will assign Gaussian Random Numbers to each gate
void circuit_class::gen_gaussian_random_numbers(int no_of_iterations){
    
    for(std::map<string, gate_class*,string_compare>::iterator it  = this->Gates_Map.begin();
                                                it != this->Gates_Map.end()  ;  ++it){
        (it->second)->assign_gaussian_random_numbers(no_of_iterations); 
    }
    
}

void circuit_class::calculate_required_arrival_time(double required_time_at_output){

		//cout << "In required arrival time calculations " << endl;
    for(vector<gate_class*>::iterator it  = Primary_Outputs.begin();
                                      it != Primary_Outputs.end()  ; ++it){
        update_required_arrival_time(*it);
    }

		//cout << "FINISHED WITH ALL THE REQUIRED ARRIVAL TIME CALCULATIONS !@!@ " << endl;
}

void circuit_class::update_required_arrival_time(gate_class* gate){
		double delay;
		double adapt;
    for(vector<gate_class*>::iterator it_i  = gate->inputs.begin();
                                      it_i != gate->inputs.end()  ; ++it_i ){

					if(this->high_low[gate->adaptive_block] == true){
						adapt = 1.1;
					} else {
						adapt = 1;
					}
					//cout << "here in update required arrival time "	 << gate->gate_name << "  " << (*it_i)->gate_name << endl;
					//cout << Gates_Map["n_131"]->gate_name << endl;
					 delay = gate->required_arrival_time - gate->mu_sigma.first/adapt - 3*gate->mu_sigma.second/adapt;
           if((*it_i)->required_arrival_time > delay){
              (*it_i)->required_arrival_time = delay;//gate->required_arrival_time - gate->delay_of_gate;
							if((*it_i)->gate_type_add != "DFF"){
              	update_required_arrival_time(*it_i);
							}
           }
    }
	//cout << "FINISHED UPDATING REQUIRED ARRIVAL TIME " << endl;
}


//This function will print the time at the primary outputs of the circuit
double circuit_class::print_time_at_outputs(){

    gate_class *temp_pointer;
    float delay = 0;
		double mu, sigma;
		//FILE *fp2;
		//fp2 = fopen("musigma.txt","w");

    for(std::map<string, gate_class*,string_compare>::iterator it  = this->Gates_Map.begin();
                                                it != this->Gates_Map.end()  ;  ++it){
        //if((it->first).compare("G17") == 0){
        //if((it->first).compare("G550") == 0){
				if(delay < (it->second)->output_time_of_gate){
					delay = (it->second)->output_time_of_gate;
					temp_pointer = it->second;
				}
    } 

  	//temp_pointer = Gates_Map["G9"];
  	//temp_pointer = Gates_Map["N6288"];
  	//temp_pointer = Gates_Map["N11334"];

		//cout << endl << endl << "The gate with the larget delay is " << temp_pointer->gate_name << "    " << delay << endl;
		//cout << "The mu is " << temp_pointer->output_mu_sigma.first << endl;
		//cout << "The sigma is " << temp_pointer->output_mu_sigma.second << endl;

    //for(std::vector<gate_class*>::iterator it  = Primary_Outputs.begin();
    //                                       it != Primary_Outputs.end() ;  ++it){
        //delay = temp_pointer->get_output_time();
        //mu    = temp_pointer->output_mu_sigma.first;
        //sigma = temp_pointer->output_mu_sigma.second;

      //  if(*it == temp_pointer){
            //fprintf(fp,"%f  %f  %f\n",delay, mu, sigma);
            //fprintf(fp2,"%f  %f  %f\n",delay, mu, sigma);            
            //fprintf(fp,"%f\n",delay);            
						return delay;
        //}
    //}

		//fclose(fp2);
    
    return 0;

}

void circuit_class::update_gate_delays(map<string, vector<float> >  gate_parameters, double required_output_time){
    
    for(std::map<std::string, gate_class*,string_compare>::iterator it  = Gates_Map.begin();
                                                                    it != Gates_Map.end()  ;  ++it){
		/*
		(it->second)->gate_resistance = (gate_parameters[(it->second)->gate_name])[1];
        (it->second)->gate_capacitance = (gate_parameters[(it->second)->gate_name])[2];
        (it->second)->leakage_power = gate_parameters[(it->second)->gate_name][4];
				if(it->second->gate_type == "PI"){
					it->second->gate_resistance = 0.0;
					it->second->gate_capacitance = 0.0;
				}
				*/
		if (it->second->gate_type == "PI"){
			it->second->gate_resistance = 0.0;
			it->second->gate_capacitance = 0.0;
		}
		else{
			(it->second)->gate_resistance = (gate_parameters[(it->second)->gate_name])[1];
			(it->second)->gate_capacitance = (gate_parameters[(it->second)->gate_name])[2];
			(it->second)->leakage_power = gate_parameters[(it->second)->gate_name][4];
		}
				//cout << it->second->gate_name <<"  "  <<it->second->gate_resistance << "  " << it->second->gate_capacitance << endl ;

				if(abs(it->second->gate_resistance) > 10000 ||  abs(it->second->gate_capacitance) > 10000 ) {
					cout << "HEY JIAFAN.. THERE IS A PROBLEM FROM THE DATA RECEIVED ... " << endl;
					cout << "THE RESISTANCE AND CAPACITANCE FOR THE GATE " << it->second->gate_name << " ARE " << it->second->gate_resistance << "   " << it->second->gate_capacitance << "RESPECTIVELY "<< endl;
					//int wait;
					//cin >> wait;
				}
        (it->second)->required_arrival_time = required_output_time;
        //(it->second)->threshold_voltage = required_output_time;
    }

    for(std::map<std::string, gate_class*,string_compare>::iterator it  = Gates_Map.begin();
                                                                    it != Gates_Map.end()  ;  ++it){
        (it->second)->update_delay();
		}
    //calculate_required_arrival_time(required_output_time);
}

void circuit_class::determine_boundary_crossing_nets(){
    for(std::map<std::string, gate_class*>::iterator it  = Gates_Map.begin();
                                                                    it != Gates_Map.end()  ;  ++it){
	if(it->second->gate_type == "PO"){
	    this->adp_bnd_crossing_nets[(it->second)->adaptive_block].push_back( (it->second) );
	}
  for(vector<gate_class*>::iterator it_out  = (it->second)->outputs.begin();
                                    it_out != (it->second)->outputs.end()  ; ++it_out){

		if( (*it_out)->adaptive_block != (it->second)->adaptive_block ){
    	this->adp_bnd_crossing_nets[(it->second)->adaptive_block].push_back( (it->second) );
      	break;
      }
    }
	}
}


void circuit_class::calculate_probability_low_high(int no_of_adaptive_blocks, bool only_FBB){
    gate_class* max_input_gate;
    pair<double,double> prob_mu_sigma;
    double max_prob = 0.0; // This is the probability of a adaptive block going HIGH !!
    double prob;
    double slack;
    double threshold;
		double delay;
    double max_slack;
		double max_delay;

	  //cout << "IN calculate probability function " << endl;

		//cout << "adaptive blocks = " << no_of_adaptive_blocks << endl;

    for(int i=0; i < no_of_adaptive_blocks; i++){
        max_prob = 0.0;
				prob_high[i] = 0.0;
        max_slack = 0;
        max_delay = 0;
        for(vector<gate_class*>::iterator it  = adp_bnd_crossing_nets[i].begin();
                                          it != adp_bnd_crossing_nets[i].end()  ; ++it){

            // The following if statement will be true in case the gate fails the timing
            if( (*it)->required_arrival_time < (*it)->output_mu_sigma.first + 3*((*it)->output_mu_sigma.second) ){
                // if the timing has failed then we have to trace back the input from where the 
                // delay has came and it should be of the same adaptive block
                max_input_gate = trace_input_of_adaptive_block( *it );
                prob_mu_sigma.first  = (*it)->output_mu_sigma.first  - max_input_gate->output_mu_sigma.first;
                prob_mu_sigma.second = pow((*it)->output_mu_sigma.second, (double)2) + pow(max_input_gate->output_mu_sigma.second, (double)2);
		
      if(prob_mu_sigma.second != 0){
                slack = ((*it)->required_arrival_time) - (*it)->output_mu_sigma.first - 3*(*it)->output_mu_sigma.second;
                threshold = prob_mu_sigma.first + 3*sqrt(prob_mu_sigma.second) + abs(slack);
                delay = prob_mu_sigma.first + 3*sqrt(prob_mu_sigma.second);
                //threshold = prob_mu_sigma.first + 3*sqrt(prob_mu_sigma.second) + slack;
                //prob = integrate( threshold, prob_mu_sigma.first, sqrt(prob_mu_sigma.second) );
								if(slack < 0){
									prob = abs(slack)/(threshold);
								} else {
									prob = 0;
								}
			} else {
				prob = 0;
			}

                if(prob > 1)
                    prob = 1.0;
            
                if(max_prob <  prob){
                    max_prob =  prob;
                    max_slack = slack;
                    max_delay = delay;
                }
            }
            //cout << "Max Prob = " << i << "    " << max_prob << endl;
						if(abs(max_slack) > max_delay - max_delay/1.1){
                if(0){
        			      prob_low[i] = max_prob*((max_delay - max_delay/1.1)/(max_delay - max_delay/1.2));
        			      prob_high[i] = max_prob*((max_delay/1.1 - max_delay/1.2)/(max_delay - max_delay/1.2)); 
							      high[i] = 1;
							      low[i] = 1;
                }  else {
        			      prob_low[i] = max_prob;
                    prob_high[i] = 0;
							      high[i] = 0;
							      low[i] = 1;
                }
						} else {
        			prob_low[i] = max_prob;	
        			prob_high[i] = 0;	
							low[i] = 1;
							high[i] = 0;
        			//prob_high[i] = max_prob*((delay/1.1 - delay/1.2)/(delay - delay/1.2));
						}
			}
    }
		//cout << "finished the FBB calc " << endl;
    //calculate_probability_RBB(no_of_adaptive_blocks, only_FBB);
}

void circuit_class::calculate_probability_RBB(int no_of_adaptive_blocks, bool only_FBB){
    gate_class* max_input_gate;
    pair<double,double> prob_mu_sigma;
    double min_prob = 0.0; // This is the probability of a adaptive block going HIGH !!
    double prob;
    double slack;
    double threshold;

	  //cout << "IN calculate probability RBB function " << only_FBB <<endl; 

		if(true == only_FBB) {
			for(int i=0; i < no_of_adaptive_blocks; i++){
				prob_low[i] = 0;
			}
			return;
		}

    for(int i=0; i < no_of_adaptive_blocks; i++){
        min_prob = 0.0;
				prob_low[i] = 0.0;
        for(vector<gate_class*>::iterator it  = adp_bnd_crossing_nets[i].begin();
                                          it != adp_bnd_crossing_nets[i].end()  ; ++it){
            // The following if statement will be true in case the gate fails the timing
            if( (*it)->required_arrival_time > (*it)->output_mu_sigma.first + 3*((*it)->output_mu_sigma.second) && prob_high[i] <= 0 ){
                // if the timing has failed then we have to trace back the input from where the 
                // delay has came and it should be of the same adaptive block
                max_input_gate = trace_input_of_adaptive_block( *it );
                prob_mu_sigma.first  = (*it)->output_mu_sigma.first  - max_input_gate->output_mu_sigma.first;
                prob_mu_sigma.second = pow((*it)->output_mu_sigma.second, (double)2) + pow(max_input_gate->output_mu_sigma.second, (double)2);
		
                slack = ((*it)->required_arrival_time) - (*it)->output_mu_sigma.first - 3*(*it)->output_mu_sigma.second;
                if(slack < 0){
                    cout << "there is an error in probability calculation !!" << endl;
                    slack = 0;
                }
                threshold = prob_mu_sigma.first + 3*prob_mu_sigma.second + slack;
                prob = slack/threshold;
								//cout << "LOW PROB = " << prob << endl;
								//int wait;
								//cin >> wait;
                //prob = integrate( threshold, prob_mu_sigma.first, sqrt(prob_mu_sigma.second) );

                if(prob > 1)
                    prob = 1.0;
            
                if(min_prob < prob){
                    min_prob = prob;
                }
            }
        prob_low[i] = min_prob;
			}
    }
 }

// This function will return the input time of the adaptive block which lies on the critical path to this output!!
gate_class* circuit_class::trace_input_of_adaptive_block(gate_class* gate){
    pair<double, double> temp_mu_sigma;
    gate_class* temp_gate;
    temp_mu_sigma = make_pair(0.0,0.0);


		//cout << "IN TRACE ADAPTIVE BLOCK " << endl;
    
		if(gate->gate_type_add == "DFF"){
			return gate;
		}

		temp_gate = gate;
    for(vector<gate_class*>::iterator it  = gate->inputs.begin();
                                      it != gate->inputs.end()  ; ++it){
		//cout<< "The name of the net is " << (*it)->gate_name << " and type  is " << (*it)->gate_type << endl;
	// Here I have to select the largest input
        if( temp_mu_sigma.first + 3*temp_mu_sigma.second <= 
						(*it)->output_mu_sigma.first + 3*((*it)->output_mu_sigma.second)){
            temp_mu_sigma = (*it)->output_mu_sigma;
            temp_gate     = (*it);
        }
    }
 
   // cout << endl << endl;
   // cout << "In the Trace input Function" << endl;
   // cout << "The gate name is " << gate->gate_name << endl;
   // cout << "The input gate type is " <<  temp_gate->gate_name << " " <<temp_gate->gate_type << endl << endl;
    
    if(temp_gate->adaptive_block != gate->adaptive_block || ((temp_gate->gate_type).compare("PI") == 0) ){ // This means we have reached an input of a adaptive block
                                                           //  And we have found the input on this critical path
			 //cout << "Returning the input gate from trace back the input " << temp_gate->gate_name << endl;
			 //cout << "Returning the input gate from trace back the input " << temp_gate << endl;
			 //cout << "Returning the input gate from trace back the input " << temp_gate->adaptive_block << endl << endl;
       return temp_gate; 
    } else {
        return trace_input_of_adaptive_block(temp_gate);
    }
}



// Writing these new function for a different implementation of Monte Carlo

 void circuit_class::Monte_Carlo1(correlation_model *corr, int no_of_iterations, bool adaptive, int blocks, double threshold){

    FILE *fp;
    gate_class* temp_gate = NULL;
    fp = fopen("delay_monte.txt", "w");
    double prob;
    float bin_low = 0;
    float bin_high = 0;

    if(adaptive == true){
        ran = new double*[blocks];
        for(int y=0; y<blocks; y++){
            ran[y] = new double[no_of_iterations];
            srand(y*8);
            for(int g=0; g<no_of_iterations; g++){
                ran[y][g] = (float)rand()/RAND_MAX;
            }
        }
    }

    for(int i=0; i<no_of_iterations; i++){
        for(std::map<string, gate_class*>::iterator it1  = Gates_Map.begin();
                                                    it1 != Gates_Map.end()  ; it1++){
            (it1->second)->output_time_of_gate = 0.0;
        }

        for(vector<gate_class*>::iterator it  = Primary_Outputs.begin();
                                          it != Primary_Outputs.end()  ; ++it){
            Run_timing(*it, corr, i, adaptive, 0);
        }

        prob = 1; 
        if(adaptive ==  true){
            for(int y=0; y<blocks; y++){
                if( ran[y][i] <= prob_low[y] ){
                    prob *= prob_low[y];
                } else {
                    prob *= (1-prob_low[y]);
                }
            }
    /*
            if( ran[1][i] <= prob_high[1] ){
                prob *= prob_high[1];
            } else {
                prob *= (1-prob_high[1]);
            }
            if( ran[2][i] <= prob_high[2] ){
                prob *= prob_high[2];
            } else {
                prob *= (1-prob_high[2]);
            }
            if( ran[3][i] <= prob_high[3] ){
                prob *= prob_high[3];
            } else {
                prob *= (1-prob_high[3]);
            }*/
        }

        double delay;

        delay = 0.0;

        for(std::map<string, gate_class*,string_compare>::iterator it  = this->Gates_Map.begin();
                                                it != this->Gates_Map.end()  ;  ++it){
        //if((it->first).compare("G17") == 0){
        //if((it->first).compare("G550") == 0){

				if(delay < (it->second)->output_time_of_gate){
					delay = (it->second)->output_time_of_gate;
					temp_gate = it->second;
				}
    } 


        //temp_gate = Gates_Map["n_103688"];
        //temp_gate = Gates_Map["n_4000"];
        //temp_gate = Gates_Map["N6288"];
        //cout << i <<" The Delay for the gate " << temp_gate->gate_name << "   " << delay <<endl;
        //cout << "The mu for gate " << temp_gate->gate_name << "   " << temp_gate->output_mu_sigma.first <<endl;
        //cout << "The sigma for gate " << temp_gate->gate_name << "   " << temp_gate->output_mu_sigma.second <<endl;
        if(temp_gate->output_time_of_gate < threshold){
            bin_low++;// = bin_low + 1*prob;
        } else {
            bin_high++;
        }

        if(delay >= 1200){
            //cout << "HELLO" << endl;
        }

        //fprintf(fp,"%f\n",temp_gate->output_time_of_gate);
    }
    //cout << "The Timing Yield for Monte Carlo = " << bin_low/10000;
    cout << "The Timing Yield for Monte Carlo = " << bin_low/no_of_iterations << endl;
    fclose(fp);
 }

 double circuit_class::Critical_Path_Delay(correlation_model *corr, int x){
    double delay  = 0;
    double delay_max = 0;
		gate_class* gate = NULL;

    for(std::map<string, gate_class*>::iterator it1  = Gates_Map.begin();
                                                it1 != Gates_Map.end()  ; it1++){
        (it1->second)->output_time_of_gate = 0.0;
    }

    for(vector<gate_class*>::iterator it  = Primary_Outputs.begin();
                                      it != Primary_Outputs.end()  ; ++it){
        delay = Run_timing(*it, corr, -1, false, x);
        if(delay > delay_max){
            delay_max = delay;
						gate = *it;
        }
   }

    if(x == 1){
        //cout << gate->gate_name << endl;
		    //cout << "The BEST case mu and sigma of the output is " << gate->output_mu_sigma.first << "  " << gate->output_mu_sigma.second << endl;
    } else {
        //cout << gate->gate_name << endl;
		    cout << "The WORST case mu and sigma of the output is " << gate->output_mu_sigma.first << "  " << gate->output_mu_sigma.second << endl;
    }
    return delay_max;
 }

 double circuit_class::Run_timing(gate_class* gate, correlation_model *corr, int iteration_number, bool adaptive, int x){
    int   row, column;
    int   random_number;
    row    = floor((gate->y_grid + corr->min_y)/GRID_SIZE); //corr->rows);
    column = floor((gate->x_grid + corr->min_x)/GRID_SIZE);    //corr->columns);
    random_number = row*(corr->columns) + column;

/*		if(gate->gate_type_add == "DFF"){
			return 
		}*/

    double temp_max = -10000;
    if(gate->gate_type == "PI"){
    //    return corr->random_numbers[random_number][iteration_number];
        temp_max = 0;
    }

		//cout << "inside the run timing fucntion !@!@ " << endl;
    //double temp_max1 = 0;

		if(gate->gate_type != "DFF") {
    for(vector<gate_class*>::iterator it  = gate->inputs.begin();
                                      it != gate->inputs.end()  ; ++it){
			 if((*it)->gate_type_add == "DFF"){
					(*it)->output_time_of_gate = (*it)->delay_of_gate;
			 } else if((*it)->output_time_of_gate == 0.0){
           Run_timing(*it, corr, iteration_number, adaptive, x);
       }

       if( (*it)->output_time_of_gate > temp_max){
           temp_max = (*it)->output_time_of_gate;
       }
       
    }
		}

    double adapt;
    if(adaptive == true){
        if(ran[gate->adaptive_block][iteration_number] <= prob_low[gate->adaptive_block] ){
            adapt = 1.2;
        }else{
            adapt = 1;
        }
    } else{
        adapt = 1;
    }

    if( x == 1 ){
        adapt = 1.2;
    }
    //cout << "ADAPT = " << adapt << endl;
    if(iteration_number == -1){
        gate->output_time_of_gate =  temp_max + (gate->mu_sigma.first + 3*(gate->mu_sigma.second))/adapt ;
        //cout << "DELAY = " << (gate->mu_sigma.first + 3*(gate->mu_sigma.second))/adapt << "   " << temp_max << endl;
    } else {
        gate->output_time_of_gate = temp_max + gate->delay_of_gate/adapt + gate->delay_of_gate*(0.15/3)*(corr->random_numbers[random_number][iteration_number])/adapt - 1*(corr->random_numbers_width[random_number][iteration_number])*gate->delay_of_gate*(0.08/3)/adapt;
    }
	return gate->output_time_of_gate;

 }


 /********************************************************************************/
 // This Function calclates the power of the circuit 
 /********************************************************************************/

 pair<double, double> circuit_class::Calculate_Power_Monte_Carlo(int no_of_iterations, correlation_model *corr, map<string, vector<float> > gate_parameters ){
     double *power_c;
     //double power_sigma;
     double mean = 0.0;
     double adapt;

		for(int j=0; j < this->no_of_adaptive_blocks; j++){
			cout << j <<" " << prob_high[j] << "   " << 1 - prob_high[j] << endl;
		}
		//int wait;
		//cin >> wait;


    int   row, column;
    int   random_number;

     pair<double, double> power_mu_sigma;
		 //power_mu_sigma = make_pair<0.0,0.0>;

     power_c = new double[no_of_iterations];

     for(int i=0; i<no_of_iterations; i++){

         power_c[i] = 0.0;
        for(std::map<string, gate_class*>::iterator it  = this->Gates_Map.begin();
                                                     it != this->Gates_Map.end()  ; ++it){
            if(this->ran[it->second->adaptive_block][i] <= prob_high[it->second->adaptive_block] && this->high_low[it->second->adaptive_block] == true ){
                adapt = -0.05;
            }else if (this->ran[it->second->adaptive_block][i] <= prob_low[it->second->adaptive_block] && this->high_low[it->second->adaptive_block] == true) {
                adapt = 0.05;
            } else {
                adapt = 0;
            }

						it->second->leakage_power = gate_parameters[it->second->gate_name][4];

						if(it->second->leakage_power <= 0 && it->second->gate_type != "PI" && it->second->gate_type != "PO"){
							cout << "The leakage power received from primal code is zero or less than 0" << it->second->leakage_power << endl;
						}

						if(it->second->gate_type == "PI"){
							it->second->leakage_power = 0;
						}
            
            row    = floor((it->second->y_grid + corr->min_y)/GRID_SIZE); //corr->rows);
            column = floor((it->second->x_grid + corr->min_x)/GRID_SIZE);    //corr->columns);
            random_number = row*(corr->columns) + column;

            power_c[i] += (it->second)->leakage_power*pow(2.718, -adapt/(1.5*0.026))*(1 - 0.05*corr->random_numbers[random_number][i] + 0.08*corr->random_numbers_width[random_number][i]/3); 
            //power_c[i] += (it->second)->leakage_power; 
        }
        mean += power_c[i];
        //cout<<"POWER_c[i] = "<<power_c[i]<<endl;
     }

     mean = mean/no_of_iterations;

     power_mu_sigma.first = mean;
     power_mu_sigma.second  = get_sigma(power_c, mean, no_of_iterations);
     delete power_c;
     return power_mu_sigma;
 }

 double circuit_class::get_sigma(double *power_c, double mean, int no_of_iterations){
     double st_dev = 0;
       for(int i=0; i<no_of_iterations; i++){
           st_dev += pow((power_c[i] - mean), 2.0);
       }
       st_dev = st_dev/no_of_iterations;
       return sqrt(st_dev);
 }

 /******************************************************************************/
 /******************************************************************************/




// ****************************************************************************** 
// ***************** _______       _______                 ___
// ***************** _______\     /_______\               /___\
// ***************** ||     \\   //       \\             //   \\
// ***************** ||      || //         \\           //     \\
// ***************** ||      || ||         ||          //       \\
// ***************** ||      // ||                    //         \\
// ***************** ||_____//  ||                   //           \\
// ***************** ||_____/   ||                  //=============\\
// ***************** ||         ||                 //               \\
// ***************** ||         ||         ||     //                 \\
// ***************** ||         \\         //    //                   \\
// ***************** ||          \\ ______//    //                     \\
// ***************** ||           \\______/    //                       \\
// This code below will generate Timing Analysis for the PCA method

void circuit_class::SSTA(correlation_model *corr, int no_of_prin_comp){
    //First I will make sure that the output time of every gate is 0

		//ran_1 = gaussian_random_numbers(10000);
		//ran_2 = gaussian_random_numbers(10000);

    for(std::map<string, gate_class*>::iterator it  = this->Gates_Map.begin();
                                                it != this->Gates_Map.end()  ;  ++it){
        (it->second)->output_time_of_gate = 0.0;

        // here We will call the functions to calculate the mean and the variance of the gates
        // depending upon the inherent delay of the gate which is the mean of the delay and the
        // random variable in which it lies will give us the variance !!
        // Since mean is the the normal delay so we will assign it straight forwardly
        (it->second)->output_mu_sigma.first  = 0.0;
        (it->second)->output_mu_sigma.second = 0.0;

				(it->second)->k_parameters.clear();
				(it->second)->output_k_parameters.clear();

        double adapt;

        if(high_low[(it->second)->adaptive_block] == true){
            if(high[(it->second)->adaptive_block] == 1){
                //adapt = -0.05;
                adapt = 1.2;
            } else if (low[(it->second)->adaptive_block] == 1) {
                //adapt = 0.05;
                adapt = 1.2;
            } else {
                //adapt = 0;
                adapt = 1;
            }
        } else {
            //adapt = 0;
            adapt = 1;
        }
        //cout << "ADAPT = " << adapt << endl;
        //(it->second)->mu_sigma.first  = (it->second)->delay_of_gate/pow(1.2-(0.2+adapt),2);
        (it->second)->mu_sigma.first  = (it->second)->delay_of_gate/adapt;
        (it->second)->mu_sigma.second = calculate_sigma_of_gate( (it->second), corr, no_of_prin_comp, adapt );

    }

    //Now we have to iterate over all the gates starting from primary inputs till primary outputs
    //or DFF i.e. D Flips Flops
    //for(std::vector<gate_class*>::iterator it  =  this->Primary_Outputs.begin();
        //                                   it !=  this->Primary_Outputs.end()  ; ++it){
        // First the max operation and then the addition operation 
        // For these we just have to calulate the mu and sigma 
        // Here the timing Analysis is not the stright forward way as we have done
        // in the normal way. Here I want that all the inputs of a gate are evaluated 
        // before i can progress to the new gate 
        //
        // So for a recursive analysis I whave to start from the primary outputs and the D Flip FLops
        // Since I have only the list of the primary outputs I need to generate another list of the D Flip Flops <- TODO
        //calculate_mu_and_sigma(*it, no_of_prin_comp);
        calculate_mu_and_sigma(*(Primary_Inputs.begin()), no_of_prin_comp);	
		
				//cout << "FINSIHED WITH THE DELAY CALCULATIONS !!" << endl;
    //}    
}

double circuit_class::calculate_sigma_of_gate(gate_class* gate, correlation_model *corr, int no_of_prin_comp, double adapt){
    //TODO: Keep all these sigma variations ready for the in the PAC file. We should not calcualte them again and again here
    int   row, column;
    int   random_number;
    row    = floor((gate->y_grid + corr->min_y)/GRID_SIZE); //corr->rows);
    column = floor((gate->x_grid + corr->min_x)/GRID_SIZE);    //corr->columns);
    
    random_number = row*(corr->columns) + column;
/*
    float tox = 1.5e-9;
		float mu =  0.01186;
		float L =  45e-9;
		float W =  135e-9;
		float E =  3.9*8.8542e-12;
		float kdL;
		//kdL =  2*tox*1e-9/(mu*W*E);
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


    DdL = (kdL/Vdd)*p_cal*gate->load_capacitance;
    DdW = (kdW/Vdd)*p_cal*gate->load_capacitance;
  */  
    double sigma = 0.0;
    double temp;


    //cout << gate->gate_name << "  random number SSTTA " <<  random_number << endl;
    //int halt1;
    //cin >> halt1;

    // This one for Length
    for(int i = 0; i < no_of_prin_comp; i++){
				if(corr->eigen_values[i] < 0){
					//corr->eigen_values[i] = -corr->eigen_values[i];
					corr->eigen_values[i] = 0;
				}
        //temp = DdL*(0.05*45/3)*sqrt(corr->eigen_values[i]) * corr->eigen_vectors[i][random_number]/100;
        //gate->delay_of_gate = gate->delay_of_gate/adapt;
        
        //temp = (0.15/3)*gate->delay_of_gate*sqrt(corr->eigen_values[i]) * corr->eigen_vectors[random_number][i]/pow(1.2-(0.2+adapt),2);
        temp = (0.15/3)*gate->delay_of_gate*sqrt(corr->eigen_values[i]) * corr->eigen_vectors[random_number][i]/adapt;

        //temp = sqrt(corr->eigen_values[i]) * corr->eigen_vectors[random_number][i];
        //temp = (0.10)*gate->gate_resistance*sqrt(corr->eigen_values[i]) * corr->eigen_vectors[i][random_number] * gate->load_capacitance;//[random_number][i];
        sigma += temp * temp;
        gate->k_parameters.push_back( temp );
    }

    // This one for Width
    for(int i = 0; i < no_of_prin_comp; i++){
				if(corr->eigen_values_width[i] < 0){
					//corr->eigen_values_width[i] = -corr->eigen_values_width[i];
					corr->eigen_values_width[i] = 0;
				}

        //cout << corr->eigen_values_width[i] << endl;
        //temp = -DdW*(0.01*135/3)*sqrt(corr->eigen_values_width[i]) * corr->eigen_vectors_width[i][random_number]/100;

        //temp = -(0.08/3)*gate->delay_of_gate*sqrt(corr->eigen_values_width[i]) * corr->eigen_vectors_width[random_number][i]/pow(1.2-(0.2+adapt),2);
        temp = -(0.08/3)*gate->delay_of_gate*sqrt(corr->eigen_values_width[i]) * corr->eigen_vectors_width[random_number][i]/adapt;

        //cout << temp << endl;
        //temp = -DdW*(0.08)*sqrt(corr->eigen_values_width[i]) * corr->eigen_vectors_width[i][random_number]* gate->gate_resistance* gate->load_capacitance;//[random_number][i];
        sigma += temp * temp;
        gate->k_parameters.push_back( temp);
    }
    //cout << "The sigma of the gate " << gate->gate_name << " is   " << sqrt(sigma) << "of which the square is "<< sigma << endl;

    return sqrt(sigma);
}

void circuit_class::calculate_mu_and_sigma(gate_class* gate, int no_of_prin_comp){
    //First I will always check whether a particular gate has all the inputs a value or not
    //
    //int wait;
/*
		if(gate->gate_type == "DFF"){ // for dff we need to stop here
			return;
		}

    for(std::vector<gate_class*>::iterator it  = gate->inputs.begin();
                                           it != gate->inputs.end()  ; ++it){
        if( ((*it)->gate_type) == "PI" ){
            for(int j = 0 ; j < 2*no_of_prin_comp; j++){
                (*it)->output_k_parameters.push_back((*it)->k_parameters[j]);
            }
            (*it)->output_mu_sigma.first  = (*it)->mu_sigma.first;
            (*it)->output_mu_sigma.second = (*it)->mu_sigma.second;
            continue;
        } else {
            if( (*it)->output_mu_sigma.first == 0.0 && (*it)->output_mu_sigma.second == 0.0 ){
                calculate_mu_and_sigma(*it, no_of_prin_comp);
            }
        }
    }
*/
 
 //cout << levels.size() << endl;
	//int wait;
	//cin >> wait;
  for (map<int, vector<gate_class*> >::iterator ii=levels.begin(); ii != levels.end(); ++ii ) {
		//cout << ii->first << endl;
		for(vector<gate_class*>::iterator it_vec = (ii->second).begin(); it_vec != (ii->second).end(); ++it_vec){
		//cout << "hello1 " << ii << "   "<< levels[ii].size() <<endl;
		gate = (*it_vec);
    mu_sigma_struct s1, s2;

		//cout << gate->gate_name << " type = " << gate->gate_type << "   " << gate->gate_level<< endl; 

    std::vector<gate_class*>::iterator it1 = gate->inputs.begin();


    if(it1 != gate->inputs.end() && gate->gate_type_add != "DFF"){
        s1.mu           = (*it1)->output_mu_sigma.first;
        s1.sigma        = (*it1)->output_mu_sigma.second;
        s1.k_parameters = (*it1)->output_k_parameters;
		} else {
        s1.mu           = 0;
        s1.sigma        = 0;
        for(int j = 0 ; j < 2*no_of_prin_comp; j++){
            s1.k_parameters.push_back(0.0);
        }
    }
    //Now I have evaluated all the inputs of the gate... Now I will evaluate the max fucntion
    for(std::vector<gate_class*>::iterator it  = gate->inputs.begin();
                                           it != gate->inputs.end()  ; ++it){
        //if( (*it)->gate_type == "PI" ){
        //    continue;
        //}
        if(*it1 == *it){
            continue;
        }

        s2.mu           = (*it)->output_mu_sigma.first;
        s2.sigma        = (*it)->output_mu_sigma.second;
        s2.k_parameters = (*it)->output_k_parameters;
	
				//cout << endl << endl;
			  //cout << "Before calling the Max Function " << endl;
				//cout << (*it)->gate_name << "    " << (*it)->gate_type << "  " << gate->gate_name << endl;
				//cout << (*it)->outputs.size() << "    " << (*it)->inputs.size() << "  " << gate->gate_level << endl;
				//cout << (*it1)->gate_name << "    " << (*it)->gate_name << "  " << gate->gate_name << endl;
				//cout << (*it1)->adaptive_block << "    " << (*it)->adaptive_block << "  " << gate->adaptive_block << endl;
				//int wait;
				//cin >> wait;
        
        s1 = SSTA_Max_Function(s1,s2, no_of_prin_comp);
				//cout << endl << endl << s1.mu << "    " << s1.sigma  << endl;


        if(s1.mu > 10000000 || s1.mu < -10000000 || s1.sigma > 10000000 || s1.sigma < -10000000){
            cout << (*it1)->gate_name  << "       " <<(*it)->gate_name << "      " << gate->gate_name << endl;
            cout << "Either the Mean or Sigma is infinity please check " << endl;
            cout << s1.mu << "    " << s1.sigma << "   "  << endl;
            //exit(1);
        }
    }

    //After we have calculated the max function at the inputs we will add this max fucntion
    //to the delay of the gate.. This means that we have to calulate the mu and the sigma
        if(s1.mu > 10000000 || s1.mu < -10000000 || s1.sigma > 10000000 || s1.sigma < -10000000){
            if(it1 != gate->inputs.end())
                cout << (*it1)->gate_name << "      " << gate->gate_name << endl;
            cout << gate->gate_name << endl;
            cout << gate->gate_type  << "       " << gate->gate_level << endl;
            cout << "Either the Mean or Sigma of S1 is infinity please check " << endl;
            cout << s1.mu << "    " << s1.sigma << "   "  << endl;
            /*int wait;
            cin >> wait;*/
            exit(1);
        }
        if(gate->mu_sigma.first > 10000000 || gate->mu_sigma.first < -10000000 || gate->mu_sigma.second > 10000000 || gate->mu_sigma.second < -10000000){
            cout <<  gate->gate_name << endl;
            cout << "Either the Mean or Sigma of GATE is infinity please check " << endl;
            cout << gate->mu_sigma.first << "    " << gate->mu_sigma.second << "   "  << endl;
            /*int wait;
            cin >> wait;*/
            exit(1);
        }

    SSTA_Addition_Function(gate, s1, no_of_prin_comp);
		}
	}
}

void circuit_class::SSTA_Addition_Function(gate_class* gate, mu_sigma_struct s1, int no_of_prin_comp){
    // This function will calculate the mean for the gate after the delay has been calculated
    // The sum operation will  use the eigen values and the eigen vectors which are available in the 
    // correlation class
    double mu;
    double sigma_sq = 0.0;

    mu = s1.mu + gate->mu_sigma.first;

    double temp;
    for(int i=0; i < 2*no_of_prin_comp; i++){
        gate->output_k_parameters.push_back( gate->k_parameters[i] + s1.k_parameters[i]);

        // The calculation of sigma is given wrong in the paper. The correct form is as follow
        temp = gate->k_parameters[i] + s1.k_parameters[i];
        //temp = pow(gate->k_parameters[i],2) + pow(s1.k_parameters[i],2);
        sigma_sq += pow( temp, 2);
        //sigma_sq += temp;
        //gate->output_k_parameters.push_back( sqrt(temp) );
    }

    if(mu > 10000000 || mu < -10000000 || sigma_sq > 10000000){
        cout << "in add function the mu and sigma are not correct " << mu << " and  sigma square is " << sigma_sq << endl;
        //exit(1);
    }

    //cout << "Input sigmas are " << gate->mu_sigma.second << "    " << s1.sigma << endl;
    //cout << "IN addition function " << sqrt(sigma_sq) << endl;

    gate->output_mu_sigma.first  = mu;
    gate->output_mu_sigma.second = sqrt(sigma_sq);

}

mu_sigma_struct circuit_class::SSTA_Max_Function(mu_sigma_struct gate_1, mu_sigma_struct gate_2, int no_of_prin_comp){
    double mu;
    double sigma;
    double co_variance = 0.0;
    double correlation;

    double alpha;
    double beta;

    double phi_plus;
    double phi_minus;
    double zhi;

    double mu1, mu2, sigma1, sigma2;

    vector<double> max_k_parameters;
    vector<double> a_r;

    mu_sigma_struct s_ret;

    mu1 = gate_1.mu;
    mu2 = gate_2.mu;

    sigma1 = gate_1.sigma;
    sigma2 = gate_2.sigma;

    // if one of them is the primary input the the mu and sigma will be zero for that gate
    // So in this condition out equations will produce wrong results as division by zero will
    // always produce infinity (NAN). So we should return from here without going further.
   
		if(mu1 - 3*sigma1 > mu2 + 3*sigma2){
        s_ret.mu = mu1;
        s_ret.sigma = sigma1;
        s_ret.k_parameters = gate_1.k_parameters;
				//cout << "Returning Early " << endl;
        return s_ret;
		}


		if(mu2 - 3*sigma2 > mu1 + 3*sigma1){
        s_ret.mu = mu2;
        s_ret.sigma = sigma2;
        s_ret.k_parameters = gate_2.k_parameters;
				//cout << "Returning Early " << endl;
        return s_ret;
		}

    
    // 1) The mu and the sigma of all the gates has been calculated
    // 2) Find the correlation between the delay of gate_1 and gate_2
    // 2.1) Find the co-variance of delays of gate_1 and gate_2
    double dummy_sigma1 = 0.0;
    double dummy_sigma2 = 0.0;
    co_variance = 0;
    correlation = 0;
    for(int i = 0; i < 2*no_of_prin_comp; i++){
        co_variance += gate_1.k_parameters[i] * gate_2.k_parameters[i];

        //this part of code is only for testing purpose and will be removed later on <- TODO
        dummy_sigma1 += gate_1.k_parameters[i] * gate_1.k_parameters[i];
        dummy_sigma2 += gate_2.k_parameters[i] * gate_2.k_parameters[i];
        
    }
    

    char inp;
    if(abs(gate_1.sigma - sqrt(dummy_sigma1)) > 1){
        cout << "The gate 1 sigmas are not equal " << sqrt(dummy_sigma1) << "   " <<  gate_1.sigma << endl;
        cout << "Enter q to quit " << endl;
        /*cin >> inp;
        if(inp == 'q'){
            exit(1);
        }*/
    }


    if(abs(gate_2.sigma - sqrt(dummy_sigma2)) > 1){
        cout << "The gate 1 sigmas are not equal " << sqrt(dummy_sigma2) << "   " <<  gate_2.sigma << endl;
        cout << "Enter q to quit " << endl;
        /*cin >> inp;
        if(inp == 'q'){
            exit(1);
        }*/
    }

    if(sigma1 == 0 && sigma2 != 0 || sigma2 == 0 && sigma1 != 0 ) {// || sigma2 == 0){ // It is a primary input and most likely the mean will be zero
        correlation = 0;
    } else {
        correlation = co_variance/(gate_1.sigma * gate_2.sigma);
    }
    
    if(sigma1 == 0 && sigma2 == 0){
        s_ret.mu = mu1 > mu2 ? mu1 : mu2;
        s_ret.sigma = 0;
        for(int m=0; m < 2*no_of_prin_comp; m++){
            s_ret.k_parameters.push_back(0.0);
        }
        return s_ret;
    }

    int corr_int;
    corr_int = (abs(correlation) + 0.01);

    if(correlation > 0.99 && abs(sigma1 - sigma2) < 0.1){ // && gate_1.sigma == gate_2.sigma){
        s_ret.mu = mu1 > mu2 ? mu1 + 0.0 : mu2 +0.0;
        s_ret.sigma = mu1 > mu2 ? sigma1 - 0.0 : sigma2 - 0.0;
        s_ret.k_parameters = mu1 > mu2 ? gate_1.k_parameters : gate_2.k_parameters;
        return s_ret;
    } //else { 
    
    // 3) Calculate the mean and the sigma for the max function using the equations as mentioned in the paper

    alpha = sqrt( pow(sigma1,2) + pow(sigma2,2) - 2*co_variance );
    beta  = (mu1-mu2)/alpha;
    
    zhi       = pow(2.718281828, -beta*beta/2)/sqrt(2*3.141592654);
    phi_plus  = integrate(beta);

    if(phi_plus > 1){
        phi_plus = 1.0;
				phi_minus = 0;
    }

		phi_minus = 1 - phi_plus;

		double sigma3;
    mu    = mu1*phi_plus + mu2*phi_minus + alpha*zhi;
    sigma3 = (pow(mu1,2) + pow(sigma1,2))*phi_plus + (pow(mu2,2) + pow(sigma2,2))*phi_minus + (mu1+mu2)*alpha*zhi;
		sigma = sigma3 - mu*mu;
		if(sigma < 0){
        sigma = 0;
		}
    sigma = sqrt(sigma);


    // 4) Calculate the a_r co-efficients
    // Here we will use the equation 32 in the paper to calculate the co-efficients
    double rho1;
    double rho2;

    double temp_cor;
		double temp_cov;
    double temp_a_r;
    double S0_2 = 0.0;
    for(int i = 0; i < 2*no_of_prin_comp; i++){
        rho1      = gate_1.k_parameters[i];
        rho2      = gate_2.k_parameters[i];
        temp_cov  = ( rho1*phi_plus + rho2*phi_minus );
        temp_a_r  = temp_cov ;//* sigma;
        S0_2     += temp_a_r * temp_a_r;  
        
        a_r.push_back(temp_a_r);
    }

    for(int i = 0; i < 2*no_of_prin_comp; i++){
        a_r[i] = a_r[i]*sigma/sqrt(S0_2);
    }

/*
  	cout << endl << endl;
		cout << "zhi = " << zhi << endl;
		cout << "alpha = " << alpha << endl;
		cout << "phi_plus = " << phi_plus << endl;
		cout << "phi_minus = " << phi_minus << endl << endl;

		cout << "Correlation = " << correlation << endl;
		cout << "Mu  = " << mu << endl;
		cout << "Sigma = " << sigma << endl;

		cout << "The input mu and sigma are as follow " << endl;
		cout << "mu1 \t\t mu2 \t\t sigma1 \t\t sigma2 " << endl; 
		cout << mu1 << " \t\t " << mu2 << " \t\t " << sigma1 << " \t\t " << sigma2 << endl;
		cout << "the calculated mu and sigma for the max operation = " << mu  << "   " << sigma << endl;

*/


		if(sigma > sigma1+0.2 && sigma > sigma2+0.2){
			cout << "CRISIS " << endl;
			//int wait;
			//cin >> wait;
		}

    s_ret.mu = mu;
    s_ret.sigma = sigma;
    s_ret.k_parameters = a_r;
    
    return s_ret;
}

double circuit_class::get_max_delay_prob(double threshold){
    double delay = 0.0;
    gate_class* gate = NULL;
    for(std::vector<gate_class*>::iterator it  = Primary_Outputs.begin();
                                           it != Primary_Outputs.end()  ; ++it){
        if(delay <  (*it)->output_mu_sigma.first + 3*(*it)->output_mu_sigma.second){
            delay = (*it)->output_mu_sigma.first + 3*(*it)->output_mu_sigma.second;
            gate = *it;
        }
    }

    //cout << "Delay values are " << gate->output_mu_sigma.first << "   " << gate->output_mu_sigma.second << endl;

    return integrate(threshold, gate->output_mu_sigma.first, gate->output_mu_sigma.second);
}


/******************************************************************************************************************/
/**               THIRD METHOD                                                                                   **/
/******************************************************************************************************************/

void circuit_class::SSTA_Minus_Function(gate_class* gate, mu_sigma_struct s1, int no_of_prin_comp){
    double mu;
    double sigma_sq = 0.0;

    mu = s1.mu - gate->mu_sigma.first;

    double temp;
    for(int i=0; i < 2*no_of_prin_comp; i++){
        temp = s1.k_parameters[i] - gate->k_parameters[i];
        gate->rat_input_k_parameters.push_back(temp);

        sigma_sq += pow( temp, 2);
    }

    gate->rat_input_mu_sigma.first  = mu;
    gate->rat_input_mu_sigma.second = sqrt(sigma_sq);
}

mu_sigma_struct circuit_class::SSTA_Minimum_Function(mu_sigma_struct gate_1, mu_sigma_struct gate_2, int no_of_prin_comp){
    double mu;
    double sigma;
    double co_variance = 0.0;
    double correlation;

    double alpha;
    double beta;

    double phi_plus;
    double phi_minus;
    double zhi;

    double mu1, mu2, sigma1, sigma2;

    vector<double> max_k_parameters;
    vector<double> a_r;

    mu_sigma_struct s_ret;

    // Since we are taking the negative of the delays. it only make sense that the mean is negative
    // I cannot comprehend the variance to be negative. The variance will remain as it is. Now lets say
    // that we have a gaussian numbers now what we do is take a negative of all these numbers. This will 
    // shift all the numbers on the number line towards the left side of the origin. Now if we calculate the mean
    // it will be negative but the variance will still remain the same !!
    mu1 = -gate_1.mu;
    mu2 = -gate_2.mu;

    sigma1 = gate_1.sigma;
    sigma2 = gate_2.sigma;
   
		if(mu1 - 3*sigma1 > mu2 + 3*sigma2){
        s_ret.mu = -mu1;
        s_ret.sigma = sigma1;
        s_ret.k_parameters = gate_1.k_parameters;
        return s_ret;
		}


		if(mu2 - 3*sigma2 > mu1 + 3*sigma1){
        s_ret.mu = -mu2;
        s_ret.sigma = sigma2;
        s_ret.k_parameters = gate_2.k_parameters;
        return s_ret;
		}

    
    double dummy_sigma1 = 0.0;
    double dummy_sigma2 = 0.0;
    co_variance = 0;
    correlation = 0;
    for(int i = 0; i < 2*no_of_prin_comp; i++){
        co_variance += gate_1.k_parameters[i] * gate_2.k_parameters[i];

        dummy_sigma1 += gate_1.k_parameters[i] * gate_1.k_parameters[i];
        dummy_sigma2 += gate_2.k_parameters[i] * gate_2.k_parameters[i];
        
    }
    

    char inp;
    if(abs(gate_1.sigma - sqrt(dummy_sigma1)) > 1){
        cout << "The gate 1 sigmas are not equal " << sqrt(dummy_sigma1) << "   " <<  gate_1.sigma << endl;
        cout << "Enter q to quit " << endl;
        /*cin >> inp;
        if(inp == 'q'){
            exit(1);
        }*/
    }


    if(abs(gate_2.sigma - sqrt(dummy_sigma2)) > 1){
        cout << "The gate 1 sigmas are not equal " << sqrt(dummy_sigma2) << "   " <<  gate_2.sigma << endl;
        cout << "Enter q to quit " << endl;
        /*cin >> inp;
        if(inp == 'q'){
            exit(1);
        }*/
    }

    if(sigma1 == 0 && sigma2 != 0 || sigma2 == 0 && sigma1 != 0 ) {// || sigma2 == 0){ // It is a primary input and most likely the mean will be zero
        correlation = 0;
    } else {
        correlation = co_variance/(gate_1.sigma * gate_2.sigma);
    }
    
    if(sigma1 == 0 && sigma2 == 0){
        s_ret.mu = -(mu1 > mu2 ? mu1 : mu2);
        s_ret.sigma = 0;
        for(int m=0; m < 2*no_of_prin_comp; m++){
            s_ret.k_parameters.push_back(0.0);
        }
        return s_ret;
    }

    int corr_int;
    corr_int = (abs(correlation) + 0.01);

    if(correlation > 0.99 && abs(sigma1 - sigma2) < 0.1){ // && gate_1.sigma == gate_2.sigma){

        s_ret.mu = -(mu1 > mu2 ? mu1 + 0.0 : mu2 +0.0);
        s_ret.sigma = mu1 > mu2 ? sigma1 - 0.0 : sigma2 - 0.0;
        s_ret.k_parameters = mu1 > mu2 ? gate_1.k_parameters : gate_2.k_parameters;
        return s_ret;
    } 
    
    alpha = sqrt( pow(sigma1,2) + pow(sigma2,2) - 2*co_variance );
    beta  = (mu1-mu2)/alpha;
    

    zhi       = pow(2.718281828, -beta*beta/2)/sqrt(2*3.141592654);
    phi_plus  = integrate(beta);

    if(phi_plus > 1){
        phi_plus = 1.0;
				phi_minus = 0;
    }

		phi_minus = 1 - phi_plus;

		double sigma3;
    mu    = mu1*phi_plus + mu2*phi_minus + alpha*zhi;
    sigma3 = (pow(mu1,2) + pow(sigma1,2))*phi_plus + (pow(mu2,2) + pow(sigma2,2))*phi_minus + (mu1+mu2)*alpha*zhi;
		sigma = sigma3 - mu*mu;
		if(sigma < 0){
        sigma = 0;
		}
    sigma = sqrt(sigma);

    double rho1;
    double rho2;

    double temp_cor;
		double temp_cov;
    double temp_a_r;
    double S0_2 = 0.0;
    for(int i = 0; i < 2*no_of_prin_comp; i++){
        rho1      = gate_1.k_parameters[i];
        rho2      = gate_2.k_parameters[i];
        temp_cov  = ( rho1*phi_plus + rho2*phi_minus );
        temp_a_r  = temp_cov ;//* sigma;
        S0_2     += temp_a_r * temp_a_r;  
        a_r.push_back(temp_a_r);
    }

    for(int i = 0; i < 2*no_of_prin_comp; i++){
        a_r[i] = a_r[i]*sigma/sqrt(S0_2);
    }

    s_ret.mu = -mu;
    s_ret.sigma = sigma;
    s_ret.k_parameters = a_r;
    return s_ret;
}

void circuit_class::Calculate_RAT(int no_of_prin_comp, double required_output_time){
    for(map<string, gate_class*>::iterator it  = Gates_Map.begin();
                                           it != Gates_Map.end()  ; ++it){
        it->second->rat_mu_sigma.first  = required_output_time;
        it->second->rat_mu_sigma.second = 0;
        it->second->rat_input_mu_sigma.first  = required_output_time;
        it->second->rat_input_mu_sigma.second = 0;
    }

    for(std::vector<gate_class*>::iterator it  = Primary_Inputs.begin();
                                           it != Primary_Inputs.end()  ; ++it){
        calculate_rat_recursively(*it, no_of_prin_comp, required_output_time);
    }
}

void circuit_class::calculate_rat_recursively(gate_class *gate, int no_of_prin_comp, double required_output_time){
    mu_sigma_struct s1, s2;

    if(gate->gate_type == "PO"){
        s1.mu  = required_output_time;
        s1.sigma = 0;
        for(int i=0; i<2*no_of_prin_comp; i++){
            s1.k_parameters.push_back(0);
        }
    }

    for(vector<gate_class*>::iterator it  = gate->outputs.begin();
                                      it != gate->outputs.end()  ; ++it){
        if((*it)->rat_mu_sigma.first == required_output_time && (*it)->gate_type_add != "DFF"){
            calculate_rat_recursively(*it, no_of_prin_comp, required_output_time);
        }

        if((*it)->gate_type_add == "DFF"){
            (*it)->rat_input_mu_sigma.first  = required_output_time;
            (*it)->rat_input_mu_sigma.second = 0;
            for(int i=0; i<2*no_of_prin_comp; i++){
                (*it)->rat_input_k_parameters.push_back(0);
            }
        }

        if(it == gate->outputs.begin()){
            s1.mu = gate->outputs[0]->rat_input_mu_sigma.first;
            s1.sigma = gate->outputs[0]->rat_input_mu_sigma.second;
            s1.k_parameters = gate->outputs[0]->rat_input_k_parameters;
            continue;
        }

        s2.mu = (*it)->rat_input_mu_sigma.first;
        s2.sigma = (*it)->rat_input_mu_sigma.second;
        s2.k_parameters = (*it)->rat_input_k_parameters;

        s1 = SSTA_Minimum_Function(s1, s2, no_of_prin_comp);
    }

    gate->rat_mu_sigma.first  = s1.mu;
    gate->rat_mu_sigma.second = s1.sigma;
    gate->rat_k_parameters    = s1.k_parameters;

    SSTA_Minus_Function(gate, s1, no_of_prin_comp);
}



/***********************************************************************************************/
/**                This function will remove those gates which will meet timing         **/
/***********************************************************************************************/

void circuit_class::Prune_nodes(){
    //for(vector<gate_class*>::iterator it  = Primary_Inputs.begin();
    //                                  it != Primary_Inputs.end()  ; ++it){
        prune_nodes_recursively(*(Primary_Inputs.begin()));
    //}
}

void circuit_class::prune_nodes_recursively(gate_class* gate){
    vector<gate_class*> temp_Gates_Vec;
    map<string, gate_class*>::iterator it1;
    vector<gate_class*>::iterator it_prune;

    for(it1 = Gates_Map.begin(); it1 != Gates_Map.end(); ++it1){
        gate = it1->second;
        if( gate->output_mu_sigma.first + 3*(gate->output_mu_sigma.second) <= gate->rat_mu_sigma.first -3*(gate->rat_mu_sigma.second) ){
            cout << gate->gate_name << "  " << gate->output_mu_sigma.first + 3*(gate->output_mu_sigma.second) << "  " << gate->rat_mu_sigma.first -3*(gate->rat_mu_sigma.second) << endl; 
            temp_Gates_Vec.push_back(gate);
        }
    }

    cout << "TO BE DELETED = " << temp_Gates_Vec.size() << endl;
    //int wait;
    //cin >> wait;

    while(temp_Gates_Vec.size() != 0){
        gate = temp_Gates_Vec[0];
        for(vector<gate_class*>::iterator it  = gate->inputs.begin();
                                          it != gate->inputs.end()  ; ++it){
            //cout << "the gate is " << gate->gate_name << "  and the input name is " << (*it)->gate_name << endl;
            it_prune = find((*it)->outputs.begin(), (*it)->outputs.end(), gate);

            (*it)->outputs.erase(it_prune);
            if((*it)->outputs.size() == 0 && ((*it)->gate_type_add != "DFF" && gate->gate_type_add != "DFF") ){
                it_prune = find( temp_Gates_Vec.begin(), temp_Gates_Vec.end(), (*it));
                if(it_prune == temp_Gates_Vec.end() ){
            cout << "TT " << (*it)->gate_name << "  " << (*it)->output_mu_sigma.first + 3*((*it)->output_mu_sigma.second) << "  " << (*it)->rat_mu_sigma.first -3*((*it)->rat_mu_sigma.second) << endl; 
                    temp_Gates_Vec.push_back(*it);
                }
            }
        }

        for(vector<gate_class*>::iterator it  = gate->outputs.begin();
                                          it != gate->outputs.end()  ; ++it){
            it_prune = find((*it)->inputs.begin(), (*it)->inputs.end(), gate);
            (*it)->inputs.erase(it_prune);
            if((*it)->inputs.size() == 0 && (*it)->gate_type_add != "DFF" && gate->gate_type_add != "DFF" ){
                it_prune = find( temp_Gates_Vec.begin(), temp_Gates_Vec.end(), (*it));
                if(it_prune == temp_Gates_Vec.end() ){
            cout << "PP " << (*it)->gate_name << "  " << (*it)->output_mu_sigma.first + 3*((*it)->output_mu_sigma.second) << "  " << (*it)->rat_mu_sigma.first -3*((*it)->rat_mu_sigma.second) << endl; 
                    temp_Gates_Vec.push_back(*it);
                }
            }
        }
        //cout << "Going to delete a gate " << gate->gate_name << endl;
        //int wait;
        //cin >> wait;
        it_prune = find(Primary_Inputs.begin(), Primary_Inputs.end(), gate);
        if(it_prune != Primary_Inputs.end() ){
            //cout << "MIL GAYA GATE I " << gate->gate_name<<endl;
            Primary_Inputs.erase(it_prune);
            //int wait;
            //cin >> wait;
        }
        it_prune = find(Primary_Outputs.begin(), Primary_Outputs.end(), gate);
        if(it_prune != Primary_Outputs.end() ){
            //cout << "MIL GAYA GATE O " << gate->gate_name<<endl;
            Primary_Outputs.erase(it_prune);
            //int wait;
            //cin >> wait;
        }
        int dd;
        dd = gate->gate_level;
        it_prune = find(levels[dd].begin(), levels[dd].end(), gate);
        if(it_prune != levels[dd].end()) {
            levels[dd].erase(it_prune);
        }
        this->Gates_Map.erase(gate->gate_name);
        temp_Gates_Vec.erase(temp_Gates_Vec.begin());
        delete gate;
    }
}

/*void circuit_class::prune_nodes_recursively(gate_class* gate){

    for(vector<gate_class*>::iterator it = gate->outputs.begin(); it != gate->outputs.end(); ++it){
        prune_nodes_recursively(*it);
    }
    
    vector<gate_class*>::iterator it_prune;

    if( gate->output_mu_sigma.first + 3*gate->output_mu_sigma.second <= gate->rat_mu_sigma.first -3*gate->rat_mu_sigma.second ||
        gate->outputs.size() == 0 && gate->gate_type != "PO"){
        for(vector<gate_class*>::iterator it  = gate->inputs.begin();
                                          it != gate->inputs.end()  ; ++it){
            it_prune = find((*it)->outputs.begin(), (*it)->outputs.end(), gate);
            //cout << "THe gate in question is " << gate->gate_name << " and the problematic input is " << (*it)->gate_name << endl;
            (*it)->outputs.erase(it_prune);
        }
        
        for(vector<gate_class*>::iterator it  = gate->outputs.begin();
                                          it != gate->outputs.end()  ; ++it){
            it_prune = find((*it)->inputs.begin(), (*it)->inputs.end(), gate);
            //cout << "THe gate in question is " << gate->gate_name << " and the problematic input is " << (*it)->gate_name << endl;
            (*it)->inputs.erase(it_prune);
            if((*it)->inputs.size() == 0){
                prunes_gate_in_line(*it);
            }
        }
        cout << "Going to delete a gate " << gate->gate_name << endl;
        this->Gates_Map.erase(gate->gate_name);
        delete gate;
    }

    //if(gate->outputs.size() == 0 && gate->gate_type)
}
*/

void circuit_class::prunes_gate_in_line(gate_class* gate){
    vector<gate_class*>::iterator it_prune;

    for(vector<gate_class*>::iterator it  = gate->outputs.begin();
                                      it != gate->outputs.end()  ; ++it){
            it_prune = find((*it)->inputs.begin(), (*it)->inputs.end(), gate);
            (*it)->inputs.erase(it_prune);
            if((*it)->inputs.size() == 0){
                prunes_gate_in_line(*it);
            }
    }
    //cout << "Hello Something here " << endl;
    this->Gates_Map.erase(gate->gate_name);
    delete gate;
}

void circuit_class::compress_gates(){
    vector<gate_class*> del_gates;
    vector<gate_class*>::iterator it1;
    vector<string> del_gates_2;
    for(map<string,gate_class*>::iterator it = Gates_Map.begin(); it != Gates_Map.end(); ++it){
        it1 = find(del_gates.begin(), del_gates.end(), it->second);
        if(it1 == del_gates.end()){
            compress_gates_recursively(it->second, del_gates);
        } else {
            del_gates_2.push_back(it->first);
        }
    }
    
    while(del_gates_2.size() != 0){
        Gates_Map.erase(del_gates_2[0]);
        //cout << del_gates_2[0] << endl;
        del_gates_2.erase(del_gates_2.begin());
    }
}

void circuit_class::compress_gates_recursively(gate_class* gate, vector<gate_class*> &del_gates){
    vector<gate_class*> it_rep;
    vector<gate_class*>::iterator it_rep1;

    if(gate->outputs.size() == 1 && gate->outputs[0]->adaptive_block == gate->adaptive_block && gate->outputs[0]->inputs.size() == 1 && gate->outputs[0]->gate_type != "PO"){
        //cout << gate->gate_name << endl;
        compress_gates_recursively(gate->outputs[0], del_gates);
        it_rep = gate->outputs[0]->outputs;
        for(vector<gate_class*>::iterator it = it_rep.begin(); it != it_rep.end(); ++it){
            it_rep1 = find((*it)->inputs.begin(), (*it)->inputs.end(), gate->outputs[0]);
            *it_rep1 = gate;
        }
        //TODO: add the delays in such cases
        del_gates.push_back(gate->outputs[0]);
        gate->outputs = it_rep;
    }
}

/*********************************************************************************************/
/*********   Monte Carlo to check R.A.T is correct or not !!             *********************/
/*********************************************************************************************/

void circuit_class::Monte_Carlo_RAT(correlation_model *corr, int no_of_iterations, bool adaptive, double required_output_time){

    FILE *fp;
    gate_class* temp_gate = NULL;
    fp = fopen("delay_RAT.txt", "w");
    double prob;
    float bin_low = 0;
    float bin_high = 0;

    for(int i=0; i<no_of_iterations; i++){
        for(std::map<string, gate_class*>::iterator it1  = Gates_Map.begin();
                                                    it1 != Gates_Map.end()  ; it1++){
            (it1->second)->rat = required_output_time;
        }

        for(vector<gate_class*>::iterator it  = Primary_Inputs.begin();
                                          it != Primary_Inputs.end()  ; ++it){
            Run_timing_rat(*it, corr, i, required_output_time);
        }

        double delay;

        delay = 0.0;

        //temp_gate = Gates_Map["G16"];
        //temp_gate = Gates_Map["N6288"];
        //temp_gate = Gates_Map["N10014"];
        temp_gate = Gates_Map["n_103688"];
        //cout << "The Delay for the gate " << temp_gate->gate_name << "   " << delay <<endl;
        cout << endl << "HELLO " << endl;
        cout << "The mu for gate " << temp_gate->gate_name << "   " << temp_gate->rat_mu_sigma.first <<endl;
        cout << "The sigma for gate " << temp_gate->gate_name << "   " << temp_gate->rat_mu_sigma.second <<endl;

        fprintf(fp,"%f\n",temp_gate->rat);
    }
    cout << "The RAT mu for gate " << temp_gate->gate_name << "   " << temp_gate->rat_mu_sigma.first <<endl;
    cout << "The RAT sigma for gate " << temp_gate->gate_name << "   " << temp_gate->rat_mu_sigma.second <<endl;
    cout << "The mu for gate " << temp_gate->gate_name << "   " << temp_gate->output_mu_sigma.first <<endl;
    cout << "The sigma for gate " << temp_gate->gate_name << "   " << temp_gate->output_mu_sigma.second <<endl;
    cout << temp_gate->output_mu_sigma.first + 3*temp_gate->output_mu_sigma.second << "       " << temp_gate->rat_mu_sigma.first -3*temp_gate->rat_mu_sigma.second  << endl;
    fclose(fp);
 }


 double circuit_class::Run_timing_rat(gate_class* gate, correlation_model *corr, int iteration_number, double required_output_time){
    int   row, column;
    int   random_number;
    row    = floor((gate->y_grid + corr->min_y)/GRID_SIZE); //corr->rows);
    column = floor((gate->x_grid + corr->min_x)/GRID_SIZE);    //corr->columns);
    random_number = row*(corr->columns) + column;

    int adapt;

    double temp_min;
    temp_min = required_output_time;
    //if(gate->gate_type == "PI"){
    //    return corr->random_numbers[random_number][iteration_number];
    //    temp_min = 0;
    //}

    //double temp_max1 = 0;

    for(vector<gate_class*>::iterator it  = gate->outputs.begin();
                                      it != gate->outputs.end()  ; ++it){
       if((*it)->rat == required_output_time && (*it)->gate_type_add != "DFF" ){
           Run_timing_rat(*it, corr, iteration_number, required_output_time);
       }
       if((*it)->gate_type_add == "DFF"){
            (*it)->rat_input = required_output_time;
       }
       if( (*it)->rat_input < temp_min){
           temp_min = (*it)->rat_input;
       }
       
    }

    gate->rat = temp_min;

    adapt = 1;

    if(iteration_number == -1){
        //gate->rat =  temp_min + gate->mu_sigma.first + 3*gate->mu_sigma.second ;
    } else {
        gate->rat_input = gate->rat - gate->delay_of_gate/adapt - gate->delay_of_gate*(0.15/3)*(corr->random_numbers[random_number][iteration_number])/adapt + 1*(corr->random_numbers_width[random_number][iteration_number])*gate->delay_of_gate*(0.08/3)/adapt;
    }
	return 0.0;
 }


#endif
