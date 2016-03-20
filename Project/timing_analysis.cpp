#ifndef TIMING_ANALYSIS_CPP
#define TIMING_ANALYSIS_CPP

#include "define.h"
#include "parser.h"
//#include "correlation_model.h"
//#include "circuit_class.h"
#include "timing_analysis.h"
#include "iostream"
#include "stdlib.h"
#include "string.h"
#include "ispd_parse.h"
#include <stdio.h>
#include <math.h>

//#include "matrix_multiplication.h"
//#include "gaussian_random_numbers.h"

#include "Eigen/Eigenvalues"

double req_time = 0.0;

using namespace Eigen;

extern int parse_netlist(const char*);

void Timing_Analysis::Init_Timing_Analysis(string benchmark, string placement_file, map<string, vector<float> >  gate_parameters, double required_output_time, int x_adapt_blocks, int y_adapt_blocks, int *actual_blocks, bool only_FBB1){
    ckt  = new circuit_class;
    corr = new correlation_model;
	
		only_FBB = only_FBB1;

    // The correlation levels are fixed as of now and I have to see what could be done with these !! 
    //const float corr_levels[] = {0.80, 0.70, 0.52, 0.28, 0.14,0.11, 0.06,0.02,0.01 };
    //const float corr_levels1[] = {0.80, 0.70, 0.52, 0.28, 0.14,0.11, 0.06,0.02,0.01 };
    const float corr_levels[] = {0.3, 0.2, 0.15, 0.01 };
    const float corr_levels1[] = {0.3, 0.2, 0.15, 0.01 };
    //const float corr_levels[] = {0.0, 0.0, 0.0, 0.00 };
    //const float corr_levels1[] = {0.0, 0.0, 0.0, 0.00 };

	// Parse the Netlist
#ifdef IS_ISCAS85_TEST
  	parse_netlist(ckt, benchmark.c_str());
#else
  ispd_parse_verilog(benchmark, ckt);
	map<string, double> library;
	//ispd_parse_lib("contest.lib", library);
	//generate_placement(ckt, benchmark, library);

	cout << "The  gates are " << ckt->Gates_Map.size() << endl;
#endif

/*
 std::vector<gate_class*>::iterator it1;
   for(map<string, gate_class*>::iterator it = ckt->Gates_Map.begin(); it != ckt->Gates_Map.end(); ++it){
			//cout << it->second->gate_name << "   " << it->second->gate_type << "    " << it->second->gate_type_add << endl;
			//cout << it->second->inputs[0] << "   " << it->second->outputs[0] << endl;
			for(it1 = it->second->inputs.begin(); it1 != it->second->inputs.end(); ++it1){
				//cout << "Inputs = " << (*it1)->gate_name << endl;
			}

			for(it1 = it->second->outputs.begin(); it1 != it->second->outputs.end(); ++it1){
				//cout << "Outputs = " << (*it1)->gate_name << endl;
			}
		//cout << endl << endl;
		cout << "The  input are " << ckt->Primary_Inputs.size() << endl;
		cout << "The  outputs are " << ckt->Primary_Outputs.size() << endl;
		int wait;
		cin >> wait;
		}

  exit(1);

*/

	// After parsing the netlist, Levelize the circuit. By Levelizing I mean that all the primary inputs are at level 0 and 
	// then the gates connected to only primary inputs are at level 1.
	ckt->Levelize_Circuit(ckt);
    
    // Update gate delays of each gate using the elmore delay model. This is a part of the 
    // initialization process !!
    ckt->update_gate_delays(gate_parameters, required_output_time);
    //ckt->calculate_required_arrival_time(required_output_time);
    
    // Parsing a placement file and generating a correlation matrix !!
    parse_placement_file(placement_file.c_str(), ckt, size_of_die);

    corr->x_grids = 0;
    corr->y_grids = 0;

    if(size_of_die != NULL){
        //corr->x_grids = ceil( (size_of_die[0] - size_of_die[2]) /GRID_SIZE);
        //corr->y_grids = ceil( (size_of_die[1] - size_of_die[3]) /GRID_SIZE);
        corr->x_grids = ceil( (size_of_die[0]) /GRID_SIZE);
        corr->y_grids = ceil( (size_of_die[1]) /GRID_SIZE);
        cout << endl << "The size of the die s1196 is "<< corr->x_grids << " for x and " << corr->y_grids << " for y " << endl << endl;
    }

    corr->rows    = corr->y_grids;
    corr->columns = corr->x_grids;
		if(size_of_die[2] < 0){
    	corr->min_x   = -floor(size_of_die[2]);
		} else {
    	corr->min_x   = 0;
		}
		if(size_of_die[3] < 0){
    	corr->min_y   = -floor(size_of_die[3]);
		} else {
    	corr->min_y   = 0;
		}


#ifdef QUAD_TREE_MODEL
    corr->generate_correlation_matrix_quadtree(corr->x_grids,corr->y_grids,corr_levels,3,"length");
    corr->generate_correlation_matrix_quadtree(corr->x_grids,corr->y_grids,corr_levels1,3,"width");
#endif
//GRID_MODEL
    corr->generate_correlation_matrix_grid(corr->x_grids,corr->y_grids,corr_levels,3,"length");
    corr->generate_correlation_matrix_grid(corr->x_grids,corr->y_grids,corr_levels1,3,"width");
//#endif


// TODO: The following has to be updated from the optimization folks. For now I am using my own
// Now a new change . Earlier the block were only four now i have to make it as a parameter !!
		if(x_adapt_blocks > corr->x_grids){
			x_adapt_blocks = corr->x_grids;
		}
	
		if(y_adapt_blocks > corr->y_grids){
			y_adapt_blocks = corr->y_grids;
		}

		int x,y;

		x = ceil((float)corr->x_grids/x_adapt_blocks);
		y = ceil((float)corr->y_grids/y_adapt_blocks);
		
		//cout << "x = " << x << " y = " << y << endl;

		ckt->high_low = new bool[x_adapt_blocks*y_adapt_blocks];
    ckt->no_of_adaptive_blocks = x_adapt_blocks*y_adapt_blocks;
		*actual_blocks = ckt->no_of_adaptive_blocks;

		gate_class* gate1;
    for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
                                                it != ckt->Gates_Map.end()  ;  ++it){
			gate1 = it->second;
			gate1->adaptive_block = (int)(gate1->x_grid/GRID_SIZE)/x + (int)(gate1->y_grid/GRID_SIZE)/y * (x_adapt_blocks);
	/*		cout << (gate1->x_grid/GRID_SIZE)/x << "     " << (gate1->y_grid/GRID_SIZE)/y << endl;
			cout << (gate1->x_grid/GRID_SIZE) << "     " << (gate1->y_grid/GRID_SIZE) << endl;
			cout << corr->x_grids/2 << "     " << corr->y_grids/2 << endl; */

  //      if( ((it)->second)->x_grid/GRID_SIZE < x /*(corr->x_grids/2)*/ && ((it)->second)->y_grid/GRID_SIZE < y/*(corr->y_grids/2) */){
  //          if(((it)->second)->adaptive_block != 0) cout << "ERROR in block assignment 0  " << gate1->adaptive_block << endl;
  //      } else if ( ((it)->second)->x_grid/GRID_SIZE >= x /*(corr->x_grids/2)*/ && ((it)->second)->y_grid/GRID_SIZE  < y/*(corr->y_grids/2)*/){
  //       if ( ((it)->second)->adaptive_block != 1) cout << "ERROR in block assignment 1  " <<  gate1->adaptive_block << endl;
  //      } else if ( ((it)->second)->x_grid/GRID_SIZE < x /*(corr->x_grids/2)*/ && ((it)->second)->y_grid/GRID_SIZE >= y/*(corr->y_grids/2)*/){
  //        if ( ((it)->second)->adaptive_block != 2) cout << "ERROR in block assignment 2  " <<  gate1->adaptive_block << endl;
  //      } else if ( ((it)->second)->x_grid/GRID_SIZE >= x/*(corr->x_grids)/2*/ && ((it)->second)->y_grid/GRID_SIZE >= y/*(corr->y_grids/2)*/){
  //         if( ((it)->second)->adaptive_block != 3) cout << "ERROR in block assignment 3  " << gate1->adaptive_block << endl;
  //      }
	//	int wait1;
	//	cin >> wait1;  
		}
	
		
/*
        if( ((it)->second)->x_grid/GRID_SIZE < (corr->x_grids/2) && ((it)->second)->y_grid/GRID_SIZE < (corr->y_grids/2) ){
            ((it)->second)->adaptive_block = 0;
        } else if ( ((it)->second)->x_grid/GRID_SIZE >= (corr->x_grids/2) && ((it)->second)->y_grid/GRID_SIZE  < (corr->y_grids/2)){
            ((it)->second)->adaptive_block = 1;
        } else if ( ((it)->second)->x_grid/GRID_SIZE < (corr->x_grids/2) && ((it)->second)->y_grid/GRID_SIZE >= (corr->y_grids/2)){
            ((it)->second)->adaptive_block = 2;
        } else if ( ((it)->second)->x_grid/GRID_SIZE >= (corr->x_grids)/2 && ((it)->second)->y_grid/GRID_SIZE >= (corr->y_grids/2)){
            ((it)->second)->adaptive_block = 3;
        }
*/

		//cout << "Printing the adaptive Blocks " << endl;
	/*	for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
                                                it != ckt->Gates_Map.end()  ;  ++it){
			//cout << (it->second)->gate_name <<  "   " << (it->second)->adaptive_block << endl;
			if((it->second)->adaptive_block >= 4){
				cout << "There is some problem with the above printed gate " << endl;
				//int wait;
				//cin >> wait;
			}
		} */
		
		//int wait;
		//cin >> wait;


   ckt->determine_boundary_crossing_nets();

    //what the hell man what i am going to do??

}

void Timing_Analysis::Monte_Carlo(int no_of_iterations){

		static int init = 0;

		if (init == 0){

    corr->random_numbers       = (double**)malloc(corr->x_grids*corr->y_grids*sizeof(double*));
    corr->random_numbers_width = (double**)malloc(corr->x_grids*corr->y_grids*sizeof(double*));

    for(int allocate=0; allocate < corr->x_grids*corr->y_grids; allocate++){
        corr->random_numbers[allocate]       = gaussian_random_numbers(no_of_iterations);
        corr->random_numbers_width[allocate] = gaussian_random_numbers(no_of_iterations);
    }

    //now perform the cholesky decomposition on the correlation matrix
    // Again I need to assign more memory... :(
    double **chol;
    chol = (double**)malloc(corr->x_grids*corr->y_grids*sizeof(double*));

    
    for(int allocate=0; allocate < corr->x_grids*corr->y_grids; allocate++){
        //*(chol+allocate) = (double*)malloc(corr->x_grids*corr->y_grids*sizeof(double));
        chol[allocate] = (double*)malloc(corr->x_grids*corr->y_grids*sizeof(double));
        for(int init=0; init < corr->x_grids*corr->y_grids; init++){    
            chol[allocate][init] = 0.0;
        }
    }

    // This one for length
    cholesky(corr->corr_matrix, corr->x_grids*corr->y_grids, NULL, 0, chol, NULL, 0 );
    matrix_multiplication(corr->random_numbers, chol, corr->x_grids*corr->y_grids, no_of_iterations);
    
    // This one for width
    cholesky(corr->corr_matrix_width, corr->x_grids*corr->y_grids, NULL, 0, chol, NULL, 0 );
    matrix_multiplication(corr->random_numbers_width, chol, corr->x_grids*corr->y_grids, no_of_iterations);
		init = 1;
}
}
/*
    for(int f1 = 0; f1 < corr->x_grids*corr->y_grids; f1++){
        for(int f2 = 0; f2 < corr->x_grids*corr->y_grids; f2++){
            if(corr->eigen_values[f2] < 0){
                //corr->eigen_values[f2] = -corr->eigen_values[f2];
                corr->eigen_values[f2] = 0;
            }
            //cout << "Eigen vector = " << corr->eigen_vectors[f1][f2]  << "Eigen Value = " << corr->eigen_values[f2] << endl;
            chol[f1][f2] = corr->eigen_vectors[f1][f2] * sqrt(corr->eigen_values[f2]);
        }
    }

    //int ww;
    //cin >> ww;

    //transpose of matrix
    double te;
    for(int f1 = 0; f1 < corr->x_grids*corr->y_grids; f1++){
        for(int f2 = f1; f2 < corr->x_grids*corr->y_grids; f2++){
            te = chol[f1][f2];
            chol[f1][f2] = chol[f2][f1];
            chol[f2][f1] = te;
        }
    }
    
    for(int h=0; h< corr->x_grids*corr->y_grids; h++ ){
       cout << corr->random_numbers[0][h] << "  ";
    }
    cout << endl;

    for(int h=0; h< corr->x_grids*corr->y_grids; h++ ){
        cout << corr->random_numbers[1][h] << "  ";
    }
    cout << endl;

    cout << endl << endl;;

    matrix_multiplication(corr->random_numbers, chol, corr->x_grids*corr->y_grids, no_of_iterations);

    for(int h=0; h< corr->x_grids*corr->y_grids; h++ ){
       cout << corr->random_numbers[0][h] << "  ";
    }
    cout << endl;

    for(int h=0; h< corr->x_grids*corr->y_grids; h++ ){
       cout << corr->random_numbers[1][h] << "  ";
    }
    cout << endl;
   
    int wai;
    cin >> wai;



    for(int f1 = 0; f1 < corr->x_grids*corr->y_grids; f1++){
        for(int f2 = 0; f2 < corr->x_grids*corr->y_grids; f2++){
            if(corr->eigen_values_width[f2] < 0){
                //corr->eigen_values_width[f2] = -corr->eigen_values_width[f2];
                corr->eigen_values_width[f2] = 0;
            }
            chol[f1][f2] = corr->eigen_vectors_width[f1][f2] * sqrt(corr->eigen_values_width[f2]);
        }
    }

    //transpose of matrix
    //double te;
    for(int f1 = 0; f1 < corr->x_grids*corr->y_grids; f1++){
        for(int f2 = f1; f2 < corr->x_grids*corr->y_grids; f2++){
            te = chol[f1][f2];
            chol[f1][f2] = chol[f2][f1];
            chol[f2][f1] = te;
        }
    }

    matrix_multiplication(corr->random_numbers_width, chol, corr->x_grids*corr->y_grids, no_of_iterations); */
/*
    for(int p=0; p < corr->x_grids*corr->y_grids ;p++){
        for(int q = 0; q < corr->x_grids*corr->y_grids ; q++){
                cout << chol[p][q] << "\t";
        }
        cout << endl;
    }

    int wait34;
    cin >> wait34;
*/
    
/*
		FILE *fp1;
		fp1 = fopen("test.txt","w");
		double temp1;
		for(int i=0; i<no_of_iterations; i++){
			//cout << corr->random_numbers[0][i] << "    " << corr->random_numbers_width[0][i]  << endl;
			//if(20.1 + corr->random_numbers[1][i] + corr->random_numbers[0][i] > 20.1 + corr->random_numbers[1][i]+ 1*corr->random_numbers_width[0][i] ){
			if(corr->random_numbers[0][i] > corr->random_numbers[0][i] ){
				temp1 = 10.1 + 1*corr->random_numbers[0][i] + corr->random_numbers[0][i];
			} else {
				temp1 = 10.1 + 1*corr->random_numbers[0][i] + corr->random_numbers[0][i];
			}
			//temp1 = 10 + 2*corr->random_numbers[0][i];
			fprintf(fp1,"%f\n",temp1);
		}

		fclose(fp1);


		//int wait12;
		//cin >> wait12;
*/
  /*  
    adaptive_blocks = new int*[ckt->no_of_adaptive_blocks];

    for(int g=0; g < ckt->no_of_adaptive_blocks; g++){
        adaptive_blocks[g] = new int[no_of_iterations];
        srand(g);
        for(int t=0; t < no_of_iterations;t++){
            //adaptive_blocks[g][t] = rand();
            adaptive_blocks[g][t] = 0;
        }
    } 
  }
*/
     //   ckt->Monte_Carlo1(corr, no_of_iterations, false);
    //This is the normal timing analysis for calculating the arrival time!!
    //ckt->pert_for_timing(-1, corr, adaptive_blocks); // Here -1 specify that 
                                                     // we need to only calculate the arrival time and the slack !!

    //After we have have generated the circuit we will apply the PERT algorithm to calculate the timings	
    //for monte carlo method we iterate in a loop
    /*
    FILE *fp;
    fp = fopen("delay.txt","w");
		double delay;
    //for(int i=0; i<no_of_iterations; i++){
        //ckt->pert_for_timing(i, corr, adaptive_blocks);
        //delay = ckt->print_time_at_outputs();
				//cout << "Delay = " << delay;
				//fprintf(fp,"%f\n",delay);
    clock_t t1,t2;
    t1=clock();

        ckt->Monte_Carlo1(corr, no_of_iterations, false);
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"Time Elapsed(s) for Monte Carlo = "<<diff/CLOCKS_PER_SEC<<endl;
    //}*/

//}

double Timing_Analysis::DSTA(){
		FILE *fp;
    //ckt->Critical_Path_Delay(corr,1);
    //int wait; 
    //cin >> wait;
    return ckt->Critical_Path_Delay(corr,0); // Here -1 specify that 
		//ckt->print_time_at_outputs(fp);
}

void Timing_Analysis::SSTA(map<string, pair<double, double> > &mu_sigma, map<string, pair<double, double> > &arrival_time_slack, map<string, vector<float> >  gate_parameters, double required_output_time, map<int, bool> adapt_block_exist){

		//ckt->update_gate_delays(gate_parameters, required_output_time);
    double slack;


    only_FBB = true;
		req_time = required_output_time;

		if(only_FBB == true){
		/*  for(int i=0; i < ckt->no_of_adaptive_blocks; i++){
				ckt->prob_low[i] = 0;
        if(adapt_block_exist[i] == true){
					ckt->high_low[i] = true;
				}
			} */
			//ckt->update_gate_delays(gate_parameters, required_output_time);
    	ckt->SSTA(corr,corr->x_grids*corr->y_grids);
			//ckt->calculate_required_arrival_time(required_output_time);
    	/*for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
      	                                          it != ckt->Gates_Map.end()  ;  ++it){
        	mu_sigma[it->first] = it->second->output_mu_sigma;	
  	      slack = 0 - it->second->output_mu_sigma.first - 3*(it->second->output_mu_sigma.second) + it->second->required_arrival_time;	
      	  arrival_time_slack[it->first] = make_pair(it->second->required_arrival_time, slack);
    	}*/
      //ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);
		//int wait1;
		//cin >> wait1;
			return;
		}


		int recalculate = 0;
		gate_class* gate;
    vector<int> remaining_blocks;
    remaining_blocks.clear();

    for(int i=0; i < ckt->no_of_adaptive_blocks; i++){
        ckt->high_low[i] = false;
        if(adapt_block_exist[i] == true){
            remaining_blocks.push_back(i);
						//cout << "the adaptivity is with " << i << endl;
						//int waot;
						//cin >> waot;
        }
    }

    clock_t t1,t2;
    t1=clock();

    ckt->SSTA(corr,corr->x_grids*corr->y_grids);
		ckt->calculate_required_arrival_time(required_output_time);
    t2=clock();
    float diff ((float)t2-(float)t1);
    //cout<<"Time Elapsed(s) for SSTA = "<<diff/CLOCKS_PER_SEC<<endl;

    for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
                                                it != ckt->Gates_Map.end()  ;  ++it){
        mu_sigma[it->first] = it->second->output_mu_sigma;

        slack = 0 - it->second->output_mu_sigma.first - 3*(it->second->output_mu_sigma.second) + it->second->required_arrival_time;
/*
				if(slack <= 0){
					if(adapt_block_exist[it->second->adaptive_block] == true){
						ckt->high_low[it->second->adaptive_block] = true;
						recalculate = 1;
					}
				}
*/
				
        arrival_time_slack[it->first] = make_pair(it->second->required_arrival_time, slack);
    }

    ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);

		cout << "DONE WITH PROBABILITY CALCULATIONS " << endl;
		//int wait1;
		//cin >> wait1;

    double max_prob = 0.01;
    int block_no;
	
	int exist[16];

    for(int i=0; i < ckt->no_of_adaptive_blocks; i++){
        if( max_prob < ckt->prob_high[i]){
            max_prob = ckt->prob_high[i];
            block_no = i;
            recalculate = 1;
        }
    }

    if(recalculate == 1){
				for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
					if(*it == block_no ){
        		remaining_blocks.erase(it);
						break;
					}
				}
    }




		while(recalculate == 1){

			//cout << "Recalculate Blocks " << remaining_blocks.size() << "   " << block_no  <<endl;
			//int wait1;
			//cin >> wait1;

      ckt->high_low[block_no] = true;

      /******************************************************************/
      /* Necessary Steps */
      /******************************************************************/
			ckt->update_gate_delays(gate_parameters, required_output_time);
    	ckt->SSTA(corr,corr->x_grids*corr->y_grids);
			ckt->calculate_required_arrival_time(required_output_time);
    	for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
      	                                          it != ckt->Gates_Map.end()  ;  ++it){
        	mu_sigma[it->first] = it->second->output_mu_sigma;	
  	      slack = 0 - it->second->output_mu_sigma.first - 3*(it->second->output_mu_sigma.second) + it->second->required_arrival_time;	
      	  arrival_time_slack[it->first] = make_pair(it->second->required_arrival_time, slack);
    	}
      /***********************************************************************/
      
      ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);
      recalculate = 0;
      max_prob = 0.02;
	/*		cout << ckt->prob_high[0] << "   " << ckt->prob_low[0] << endl;
			cout << ckt->prob_high[1] << "   " << ckt->prob_low[1] << endl;
			cout << ckt->prob_high[2] << "   " << ckt->prob_low[2] << endl;
			cout << ckt->prob_high[3] << "   " << ckt->prob_low[3] << endl;
*/
			//ist = 0;
			for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
        if( max_prob < ckt->prob_high[*it] ){
            max_prob = ckt->prob_high[*it];
            block_no = *it;
						//cout << "BLOCK NO " << block_no << endl;
            recalculate = 1;
        }
			/*	if( block_no == *it ){
					exist[*it] = 1;
				} else {
					exist[*it] = 0;
				} */
			}
      //for(int i=0; i < ckt->no_of_adaptive_blocks; i++){
			/*	
        if( max_prob < ckt->prob_high[i] && exist[i] == 0 ){
            max_prob = ckt->prob_high[i];
            block_no = i;
						cout << "BLOCK NO " << block_no << endl;
            recalculate = 1;
        } */
      //}
      if(recalculate == 1){
				for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
					//cout << *it << endl;
					if(*it == block_no ){
        		remaining_blocks.erase(it);
						//cout << "Erasing " << *it << endl;
						break;
					}
				}
      }

		}

		//cout << "Assigned high FBB " << endl;
		//int wait;
		//cin >> wait;

    //Now we have calculted the probabilities of going high 
    //So we calculate the probablity of going low 
    //We will exclude those blocks which have gone high
    //Another thing to note is that all now there is no block which should
    //go high !!
    if(remaining_blocks.size() > 0){
        recalculate = 1;
    } else {
        recalculate = 0;
    }

/*
		cout << "recalculqate " << recalculate << endl;
			cout << ckt->prob_high[0] << "   " << ckt->prob_low[0] << endl;
			cout << ckt->prob_high[1] << "   " << ckt->prob_low[1] << endl;
			cout << ckt->prob_high[2] << "   " << ckt->prob_low[2] << endl;
			cout << ckt->prob_high[3] << "   " << ckt->prob_low[3] << endl;

		int t;
		cin >> t;
*/
    double min = 1;
    bool reverse_decision = false;
    while(recalculate == 1){
        recalculate = 0;
        min = 1;
        reverse_decision = false;
        for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
            // Turn it to reverese BB and if timing meets then let it be else return it to original state
            // turn to reverse BB will be on its probability
            if(ckt->prob_low[*it] < min){
            	min = ckt->prob_low[*it];
            	recalculate = 1;
            	block_no = *it;
            }
        }

        if(recalculate == 1){
            ckt->high_low[block_no] = true;
    
          /******************************************************************/
          /* Necessary Steps */
          /******************************************************************/
			    ckt->update_gate_delays(gate_parameters, required_output_time);
    	    ckt->SSTA(corr,corr->x_grids*corr->y_grids);
			    ckt->calculate_required_arrival_time(required_output_time);
    	    for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
      	                                              it != ckt->Gates_Map.end()  ;  ++it){
            mu_sigma[it->first] = it->second->output_mu_sigma;	
  	        slack = 0 - it->second->output_mu_sigma.first - 3*(it->second->output_mu_sigma.second) + it->second->required_arrival_time;	
      	    arrival_time_slack[it->first] = make_pair(it->second->required_arrival_time, slack);
            //if(slack < 0){
            //    reverse_decision = true;
            //}
    	    }
          ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);
					for(int j=0; j< ckt->no_of_adaptive_blocks; j++){
						if(ckt->prob_high[j] > 0){
							reverse_decision = true;
						}
					}
/*			cout << ckt->prob_high[0] << "   " << ckt->prob_low[0] << endl;
			cout << ckt->prob_high[1] << "   " << ckt->prob_low[1] << endl;
			cout << ckt->prob_high[2] << "   " << ckt->prob_low[2] << endl;
			cout << ckt->prob_high[3] << "   " << ckt->prob_low[3] << endl;  */
          /***********************************************************************/

        }

		/*		cout << "Reverse Decision " << reverse_decision << endl;
				int g;
				cin >> g;  */

        if(reverse_decision == true){
            ckt->high_low[block_no] = false;
			      ckt->update_gate_delays(gate_parameters, required_output_time);
    	      ckt->SSTA(corr,corr->x_grids*corr->y_grids);
			      ckt->calculate_required_arrival_time(required_output_time);
    	      for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
      	                                                  it != ckt->Gates_Map.end()  ;  ++it){
                mu_sigma[it->first] = it->second->output_mu_sigma;	
  	            slack = 0 - it->second->output_mu_sigma.first - 3*(it->second->output_mu_sigma.second) + it->second->required_arrival_time;	
      	        arrival_time_slack[it->first] = make_pair(it->second->required_arrival_time, slack);
    	       }
            ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);
        }
			//exist = 0;
			for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
				if( block_no == *it ){
				//	exist = 1;
				}
			}
        for(std::vector<int>::iterator it = remaining_blocks.begin(); it != remaining_blocks.end(); ++it){
					if(*it == block_no ){
        		remaining_blocks.erase(it);
						break;
					}
				}
    }



}

void Timing_Analysis::get_pca_params_plac_info( map<string, vector<double> > &pca_params, map<string, vector<float> > &plac_info ){

		std::vector<double> k_params;
		int s = 0;
		int size;
    for(std::map<string, gate_class*>::iterator it  = ckt->Gates_Map.begin();
                                                it != ckt->Gates_Map.end()  ;  ++it){
        int row, column;
				double res, cap;
				size = it->second->k_parameters.size();
				s = 0;
				for(std::vector<double>::iterator it1  = it->second->k_parameters.begin();
																					it1 != it->second->k_parameters.end()  ; ++it1){
					res = it->second->gate_resistance; cap = it->second->load_capacitance;
					if(s< size/2 ){
						k_params.push_back( (*it1)*3/0.15/res/cap );
					} else {
						k_params.push_back( (*it1)*3/0.08/res/cap );
					}
					s++;
				}
        pca_params[it->first] = k_params;
				k_params.clear();

	//cout << "Printing the k parameters while returning to optimization tools " << endl;
	//for(vector<double>::iterator it1 = pca_params[it->first].begin(); it1 != pca_params[it->first].end(); ++it1){
	//for(vector<double>::iterator it1 = it->second->k_parameters.begin(); it1 != it->second->k_parameters.end(); ++it1){
	//	cout << *it1 << endl;
	//}

        row    = it->second->y_grid/GRID_SIZE;
        column = it->second->x_grid/GRID_SIZE;
        plac_info[it->first].push_back(row);
        plac_info[it->first].push_back(column);
        plac_info[it->first].push_back(it->second->adaptive_block);
    }
}

void Timing_Analysis::get_probability(map<int,vector<float> > &probability ){
		for(int i = 0; i< ckt->no_of_adaptive_blocks; i++){
			ckt->high_low[i] = false;	
		}

    	      ckt->SSTA(corr,corr->x_grids*corr->y_grids);
			      ckt->calculate_required_arrival_time(req_time);
            ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks, only_FBB);
		
    //ckt->calculate_probability_low_high(ckt->no_of_adaptive_blocks);
		for(int i = 0; i < ckt->no_of_adaptive_blocks; i++){
			if(ckt->prob_high[i] > 0){
    		probability[i].push_back(1 - (float)ckt->prob_high[i]);	 //ZERO
    		probability[i].push_back((float)ckt->prob_high[i]);	     //FBB
    		probability[i].push_back(0);	                           //RBB
				//prob_type[i] = "FBB";
			} else if(ckt->prob_low[i] > 0){
    		probability[i].push_back(1 - (float)ckt->prob_low[i]);	 //ZERO
    		probability[i].push_back(0);														 //FBB
    		probability[i].push_back((float)ckt->prob_low[i]);				 //RBB
    		//probability[i] = ckt->prob_low[i];	
				//prob_type[i] = "RBB";
			}	else {
    		probability[i].push_back(1);
    		probability[i].push_back(0);	
    		probability[i].push_back(0);	
    		//probability[i] = 0;	
				//prob_type[i] = "0";
			}
		}
/*
    cout << "Printing Porbabilities " << endl;
    cout << probability[0] << "   " << prob_type[0] <<endl;
    cout << probability[1] << "   " << prob_type[1] << endl;
    cout << probability[2] << "   " << prob_type[2] << endl;
    cout << probability[3] << "   " << prob_type[2] << endl;

    int wait;
    cin >> wait;
  */  
}

//pair<double, double> Timing_Analysis::get_power(map<string, vector<float> > gate_parameters, map<int, bool> adapt_block_exist){
//    int iterations;
//    iterations = 10000;
//
//		Monte_Carlo(iterations);
//
//		ckt->adapt_block_exist = new bool[adapt_block_exist.size()];
//
//		for(int h=0; h <  ckt->no_of_adaptive_blocks; h++){
//			ckt->adapt_block_exist[h] = adapt_block_exist[h];
//		}
//
//    ckt->ran = new double*[ckt->no_of_adaptive_blocks];
//
//		for(int i=0; i<ckt->no_of_adaptive_blocks; i++){
//			ckt->ran[i] = new double[iterations];
//		}
//    // Now I will generate random numbers for probability between 0 and 1
//		for(int i=0; i<ckt->no_of_adaptive_blocks; i++){
//    	srand(2+i);
//    	for(int g=0; g<iterations; g++){
//      	  ckt->ran[i][g] = (float)rand()/RAND_MAX;
//    	}
//		}
///*
//    srand(234);
//    for(int g=0; g<iterations; g++){
//        ckt->ran[1][g] = (float)rand()/RAND_MAX;
//    }
//
//    srand(933);
//    for(int g=0; g<iterations; g++){
//        ckt->ran[2][g] = (float)rand()/RAND_MAX;
//    }
//
//    srand(19589);
//    for(int g=0; g<iterations; g++){
//        ckt->ran[3][g] = (float)rand()/RAND_MAX;
//    }
//*/
//    return ckt->Calculate_Power_Monte_Carlo(iterations, corr, gate_parameters);
//
//}

void Timing_Analysis::get_timing_yield(){
    //Calculate the number of runs required
    int runs;
    int blocks_to_evaluate;
    int bit0, bit1, bit2, bit3;
    int num;
    int blocks[4];
    double prob = 0;
    double prob2 = 0;
    double threshold;
    //runs = pow(2,no_of_adaptive_blocks);
    //
    map<int, double> probability;

    //get probabilities
    //get_probability(probability);

    //remove the cases like 0 or totally 1 cases by removing that no of blocks
    blocks_to_evaluate = ckt->no_of_adaptive_blocks;
    int j = 0;
    blocks[0] = 5; blocks[1] = 5; blocks[2] = 5; blocks[3] = 5; 
    threshold = 550;

    double prob_w;
    for(int i=0; i<4; i++){
        if(ckt->prob_high[i] < 0.01 || ckt->prob_high[i] > 0.99){ // these means they are absolutley going to be in either high or low modes
            blocks_to_evaluate--;
            if(ckt->prob_high[i] < 0.01){
                ckt->high_low[i] = false;
            } else{
                ckt->high_low[i] = true;
            }
        } else {
            blocks[j] = i;
            j++;
        }
    }

    runs = pow(2.0,(double)blocks_to_evaluate);
    cout << "runs = " << runs << endl;

    // this one will get timing yield for SSTA

    clock_t t1,t2;
    t1=clock();
    for(int i=0; i<runs; i++){
        num = i;
        bit0 = num%2; num = num/2;
        bit1 = num%2; num = num/2;
        bit2 = num%2; num = num/2;
        bit3 = num%2; num = num/2;

        prob_w = 1;

            if(blocks[0] <= 4){
                ckt->high_low[blocks[0]] = bit0 == 1 ? true : false;
                if(bit0 == true){
                    prob_w *= ckt->prob_high[blocks[0]];
                } else {
                    prob_w *= 1 - ckt->prob_high[blocks[0]];
                }
            }
            if(blocks[1] <= 4){
                ckt->high_low[blocks[1]] = bit1 == 1 ? true : false;
                if(bit1 == true){
                    prob_w *= ckt->prob_high[blocks[1]];
                } else {
                    prob_w *= 1 - ckt->prob_high[blocks[1]];
                }
            }
            if(blocks[2] <= 4){
                ckt->high_low[blocks[2]] = bit2 == 1 ? true : false;
                if(bit2 == true){
                    prob_w *= ckt->prob_high[blocks[2]];
                } else {
                    prob_w *= 1 - ckt->prob_high[blocks[2]];
                }
            }
            if(blocks[3] <= 4){
                ckt->high_low[blocks[3]] = bit3 == 1 ? true : false;
                if(bit3 == true){
                    prob_w *= ckt->prob_high[blocks[3]];
                } else {
                    prob_w *= 1 - ckt->prob_high[blocks[3]];
                }
            }

            ckt->SSTA(corr,corr->x_grids*corr->y_grids);
            prob = ckt->get_max_delay_prob(threshold);
            prob2 += prob*prob_w;
            //cout << "Probability = " << prob << endl;
            cout << "Probability 2  = " << prob2 << "   " <<  prob_w <<endl;
    //cout << "Printing working high and lows " << endl;
    //cout << ckt->high_low[3] << "  ";
    //cout << ckt->high_low[2] << "  ";
    //cout << ckt->high_low[1] << "  ";
    //cout << ckt->high_low[0] << endl;
    //int wait;
    //cin >> wait;
    // Now i know the probability for the delay of a circuit. I just have to multiply with the 
    // probability of the circuit in which it is working.

        
    }

    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"Time Elapsed(s) for Timing Yield SSTA = "<<diff/CLOCKS_PER_SEC<<endl;


    // this one will get timing yield for Monte Carlo
    // Generate random numbers based upon probabilities;

    int iterations;
    iterations = 20000;

    ckt->ran = new double*[4];

    ckt->ran[0] = new double[iterations];
    ckt->ran[1] = new double[iterations];
    ckt->ran[2] = new double[iterations];
    ckt->ran[3] = new double[iterations];

    // Now I will generate random numbers for probability between 0 and 1
    srand(1);
    for(int g=0; g<iterations; g++){
        ckt->ran[0][g] = (float)rand()/RAND_MAX;
    }

    srand(234);
    for(int g=0; g<iterations; g++){
        ckt->ran[1][g] = (float)rand()/RAND_MAX;
    }

    srand(933);
    for(int g=0; g<iterations; g++){
        ckt->ran[2][g] = (float)rand()/RAND_MAX;
    }

    srand(19589);
    for(int g=0; g<iterations; g++){
        ckt->ran[3][g] = (float)rand()/RAND_MAX;
    }

    //clock_t t1,t2;
    t1=clock();
    //ckt->Monte_Carlo1(corr, iterations, true);
    t2=clock();
    float diff1 ((float)t2-(float)t1);
    cout<<"Time Elapsed(s) for Timing Yield Monte Carlo = "<<diff1/CLOCKS_PER_SEC<<endl;



}

float wq;
double Timing_Analysis::get_timing_yield_temp(map<int, bool> adapt_block_exist, double threshold, map<string, vector<float> > gate_parameters ){

vector<int> blocks;
int a_block;



for(map<string, gate_class*>::iterator it = ckt->Gates_Map.begin(); it != ckt->Gates_Map.end(); ++it){
    a_block = it->second->adaptive_block;
    if(find(blocks.begin(), blocks.end(), a_block) == blocks.end()){
        if(ckt->prob_high[a_block] == 1){
            ckt->high[a_block] = 1;
            ckt->low[a_block] = 0;
        } else if(ckt->prob_low[a_block] == 1) {
            ckt->high[a_block] = 0;
            ckt->low[a_block] = 1;
        } else if((ckt->prob_low[a_block] + ckt->prob_high[a_block]) == 0) {
            ckt->high[a_block] = 0;
            ckt->low[a_block] = 0;
        } else {
            blocks.push_back(a_block);
            //cout << a_block << "  " << ckt->prob_low[a_block] << "  " << ckt->prob_high[a_block] << endl;
        }
    }
}

/*

//int wait;
//cin >> wait;
wq = 1;
for(int g=0; g < ckt->no_of_adaptive_blocks; g++){
    ckt->high_low[g] = false;
    //ckt->prob_high[g] = 0;
    //ckt->prob_low[g] = 0;
    wq *= 1 - (ckt->prob_high[g] + ckt->prob_low[g]);
}

ckt->SSTA(corr, corr->x_grids*corr->y_grids);
wq = ckt->get_max_delay_prob(threshold)*wq;

clock_t t11,t21;
t11=clock();
    Monte_Carlo(20000);
    ckt->Monte_Carlo1(corr, 20000, true, ckt->no_of_adaptive_blocks, threshold);
t21=clock();
float diff11 ((float)t21-(float)t11);
cout<<"Time Elapsed for Monte Carlo ==== "<<diff11/CLOCKS_PER_SEC<<endl;
//int wai1t;
//cin >> wai1t;
getchar();

ckt->Calculate_RAT(corr->x_grids*corr->y_grids, threshold);
// Monte Carlo is working fine !!
//Monte_Carlo(2000);
//ckt->Monte_Carlo1(corr, 2000, false);
//ckt->Monte_Carlo_RAT(corr, 2000, false, threshold);
//exit(1);
cout << "THE number of NODES originally are " << ckt->Gates_Map.size() << endl;
ckt->Prune_nodes();
cout << "THE number of NODES left are " << ckt->Gates_Map.size() << endl;
//ckt->compress_gates();
cout << "THE number of NODES left after compression are " << ckt->Gates_Map.size() << endl;

//int wait;
//cin >> wait;

clock_t t1,t2;
t1=clock();

third_method(blocks, 2, threshold);

t2=clock();
float diff1 ((float)t2-(float)t1);
cout<<"Time Elapsed(s) ==== "<<diff1/CLOCKS_PER_SEC<<endl;


*/


//Monte_Carlo(1000);
//ckt->Monte_Carlo_RAT(corr, 1000, false, threshold);

	//vector<int> blocks;

clock_t t11, t21;
t11 = clock();
Monte_Carlo(1000);
ckt->Monte_Carlo1(corr, 1000, true, ckt->no_of_adaptive_blocks, threshold);
t21 = clock();
float diff11((float)t21 - (float)t11);
cout << "Time Elapsed for Monte Carlo ==== " << diff11 / CLOCKS_PER_SEC << endl;




	vector<int>::iterator it;
	blocks.clear();
	int runs;
	int blocks_to_evaluate = 0;
	for(int i=0; i < ckt->no_of_adaptive_blocks; i++){
		if(adapt_block_exist[i] == true){
			ckt->high_low[i] = true;
			blocks.push_back(i);
			blocks_to_evaluate++;
		} else {
			ckt->high_low[i] = false;
		}
	}


	 double prob_w;
		//int g = blocks.size();
		double prob_w2;

    for(int i=0; i< blocks.size(); i++){
        if( (ckt->prob_high[blocks[i]] < 0.01 || ckt->prob_high[blocks[i]] > 0.99) && ckt->prob_low[blocks[i]] == 0  || (ckt->prob_low[blocks[i]] < 0.01 || ckt->prob_low[blocks[i]] > 0.99) && ckt->prob_high[blocks[i]] == 0){ // these means they are absolutley going to be in either high or low modes
            blocks_to_evaluate--;
            if(ckt->prob_high[blocks[i]] < 0.01 && ckt->prob_low[blocks[i]] == 0){
                ckt->high_low[blocks[i]] = false;
								it = find(blocks.begin(), blocks.end(), blocks[i]);
								if(it != blocks.end()) {blocks.erase(it); i--;}
            } else if(ckt->prob_high[blocks[i]] > 0.99 && ckt->prob_low[blocks[i]] == 0){
                ckt->high_low[blocks[i]] = true;
								it = find(blocks.begin(), blocks.end(), blocks[i]);
								if(it != blocks.end()) {blocks.erase(it); i--;}
            } else if(ckt->prob_low[blocks[i]] < 0.01 && ckt->prob_high[blocks[i]] == 0){
								ckt->high_low[blocks[i]] = false;
								it = find(blocks.begin(), blocks.end(), blocks[i]);
								if(it != blocks.end()) {blocks.erase(it); i--;}
            } else if(ckt->prob_low[blocks[i]] > 0.99 && ckt->prob_high[blocks[i]] == 0){
								ckt->high_low[blocks[i]] = true;
								it = find(blocks.begin(), blocks.end(), blocks[i]);
								if(it != blocks.end()) { blocks.erase(it); i--;}
            } else {
								ckt->high_low[blocks[i]] = false;
						}
        } 
    }

    runs = pow(2.0,(double)blocks_to_evaluate);
    cout << "runs = " << runs << endl;

	t11 = clock();
		int t;
		int bit;
		//double prob_w;
		double prob;
		double prob2 = 0;
		for(int i=0;i<runs;i++){
			t = i;
			prob_w = 1;
			prob_w2 = 1;
			for(int j=0; j < blocks_to_evaluate; j++){
				bit = t%2;
				t = t/2;
				if(bit == 0){
					ckt->high_low[blocks[j]] = false;
					prob_w *= 1 - ckt->prob_high[blocks[j]] - ckt->prob_low[blocks[j]];
				} else {
					ckt->high_low[blocks[j]] = true;
					prob_w *= ckt->prob_high[blocks[j]] + ckt->prob_low[blocks[j]];
				}
			}
			//prob
			ckt->update_gate_delays(gate_parameters, threshold);
			ckt->SSTA(corr,corr->x_grids*corr->y_grids);
      prob = ckt->get_max_delay_prob(threshold);
	  cout << prob << " " << prob_w << endl;
	  getchar();
      prob2 += prob*prob_w;
		}


		t21 = clock();
		diff11= (float)t21 - (float)t11;

		cout << "Time Elapsed for PCA ==== " << diff11 / CLOCKS_PER_SEC << endl;
  //return ckt->get_max_delay_prob(threshold);
  return prob2;


/////Temprarily return 0////

//return 0.0;
}


void Timing_Analysis::third_method(vector<int> &blocks, int no_of_levels, double threshold ){
    //vector<int> blocks;
    //int no_of_levels = 4;
    //blocks.push_back(0);
    //blocks.push_back(1);
    //blocks.push_back(2);
    //blocks.push_back(3);

    bitset<64> case1;
    bitset<64> temp_case;

    int no_of_bits;

    // The way i have to do is to select specified no of blocks and then evaluate all the possible
    // cases for them
    no_of_bits = blocks.size();

    map<int, vector< bitset<64> > > cases;
    map<int, vector< bitset<64> > >::iterator it_map;
    vector< bitset<64> >::iterator it_vec;
    int start_point;

    // So outer loop will be to select the no of blocks and the internal loop will evaluate all the possible no of cases !!
    for(int i=1; i <= no_of_bits; i++){ // no of blocks to select !!
        // depending upon the no of blocks that many bits will be set
        for(int z=0; z < i; z++){
            case1[z] = 1;
        }

        cases[i].push_back(case1);
        start_point = cases[i].size() - 1;

        for(int j=1; j <= i; j++){ // this loop will select the blocks
            // No of blocks to select are nCp terms.
            // The way to evaluate this is to add 1 to those positions which have zero infront of them one by one.
            // Now creating a linked list
            cout << "HELLO " << no_of_bits << endl;
            for(int t = 0; t < cases[i].size(); t++){
                for(int b=0; b < no_of_bits; b++){
                    if(b+1 < no_of_bits){
                        temp_case = cases[i][t];
                        if(temp_case[b+1] == 0 && temp_case[b] == 1 ){
                            temp_case[b+1] = 1;
                            temp_case[b]   = 0;
                            it_vec = find(cases[i].begin(), cases[i].end(), temp_case);
                            if( it_vec == cases[i].end() ){
                                cases[i].push_back(temp_case);
                            }
                        }
                    }
                }
            }
        }
    }
    /*for(it_map = cases.begin(); it_map != cases.end(); ++it_map){
        for(it_vec = (it_map->second).begin(); it_vec != (it_map->second).end(); ++it_vec){
            cout << *it_vec << "\t";
        }
        cout << endl;
    }*/
    //call the method which calculate all the possible case's for the given no of blocks
    calculate_all_possible_cases(cases, no_of_levels, blocks.size(), blocks, threshold);
}

void Timing_Analysis::calculate_all_possible_cases(map<int, vector<bitset<64> > >&cases, int no_of_levels, int no_blocks, vector<int>&act_blocks, double threshold){
//void calculate_all_possible_cases(){
    int no_of_blocks;
    //int no_of_levels;
    int runs;
    int previous_runs;
    map<int, vector< bitset<64>  > >::iterator it_map;
    vector< bitset<64> >::iterator it_vec;
    double prob = 0;
    double prob2 = 0;
    double w_prob = 1;
    int count = 0;
    //int *blocks;
    //no_of_blocks = 3;
    //no_of_levels = 3;
    int *blocks;
    int *blocks_new;
    int ex;

clock_t t1,t2;
float t3;
t3 = 0.0;
    // require storage of cases to be excluded !!
    vector<int*> exclude_cases;
    
    blocks = new int[64]();

    previous_runs = 0;
    runs = 0;

    for(it_map = cases.begin(); it_map != cases.end(); ++it_map){
        for(it_vec = (it_map->second).begin(); it_vec != (it_map->second).end(); ++it_vec){
            no_of_blocks = it_map->first;
            //previous_runs += runs;
            runs = pow(no_of_levels - 1, (double)no_of_blocks);// - no_of_blocks*(pow(no_of_levels, no_of_blocks - 1) - 1); 
            //if(no_of_blocks > 1){
            //    runs--;
            //}

            for(int i=0; i < no_of_blocks; i++){
                blocks[i] = 1;
            }
            
            for(int i=0; i < runs; i++){
                for(int j=0; j < no_of_blocks; j++){
                    if(blocks[j] >= no_of_levels){
                        blocks[j] = 1;
                        blocks[j+1]++;
                    }
                }
                //cout << i+1 << ")   " << blocks[2] << "  " << blocks[1] << "   " << blocks[0] << endl;
                count++;
                //Before calling the SSTA here we will check that is this case is to be excluded !!
                ex = 0;
                blocks_new = new int[64]();
                int k;
                k=0;
                w_prob = 1;
                for(int a=0; a<64; a++ ){
                    if((*it_vec)[a] == 1){
                        blocks_new[a] = blocks[k];
                        ckt->high_low[act_blocks[a]] = true;
                        if(blocks[k] == 1){
                            ckt->high[act_blocks[a]] = 0;
                            ckt->low[act_blocks[a]] = 1;
                            if(a < no_blocks){
                                w_prob *= ckt->prob_low[act_blocks[a]];
                            //cout << ckt->prob_high[act_blocks[a]] << "  "  << act_blocks[a] << endl;
                            }
                        } else if(blocks[k] == 2){
                            ckt->high[act_blocks[a]] = 1;
                            ckt->low[act_blocks[a]] = 0;
                            if(a < no_blocks){
                                w_prob *= ckt->prob_high[act_blocks[a]];
                            //cout << ckt->prob_low[act_blocks[a]] << "   " << act_blocks[a] <<  endl;
                            }
                        } else {
                            ckt->high[act_blocks[a]] = 0;
                            ckt->low[act_blocks[a]] = 0;
                            if(a < no_blocks){
                                w_prob *= 1 - (ckt->prob_high[act_blocks[a]] + ckt->prob_low[act_blocks[a]]);
                            }
                        }
                        k++;
                    } else {
                        blocks_new[a] = 0;
                            ckt->high[act_blocks[a]] = 0;
                            ckt->low[act_blocks[a]] = 0;
                            if(a < no_blocks){
                                w_prob *= 1 - (ckt->prob_high[act_blocks[a]] + ckt->prob_low[act_blocks[a]]);
                                //cout << ckt->prob_high[act_blocks[a]] << "  "  << act_blocks[a] << endl;
                                //cout << ckt->prob_low[act_blocks[a]] << "   " << act_blocks[a] <<  endl;
                            }
                    }
                //cout << blocks_new[a];
                }
                //cout << endl;
                for(std::vector<int*>::iterator it = exclude_cases.begin(); it != exclude_cases.end(); ++it ){
                    //cout << endl;
                  ex = 0;  //cout << "Existing CASES " << exclude_cases.size() << endl; 
                    ex = 1;
                    for(int h=0; h < no_blocks; h++){
                        //cout << (*it)[h];
                        if((*it)[h] <= blocks_new[h] ){
                            ex *= 1;
                        } else {
                            ex *= 0;
                        }
                    } 
                //cout << endl;
                    if(ex == 1){
                        prob2 += 1*w_prob;
                        cout << "1 less :) " << endl << endl;
                        break;
                    }
                }

                if(ex == 0){
                    //cout << "YES" << endl;
                    //prob += ckt->get_max_delay_prob(threshold);
                    prob = 0;
                    if(w_prob != 0){
t1=clock();
                        ckt->SSTA(corr, corr->x_grids*corr->y_grids);
                        prob = ckt->get_max_delay_prob(threshold);
                        prob2 += w_prob*prob;
t2=clock();
float diff1 ((float)t2-(float)t1);
t3 += diff1;
//cout<<"Time Elapsed(s) ==== "<<diff1/CLOCKS_PER_SEC<<endl;
                        cout << "Probability = " << prob2 << "   "  << prob << endl;
                    }
                //call the ssta here :)
                    if(prob >= 0.997){
                        //store this case for exclusion !!
                        exclude_cases.push_back(blocks_new);
                    } else {
                        delete blocks_new;
                    }
                } else {
                    delete blocks_new;
                }
                blocks[0]++;
            }
        }
    }
    cout << "COUNT = " << count << endl;
    cout << "TIME = " << t3/CLOCKS_PER_SEC << endl;
    cout << "PROB = " << prob2 + wq << endl;
}

//void Timing_Analysis::evaluate_timing_yield(){ // this will evaluate all the possible cases and also prune the cases !!
//}


#endif
