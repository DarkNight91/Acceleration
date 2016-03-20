#ifndef CORRELATION_MODEL_CPP
#define CORRELATION_MODEL_CPP

#include "correlation_model.h"
#include "colours.h"
#include "stdlib.h"
#include "fstream"
#include "string.h"
#include "map"
#include "vector"

//I willbe using only the placement model which is generated through the CAPO tool !!
//
//Reading the placement file out.pl as provided in the data
/*int* correlation_model::parse_placement_file(const char *placement_file_name, circuit_class *ckt){
    
    ifstream placement_file;
    char* word1;
    char* c_string;
    gate_class *gate_name;

    string a_line_from_file;

    int* temp_max;

    temp_max = (int *)malloc(2*sizeof(int));
    temp_max[0] = 0;
    temp_max[1] = 1;

    c_string = (char*)malloc(100);
    word1    = (char*)malloc(10);

    placement_file.open(placement_file_name);

    if(!placement_file.is_open()){
        printf("The placement file can't be opened!! \n");
        return NULL;
    }else{
        printf(KBLU "The placement file is opened!! \n" KNRM);
    }
    
    while(getline(placement_file,a_line_from_file)){
        strcpy(c_string,a_line_from_file.c_str());
        
        word1 = strtok(c_string, " :\t");
        
#ifdef DEBUG
        //printf("The word read from the placement file is %s \n", word1);
#endif
        if(word1 == NULL){
            continue;
        }

        if(ckt->Gates_Map.find(word1) == ckt->Gates_Map.end()){
            continue;
        }
        gate_name         = (ckt->Gates_Map.find(word1))->second;
        gate_name->x_grid = atof(strtok(NULL, " :\t"));
        gate_name->y_grid = atof(strtok(NULL, " :\t"));

        if(temp_max[0] < gate_name->x_grid){
            temp_max[0] = gate_name->x_grid;

        } else if(temp_max[1] < gate_name->y_grid){
            temp_max[1] = gate_name->y_grid;

        }
    }
    
    return temp_max;
}*/

void correlation_model::generate_correlation_matrix_quadtree(int x_grids, int y_grids, const float *corr_levels, int no_of_levels, string param){
    //since we know the x and y location of the gates we can simple deduce the location of the gate in the 
    //grid. So there is not any need to specifically store the gates grid as 0,1 etc.

    //For generating the correlation model I will be using the quad-tree model as it as used in the earlier tool
    //and the PCA paper.
    //
    //And the data for the correaltion between the gates comes from the same tool

    //Now there are two ways to generate the correlation matrix. On is the by the number in which the grid is divied and the
    //second one in which the number of levels. I think the grid level is more robust since this will model the correlation very well
    //In case after dividing the grid we can ask for the number of levels and the extra number of levels we can say that they have zero 
    //correlation.... Such things we can look later on.

    //Now generate a correlation matrix
    //1) The correlation matrix will be equal to the number of grids

    MatrixXd corr_matrix_eigen(x_grids*y_grids,x_grids*y_grids);
    corr_matrix1 = corr_matrix_eigen;
    double **itrtr;

    itrtr = (double **)malloc(x_grids*y_grids*sizeof(double*));

    for(int i=0; i<x_grids*y_grids; i++){
        *(itrtr+i) = (double *)malloc(x_grids*y_grids*sizeof(double));
    }

    if(param.compare("length") == 0){
        corr_matrix = itrtr;
    } else if(param.compare("width") == 0){
        corr_matrix_width = itrtr;
    } else if(param.compare("oxide") == 0){
        corr_matrix_oxide = itrtr;
    } else {
        cout << "The parameters should be only length, width and oxide. This one is not supported... Exiting ";
        exit(1);
    }

    //intializing the array to all zeros
    for(int n1=0; n1 < x_grids*y_grids; n1++){
        for(int n2=0; n2 < x_grids*y_grids; n2++){
            *(*(itrtr+n1)+n2) = 0;
        }
    }
    //this->print_matrix(x_grids, y_grids);

#ifdef DEBUG
    cout << endl;
    cout << "The size of the array is " << sizeof(corr_matrix) << endl;
    cout << endl;
#endif

    int rows      = -1; // This is the number of the Random variable from 0 to n^2 - 1
    int columns   = -1; 
    int rows_1    = -1; // This is the x_coordinate number or rows of the matrix of correlations
    int columns_1 = -1; // This is the y_coordinate number or columns of the matrix of correlations

    int row_diff;
    int column_diff;
    for(int n1=0; n1 < x_grids*y_grids; n1++){
        columns_1++;
        if((n1)%x_grids == 0){
            rows_1++;
            columns_1 = 0;
        }
        rows    = -1;
        columns = -1;

#ifdef DEBUG
        //printf("The word from n1 file is %d \n", n1);
#endif
    
        for(int n2=0; n2 < x_grids*y_grids; n2++){
            //i++;
            columns++;
            if((n2)%x_grids == 0){
                rows++;   
                columns = 0;
            }
            
            row_diff    = rows    - rows_1;
            column_diff = columns - columns_1;
    
            for(int t=1; t < no_of_levels+1; t++){
                if(row_diff == 0 && column_diff == 0){
                    *(*(itrtr+n1)+n2) = 1;
                    corr_matrix_eigen(n1,n2) = 1;

                //} else if( ((row_diff > 0 && column_diff >= 0) || (column_diff > 0 && row_diff >= 0)) 
                } else if(  ((int)(rows/(t*2)) == (int)(rows_1/(t*2))) 
                            && ((int)(columns/(t*2)) == (int)(columns_1/(t*2)))){
                    *(*(itrtr+n1)+n2) += *(corr_levels + t - 1);
                    corr_matrix_eigen(n1,n2) = *(corr_levels + t - 1);
                    //cout << rows << " " << rows_1 << " " << columns << " " << columns_1 << endl;
                    //cout << *(*(itrtr+n1)+n2) << endl;

                }
            } // End of t For loop
        } // End of n2 For Loop 
    }
        //Printing the generated co-related matrix
#ifdef DEBUG
    //this->print_matrix(x_grids, y_grids);
#endif

}

void correlation_model::generate_correlation_matrix_grid(int x_grids, int y_grids, const float *corr_levels, int no_of_levels, string param){
    //since we know the x and y location of the gates we can simple deduce the location of the gate in the 
    //grid. So there is not any need to specifically store the gates grid as 0,1 etc.

    //For generating the correlation model I will be using the grid model
    //
    //And the data for the correaltion between the gates comes from the quad-tree model with the same way it is used in the quad-tree
    //model as for a particular grid we sum all the correlation levels depending upon the distance between them

    //Now generate a correlation matrix
    //1) The correlation matrix will be equal to the number of grids
    //2) The data for the correlation will be as explained earlier.
    //3) The model however will be based on the grid model
    //4) The gird model is better than the quad tree model as it does not have the deficiency of the quad tree
    //as the two adjacent blocks will have different correlations
    //

/*
    MatrixXd matx(3,3);
    matx(0,0) = 1;
    matx(0,1) = 0.6;
    matx(0,2) = 0.3;
    matx(1,0) = 0.6;
    matx(1,1) = 1;
    matx(1,2) = 0.5;
    matx(2,0) = 0.3;
    matx(2,1) = 0.5;
    matx(2,2) = 1;

    EigenSolver<MatrixXd> es;
    //MatrixXf A = MatrixXf::Random(4,4);
    es.compute(matx, true);
    int ww;

    cout << es.eigenvalues() << endl << endl;
    cout << es.eigenvectors() << endl << endl;

    cin >> ww;


    for(int i=0; i< 3; i++){
        //eigen_values_temp[i] = es.eigenvalues()(i).real();
        for(int j=0; j< 3; j++){
            cout << "\t" <<  es.eigenvectors()(i,j).real();
        }
        cout << endl;
    }

    cin >> ww;

*/

    MatrixXd corr_matrix_eigen(x_grids*y_grids,x_grids*y_grids);
    corr_matrix1 = corr_matrix_eigen;
    //corr_matrix1 = MatrixXd::zero(x_grids*y_grids,x_grids*y_grids);

    double **itrtr;

    itrtr = (double **)malloc(x_grids*y_grids*sizeof(double*));

    for(int i=0; i<x_grids*y_grids; i++){
        *(itrtr+i) = (double *)malloc(x_grids*y_grids*sizeof(double));
    }

    if(param.compare("length") == 0){
        corr_matrix = itrtr;
    } else if(param.compare("width") == 0){
        corr_matrix_width = itrtr;
    } else if(param.compare("oxide") == 0){
        corr_matrix_oxide = itrtr;
    } else {
        cout << "The parameters should be only length, width and oxide. This one is not supported... Exiting ";
        exit(1);
    }
    

    //intializing the array to all zeros
    for(int n1=0; n1 < x_grids*y_grids; n1++){
        for(int n2=0; n2 < x_grids*y_grids; n2++){
            *(*(itrtr+n1)+n2) = 0;
            corr_matrix1(n1,n2) = 0;
        }
    }

    int rows      = -1; // This is the number of the Random variable from 0 to n^2 - 1
    int columns   = -1; 
    int rows_1    = -1; // This is the x_coordinate number or rows of the matrix of correlations
    int columns_1 = -1; // This is the y_coordinate number or columns of the matrix of correlations

    int row_diff;
    int column_diff;
    for(int n1=0; n1 < x_grids*y_grids; n1++){
        //i++;
        columns_1++;
        if((n1)%x_grids == 0){
            rows_1++;
            columns_1 = 0;

        }
        rows    = -1;
        columns = -1;

#ifdef DEBUG
        //printf("The word from n1 file is %d \n", n1);
#endif
    
        for(int n2=0; n2 < x_grids*y_grids; n2++){
            //i++;
            columns++;
            if((n2)%x_grids == 0){
                rows++;   
                columns = 0;

            }
            
            row_diff    = rows    - rows_1;
            column_diff = columns - columns_1;
    
						int t;
						t = no_of_levels+1;
            //for(int t=1; t < no_of_levels+1; t++){
                if(row_diff == 0 && column_diff == 0){
                    *(*(itrtr+n1)+n2) = 1;
                    corr_matrix1(n1,n2) = 1;

                //} else if( abs(row_diff) < pow(2,(double)t) && abs(column_diff) < pow(2,(double)t) ){
                } else if( abs(row_diff) <= t && abs(column_diff) <= t ){
										int a;
										a = abs(row_diff) > abs(column_diff) ? abs(row_diff) : abs(column_diff) ;
                    //*(*(itrtr+n1)+n2) += *(corr_levels + t - 1);
                    //corr_matrix1(n1,n2) += *(corr_levels + t - 1);
                    *(*(itrtr+n1)+n2) = *(corr_levels + a - 1);
                    corr_matrix1(n1,n2) = *(corr_levels + a - 1);
                    //corr_matrix1(n1,n2) = 0.9;  // <- debugging
                } else {
                    *(*(itrtr+n1)+n2) = 0;
                    corr_matrix1(n1,n2) = 0;
								}
            //} // End of t For loop
        } // End of n2 For Loop 
    }
        //Printing the generated co-related matrix
//#ifdef DEBUG
    //this->print_matrix();
    //int wait; cin >> wait;
//#endif

 this->get_eigen_values(param);

}


void correlation_model::print_matrix(){
    double **itrtr;
    cout<<"Prinitng the Correlation Matrix" << endl;
    itrtr = this->corr_matrix;
    for(int p1=0; p1 < x_grids*y_grids; p1++){
        for(int p2=0; p2 < x_grids*y_grids; p2++){
            cout << *(*(itrtr+p1)+p2) << "\t" ; //p1 << " " << p2 << endl;
        }
        cout<<endl<<endl;
    }
}

void correlation_model::get_eigen_values(string param){
    double *eigen_values_temp;
    double **eigen_vectors_temp;

    eigen_values_temp = new double[x_grids*y_grids];

    eigen_vectors_temp = new double*[x_grids*y_grids];

    for(int a=0; a < x_grids*y_grids; a++){
        eigen_vectors_temp[a] = new double[x_grids*y_grids];
    }


    if(param.compare("length") == 0){
        eigen_values = eigen_values_temp;
        eigen_vectors = eigen_vectors_temp;
    } else if(param.compare("width") == 0){
        eigen_values_width = eigen_values_temp;
        eigen_vectors_width = eigen_vectors_temp;
    } else if(param.compare("oxide") == 0){
        eigen_values_oxide = eigen_values_temp;
        eigen_vectors_oxide = eigen_vectors_temp;
    } else {
        cout << "The parameters should be only length, width and oxide. This one is not supported... Exiting ";
        exit(1);
    }

    EigenSolver<MatrixXd> es;
    //MatrixXf A = MatrixXf::Random(4,4);
    es.compute(corr_matrix1, /* computeEigenvectors = */ true);
    for(int i=0; i<x_grids*y_grids; i++){
        eigen_values_temp[i] = es.eigenvalues()(i).real();
        for(int j=0; j<x_grids*y_grids; j++){
            eigen_vectors_temp[i][j] = es.eigenvectors()(i,j).real();
        }
    }

    cout << "The eigenvalues of A are: " << es.eigenvalues().transpose() << endl;
    //es.compute(A + MatrixXf::Identity(4,4), false); // re-use es to compute eigenvalues of A+I
    //cout << "The eigenvalues of A+I are: " << es.eigenvalues()(0).real() << endl;
    //cout << "The eigenvalues of A+I are: " << es.eigenvectors() << endl;
    //cout << "The eigenvalues of A+I are: " << es.eigenvalues().cols() << endl;
    
}

#endif
