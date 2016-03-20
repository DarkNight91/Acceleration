#ifndef PARSER_CPP
#define PARSER_CPP

#include "parser.h"
#include "colours.h"
#include "circuit_class.h"

#include"string"
#include"iostream"
#include"fstream"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>

char* process_string(char[]);

using namespace std;

int parse_netlist(circuit_class *ckt, const char *netlist_file_name){

	ifstream netlist_file;
	netlist_file.open(netlist_file_name);

	string netlist_line;
	string rd_wd_fr_line;
	int rd_wd_st_pos;
	int rd_wd_end_pos;
	
	if(!netlist_file.is_open()){ //is_open returns true or false
		printf(KRED "ERROR: The file could not be opened.... Exiting \n" KNRM);
		exit(1);
	} else {
		printf(KBLU "The file is opened....\n" KNRM);
	}

	// Parsing the netlist file by reading line by line
	//The file format is as follow
	//**********************************************************
    // #hello I am dummy file
  	//
 	//  INPUT(G0)
 	//  OUTPUT(G1)
  	//
  	// G1 = NOT(G0)          	
	//**********************************************************
	// TODO: Write the parser in Bison
			char *word;
			char *word1;
			char *word_2;
			char *cstr;
            word   = (char *)malloc(1000*sizeof(char));
            word_2 = (char *)malloc(1000*sizeof(char));
			//cstr = new char(100);
			cstr   = (char *)malloc(1000*sizeof(char));

	while(getline(netlist_file,netlist_line)){
		if(netlist_line[0] == '#'){ // If comment ignore the line
			//printf(KBLU "A comment is found \n" KNRM);
			continue;
		}
		// if it is a Primary Input we will store the gate name in the gate_class and mark it is as a primary input  
		else if(netlist_line.find("INPUT",0) != std::string::npos){
			//printf(KBLU "A input is found \n" KNRM);

			rd_wd_st_pos = netlist_line.find("(");
			rd_wd_end_pos = netlist_line.find(")");
			rd_wd_fr_line = netlist_line.substr(rd_wd_st_pos+1,rd_wd_end_pos-rd_wd_st_pos-1);

			printf(KYEL);
			//cout<<"The input gate name is: " << rd_wd_fr_line <<  endl;
			printf(KNRM);

			ckt->add_gate_to_map(rd_wd_fr_line, "PI");
			//ckt->add_gate_to_map("G0", "PI");
		} else if(netlist_line.find("OUTPUT",0) != std::string::npos){
			//printf(KMAG "A output is found \n" KNRM);

			rd_wd_st_pos = netlist_line.find("(");
			rd_wd_end_pos = netlist_line.find(")");
			rd_wd_fr_line = netlist_line.substr(rd_wd_st_pos+1,rd_wd_end_pos-rd_wd_st_pos-1);

			printf(KYEL);
			//cout<<"The output gate name is: " << rd_wd_fr_line <<  endl;
			printf(KNRM);

			ckt->add_gate_to_map(rd_wd_fr_line, "PO");
		
		} else if (netlist_line[0] != ' ' && netlist_line[0] != '\n' && netlist_line[0] != '\r' && netlist_line.size() != 0){
			//printf(KBLU "A Gate is found \n" KNRM);

            //cout << netlist_line << endl;

			strcpy(cstr,netlist_line.c_str());

            
			word1 = strtok(cstr, " =(,)"); // Read the output gate name
            if(word1 != NULL)
                strcpy(word,process_string(word1));

			if(word1 != NULL){
				word1 = strtok(NULL, " =(,)"); // Read the gate Type
                if(word1 != NULL)
                    strcpy(word_2,process_string(word1));

                if(word1 != NULL) {
                    if (strcmp(word_2, "DFF") == 0){   // This will check whether the given gate is a DFF
                                                       // This is found by looking at the type of gate
                        assert(word_2 != NULL);
                        // The Dff is a special type of circuit 
                        // We have to mark it as a Primary input and also as a primary output
                        ckt->add_gate_to_map(string(word),"PI"); // This is we are adding it to the Gate_Maps
                        word1 = strtok (NULL, " =(,)");     // Read the input of the D Flip Flop
                        if(word1 != NULL)
                            strcpy(word_2,process_string(word1));

                    // This is one thing i might have forgotten. I think the way i have done is to enter the DFF as
                    // two type of gate. One is the output of the DFF as a primary input and the input of the DFF as
                    // DFF which is also equivalent to making it as a primary output !!! <-- TODO: Need to check this properly!!
                        assert(word1 != NULL);
                        ckt->add_gate_to_map(string(word_2),"DFF");
                    } else {
                        ckt->add_gate_to_map(string(word),word_2);
                    }
                }
			} else {
				continue;
			}

            word1 = strtok(NULL, " =(,)"); //After we have read the Gate Type we will read the first input

            if(word1 != NULL)
                strcpy(word_2,process_string(word1));

			printf(KRED);
			while (word1 != NULL)
			{
				if(word1 != NULL){
					if( (word_2[0] >= 65 && word_2[0] <= 90) || (word_2[0] >= 97 && word_2[0] <= 122)){
                        assert(word != NULL);
                        assert(word_2 != NULL);
                        //std::cout << "the input being added to " << word << " is " <<word_2 << endl;
						ckt->add_gate_to_map(string(word_2), "dummy");
						ckt->add_gates_input_output(string(word),string(word_2));
					}
				}
                
			    word1 = strtok(NULL, " =(,)");
                if(word1 != NULL)
                    strcpy(word_2,process_string(word1));
  			} //end while
			printf(KNRM);
		}
	}

    free(word);
    free(word_2);
    free(cstr);
    netlist_file.close();
    return 0;
}


char* process_string(char *word){
    //std::cout<< "The string is "<< word << endl;
    if(word == NULL) return NULL;
    for(unsigned int i = 0;i<strlen(word)+1; i++){
        if((word[i] < 48) && (word[i] != '\0' || word[i] > 57) && (word[i] < 65 || word[i] > 90)){
            //printf("In the process string before the replacement %c \n", word[i]);
            //word[i] = '\0';
            //printf("In the process string after the replacement %c \n", word[i]);
        }
    }
    return word;
}

//I willbe using only the placement model which is generated through the CAPO tool !!
//
//Reading the placement file out.pl as provided in the data
void parse_placement_file(const char *placement_file_name, circuit_class *ckt, float temp_max_min[]){
    
    ifstream placement_file;
    char* word1;
    char* word12;
    char* c_string;
    gate_class *gate_name;

    string a_line_from_file;


    word1 = NULL;
    c_string = NULL;
    //temp_max_min = NULL;

    //temp_max_min    = (float *)malloc(4*sizeof(float));
    temp_max_min[0] = 0;
    temp_max_min[1] = 0;
    temp_max_min[2] = 0;
    temp_max_min[3] = 0;

    c_string = (char*)malloc(1000*sizeof(char));
    if(!c_string){
        cout << "Memory cannot be assigned " << endl;
        exit(1);
    }
    //word1    = (char*)malloc(sizeof(100));
    //word12 = word1;

    placement_file.open(placement_file_name);

    if(!placement_file.is_open()){
        printf("The placement file can't be opened!! \n");
        return;
    }else{
        printf(KBLU "The placement file is opened!! \n" KNRM);
    }
    
    while(getline(placement_file,a_line_from_file)){
        strcpy(c_string,a_line_from_file.c_str());

        //cout << c_string << endl;
        
        word1 = strtok(c_string, " :\t");
        
#ifdef DEBUG
        //printf("The word read from the placement file is %s \n", word1);
#endif
        if(word1 == NULL){
            continue;
        }

        //cout << word1 << "   " << ckt->Gates_Map.size() << endl;

        if(ckt->Gates_Map.find(string(word1)) == ckt->Gates_Map.end()){
            continue;
        }
        gate_name         = (ckt->Gates_Map.find(word1))->second;
        gate_name->x_grid = atof(strtok(NULL, " :\t"));
        gate_name->y_grid = atof(strtok(NULL, " :\t"));

        if(temp_max_min[0] < gate_name->x_grid){
            temp_max_min[0] = gate_name->x_grid;
        }

        if(temp_max_min[1] < gate_name->y_grid){
            temp_max_min[1] = gate_name->y_grid;
        }

        if(temp_max_min[2] > gate_name->x_grid){
            temp_max_min[2] = gate_name->x_grid;
        }
        
        if(temp_max_min[3] > gate_name->y_grid){
            temp_max_min[3] = gate_name->y_grid;
        }


    }
    
    free(c_string);
    //free(word1);

    placement_file.close();
    
//    return temp_max_min;
}

#endif
