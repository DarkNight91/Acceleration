#ifndef ISPD_PARSE_CPP
#define ISPD_PARSE_CPP

#include "ispd_parse.h"

void ispd_parse_verilog(string filename, circuit_class *ckt){

		cout << endl << filename << endl;
    VerilogParser vp (filename);
    gate_class* gate;
		bool valid;

    string modulename;
    valid = vp.read_module(modulename);
    assert(valid);

    do{
        string primaryInput;
        valid = vp.read_primary_input (primaryInput);

        if(valid){
						if(primaryInput != "ispd_clk"){
            	ckt->add_gate_to_map(primaryInput, "PI");
						}
            //gate = new[gate_class];
            //ckt->Gates_Map[primaryInput] = gate;
            //gate->gate_type = "PI";
            //gate->gate_name = primatyInput;
        }

    } while(valid);

    do{
        string primaryOutput;
        valid = vp.read_primary_output(primaryOutput);
    
        if(valid){
            ckt->add_gate_to_map(primaryOutput, "PO");
            //gate = new[gate_class];
            //ckt->Gates_Map[priamryOutput] = gate;
            //gate->gate_type = "PO";
            //gate->gate_name = primaryInput;
        }

    }while(valid);

		do {
    string net ;
    valid = vp.read_wire (net) ;

  /*  if (valid)
      cout << "Net: " << net << endl ;*/

  } while (valid) ;


    do {
    string cellType, cellInst ;
    vector<std::pair<string, string> > pinNetPairs ;
    
    valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

  /*  if (valid) {
      cout << cellType << " " << cellInst << " " ;
      for (int i=0; i < pinNetPairs.size(); ++i) {
        cout << "(" << pinNetPairs[i].first << " " << pinNetPairs[i].second << ") " ;
      }
      cout << endl ;
    } */

    string gate_name;
    string in_not_clk;
    int loc;
    int flag;
    int flag2;
    flag2 = 0;
    if(valid){
        loc = pinNetPairs.size() - 1;
        gate_name = pinNetPairs[loc].second;
        
        flag = 0;
        for(int i=0; i < loc; ++i){
            if(flag == 1){
						    ckt->add_gate_to_map(pinNetPairs[i].second, "dummy");
						    ckt->add_gates_input_output(gate_name, pinNetPairs[i].second);
                ckt->add_gate_to_map(pinNetPairs[i].second, "DFF_I");
            }

            if(pinNetPairs[i].second == "ispd_clk"){
                flag2 = 1;
                ckt->add_gate_to_map(gate_name, "DFF");
                if(i != 0){
						    		ckt->add_gate_to_map(in_not_clk, "dummy");
						    		ckt->add_gates_input_output(gate_name, in_not_clk);
                    ckt->add_gate_to_map(in_not_clk, "DFF_I");
                } else {
                    flag = 1;
                }
            }

            in_not_clk = pinNetPairs[i].second;
        }

        // This is for normal circuit !!
        if(flag2 == 0){
            ckt->add_gate_to_map(gate_name, cellType);
            for(int i=0; i < loc; ++i){
						    ckt->add_gate_to_map(pinNetPairs[i].second, "dummy");
						    ckt->add_gates_input_output(gate_name, pinNetPairs[i].second);
            }
        }
    }
    
  } while (valid) ;

}

void ispd_parse_lib(string filename, map<string, double> &library){ // I need cell area form this file

    LibParser lp(filename);

    LibParserCellInfo cell;
    bool valid;

    do{
        valid = lp.read_cell_info(cell);

        if(valid){
            library[cell.name] = cell.area;
						//cout << cell.name << "   " << library[cell.name] << endl;
        }
    } while(valid);
}

void generate_placement(circuit_class *ckt, string benchmark, map<string, double> &library){
    //write nets file
    string filename;
    filename = "output/" + benchmark + ".nets";

		//cout << "In Placement Function " << endl;

    ofstream fp;
    fp.open(filename.c_str(), ios::out);
    write_nets_file(fp, ckt);
    fp.close();

    filename = "output/" + benchmark + ".nodes";
    fp.open(filename.c_str(), ios::out);
    write_nodes_file(fp, ckt, library);
    fp.close();

    filename = "output/" + benchmark + ".wts";
    fp.open(filename.c_str(), ios::out);
    write_wts_file(1, ckt, fp); // 1  for wts and 0 for pl files
    fp.close();
}

void write_nets_file(ofstream& fp, circuit_class *ckt){
		//cout << "IN nets file " << endl;

    fp << "UCLA nets 1.0\n";
    fp << "NumNets :      " << ckt->Gates_Map.size() << "\n";

    std::map<string, gate_class*>::iterator it;    
    std::vector<gate_class*>::iterator it_vec;

    int sum=0;
    for(it = ckt->Gates_Map.begin(); it != ckt->Gates_Map.end(); ++it){
        if( it->second->gate_type != "PI"){ // || it->second->gate_type != "PO"){
            sum += it->second->inputs.size() + 1;
        }
    }

    fp << "NumPins :      " << sum << "\n";

    for(it = ckt->Gates_Map.begin(); it != ckt->Gates_Map.end(); ++it){
        if((it->second)->gate_type != "PO"){
            fp << "NetDegree : " << it->second->outputs.size() + 1 <<"\n";
        		for(it_vec = it->second->outputs.begin(); it_vec != it->second->outputs.end(); ++it_vec){
            	fp << "     " << (*it_vec)->gate_name << "\t\t" << "I" << " : 0.5 0.5\n";
        		}
           		fp << "     " << (it->second)->gate_name << "\t\t" << "O" << " : 0.5 0.5\n";
        }
    }

    for(it_vec = ckt->Primary_Outputs.begin(); it_vec != ckt->Primary_Outputs.end(); ++it_vec){
			if((*it_vec)->gate_type == "PO"){
        fp << "NetDegree : " << 2 <<"\n";
        fp << "     " << "PO"+(*it_vec)->gate_name << "\t\t" << "I" << " : 0.5 0.5\n";
        fp << "     " << (*it_vec)->gate_name << "\t\t" << "O" << " : 0.5 0.5\n";
			}
    }
}

void write_nodes_file(ofstream& fp, circuit_class *ckt, map<string, double> &library){
    fp << "UCLA nodes    1.0\n";

		int sum = 0;
		int sum_out = 0;
		for(vector<gate_class*>::iterator it = ckt->Primary_Inputs.begin(); it != ckt->Primary_Inputs.end(); ++it){
			if((*it)->gate_type == "PI"){
				sum++;
			}
		}

		for(vector<gate_class*>::iterator it = ckt->Primary_Outputs.begin(); it != ckt->Primary_Outputs.end(); ++it){
			if((*it)->gate_type == "PO"){
				sum_out++;
			}
		}

    fp << "NumNodes :\t" << ckt->Gates_Map.size() + sum_out <<"\n";
    fp << "NumTerminals :\t"<< sum + sum_out << "\n";

    std::map<string, gate_class*>::iterator it_map;    
    std::vector<gate_class*>::iterator it_vec;

		string abc;
    for(it_map = ckt->Gates_Map.begin(); it_map != ckt->Gates_Map.end(); ++it_map){
        if( it_map->second->gate_type_add != "PI" ){
					if(it_map->second->gate_type_add  == "DFF"){
						abc = "ms00f80";
					} else {
						abc = it_map->second->gate_type_add;
					}
            fp << "\t" << it_map->second->gate_name << "\t" << ceil(sqrt(library[abc]));
            fp << "\t" << ceil(sqrt(library[abc])) << endl;
        }
    }
    
    for(it_vec = ckt->Primary_Inputs.begin();  it_vec != ckt->Primary_Inputs.end(); ++it_vec){
			if((*it_vec)->gate_type_add == "PI") {
            fp << "\t" << (*it_vec)->gate_name << "\t" << 2;
            fp << "\t" << 2 << " terminal"<< endl;
			}
    }

    for(it_vec = ckt->Primary_Outputs.begin();  it_vec != ckt->Primary_Outputs.end(); ++it_vec){
			if((*it_vec)->gate_type == "PO") {
            fp << "\t" << "PO"+(*it_vec)->gate_name << "\t" << 2;
            fp << "\t" << 2 << " terminal" <<endl;
			}
    }

}


void write_wts_file(int ext, circuit_class *ckt , ofstream& fp)
{
	if (ext = 0)                           //.pl
	{/*
		ofstream io("output_pl");
		io << "UCLA" << "pl" << "1.0" << endl;
		map<string, pair<int, int> > ::iterator it = PI.begin();
		map<string, pair<int, int> > ::iterator it_O = PO.begin();
		for (; it != PI.end(); ++it)
		{
			io << it->first << "/t" << 0 << "/t" << 0 << ":" << "N" <<endl;
		}

		for (; it != PI.end(); ++it)
		{
			io << it->first << "/t" << it->second.first << "/t" << it->second.second << ":" << "N"<<endl;

		}

		for (; it_O != PO.end(); ++it)
		{
			io << it->first << "/t" << it->second.first << "/t" << it->second.second << ":" << "N"<<endl;

		}
		io.close();
*/
	}
	else if (ext = 1)                   //.wts
	{
		fp << "UCLA " << "wts " << "1.0" << endl;
    std::map<string, gate_class*>::iterator it_map;    
    std::vector<gate_class*>::iterator it_vec;

    for(it_map = ckt->Gates_Map.begin(); it_map != ckt->Gates_Map.end(); ++it_map){
        //if(it->second->gate_type !=  "PI"){
            fp << it_map->second->gate_name << " 1" << "\n";
        //}
    }

		for (it_vec = ckt->Primary_Outputs.begin(); it_vec != ckt->Primary_Outputs.end(); ++it_vec)
		{
			if((*it_vec)->gate_type == "PO"){
				fp << "PO"+(*it_vec)->gate_name << " 1" << endl;
			}
		}
	}

}


#endif

