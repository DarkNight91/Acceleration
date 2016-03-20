#include "ocl_opt.h"
#include "TimingAnalysis.h"

/*------------------------------------------------------------------------*/
/*                           Macro Definitions                            */
/*------------------------------------------------------------------------*/
#define     FBB_ONLY            (YES)
#define     IS_SIMULTANEOUS     (NO)//YES = simultaneously Adaptivity +Varibility Only, NO = Sensitive method
//#define     IS_ISCAS85_TEST     (YES)//YES = ISCAS85; NO = ISPD Test
#define     IS_NO_ADAPTIVITY    (NO)//next time, to use environment variable
//#define     IS_OPTIMAL_CODE     (YES)//YES = optimal code, NO = original code
#define     IS_EQUAL_LOAD       (NO)//YES = divide the fanout object euqally; NO = HaoHe's method
#define     IS_ADAPT_NUM_MODEL  (NO)//YES = area/power overhead related to the # of gates;

#define     SSTA_PARAM_IO1      (0)//Gate size
#define     SSTA_PARAM_IO2      (1)//Resistance
#define     SSTA_PARAM_IO3      (2)//Capacitance
#define     SSTA_PARAM_IO4      (3)//Vth voltage
#define     SSTA_PARAM_IO5      (4)//Leakage Power
#define     SSTA_PARAM_MAX      (5)

#define     PLAC_PARAM_IO1      (0)//X-position
#define     PLAC_PARAM_IO2      (1)//Y-position
#define     PLAC_PARAM_IO3      (2)//adaptivity_block

#define     SUBPP_AREA_SUM      (0)
#define     SUBPP_POWER_SUM     (1)
#define     SUBPP_LEAKP_SUM     (2)
#define     SUBPP_DELAY_SUM     (3)
#define     SUBPP_AREA_LM       (4)
#define     SUBPP_DELAY_LM      (5)
#define     SUBPP_MAX_DATA_NUM  (6)

#define     MAX_TUNING_UNIT     (4)//current tuning unit
#define     MAX_TUNING_FFF      (0xFFFF)//Invalid tuning unit

#define     ADAPTIVITY_0        (0)//no adaptivity circuit
#if IS_SIMULTANEOUS
#define     ADAPTIVITY_1        (1)//adaptivity circuit exist
#else
#define     ADAPTIVITY_1        (ADAPTIVITY_0)//adaptivity circuit doesn't exist
#endif
#define     MAX_AREA_LIMIT      (1000)//this might change for different circuit
#define     RESIST_PARAM_P      (0.2)//R=P*L/W


void ocl_opt::print_edges(){
	for(int i=0; i<igraph_vector_size(&edges);i++)
	{	
		if(!(i%2))
			cout<<"("<<VECTOR(edges)[i]<<",";
		else
			cout<<VECTOR(edges)[i]<<")  ";
	}
	igraph_vector_t test;
	igraph_vector_init(&test, 0);
	int v = igraph_neighbors(&circuitGraph, &test, 536, IGRAPH_ALL);
	assert(v == 0);
	for(int i=0; i<igraph_vector_size(&test);i++)
	{
			cout<<endl<<VECTOR(test)[i];
	}
}

bool ocl_opt::init_lracc_database(/*CHAR* fpart_name*/)
{
	gGate_Index["in01"] = INV;
    gGate_Index["na02"] = NAND2;
    gGate_Index["na03"] = NAND3;
    gGate_Index["na04"] = NAND4;
    gGate_Index["na05"] = NAND5;
    gGate_Index["na06"] = NAND6;
    gGate_Index["na07"] = NAND7;
    gGate_Index["na08"] = NAND8;
    gGate_Index["na09"] = NAND9;
    gGate_Index["no02"] = NOR2;
    gGate_Index["no03"] = NOR3;
    gGate_Index["no04"] = NOR4;
    gGate_Index["no05"] = NOR5;
    gGate_Index["no06"] = NOR6;
    gGate_Index["no07"] = NOR7;
    gGate_Index["no08"] = NOR8;
    gGate_Index["no09"] = NOR9;
#ifdef IS_ISCAS85_TEST
    gGate_Index["oa22"] = XOR2;
#else
    gGate_Index["oa12"] = OA12;
    gGate_Index["oa22"] = OA22;
    gGate_Index["ao12"] = AO12;
    gGate_Index["ao22"] = AO22;
    gGate_Index["ms00"] = FFMS;
#endif
    return  TRUE;
}

int ocl_opt::find_wire(string& name, vector<wire_par* >& wirePar){
	vector<wire_par* >::iterator iter;
	for(iter=wirePar.begin(); iter!=wirePar.end(); iter++)
	{
		if((**iter).name == name)
			break;
	}
	if(iter!=wirePar.end())
		return iter-wirePar.begin();
	else
		return 0;
}

bool ocl_opt::createEdges (map<string, wire_par* >& wirePar, Gates* gates){
	int n = 0;
	map<string, wire_par* >::iterator itr;
	for(itr=wirePar.begin(); itr!=wirePar.end(); itr++)
	{	
		assert(!((itr->second->out.empty())&&(itr->second->in.empty())));
		if((!itr->second->out.empty())&&(!itr->second->in.empty())){ // primary in/output do not create edges
			for(int i=0; i<itr->second->out.size(); i++)
			for(int j=0; j<itr->second->in.size(); j++)
			{	
				assert(i<1);
				//gates[itr->second->out[i]].is_input = true;
				//gates[itr->second->out[i]].is_output = true;
				VECTOR(edges)[n++] = itr->second->out[i];
				//gates[itr->second->in[j]].is_input = true;
				//gates[itr->second->in[j]].is_output = true;
				VECTOR(edges)[n++] = itr->second->in[j];
			}

		}
		/* PI and PO assign moves to gate_levelize()
		else if(itr->second->in.empty()){
			//assert(!itr->second->out.empty());
			for(int i=0; i<itr->second->out.size();i++){
			gates[itr->second->out[i]].is_output = true;
			cout<<"gate "<<itr->second->out[i]<<"is PO."<<endl;
			//gates[itr->second->out[i]].is_input = false;
			}
		}
		else{
			//assert(!itr->second->in.empty());
			for(int i=0; i<itr->second->in.size();i++){
			cout<<"gate "<<itr->second->in[i]<<"is PI."<<endl;
			//gates[itr->second->in[i]].is_output = false;
			}
		}
		*/
	}

	assert(itr == wirePar.end());
	assert((n%2) == 0);
	assert(n<igraph_vector_size(&edges));
	igraph_vector_resize(&edges,n);    // should check the memory leak, I cannot figure it out now... what is clear now: capacity of edges is num_up, but size is n!!!
	cout<<"edges create:"<<igraph_vector_size(&edges)/2<<endl;
	return true;

}

bool  ocl_opt::load_graph_and_tiles(string  &filename)
{

    COMPONENTTYPE   type;
	map<string, wire_par* > wirePar;
	
	int v;

    {
        VerilogParser vp (filename) ;
        string moduleName ;
        bool valid = vp.read_module (moduleName);
        assert (valid);

        cout << "Module " << moduleName << endl << endl ;

        do {
            string primaryInput ;
            valid = vp.read_primary_input (primaryInput) ;
		//	cout<<primaryInput<<endl;


            /*if (valid) 
            {   //cout << "Primary input: " << primaryInput << endl ;

                // skip the ispd_clk signal, don't create Vertex for it
                if (primaryInput == "ispd_clk")
                    continue;
				wire_par* w = new wire_par;
				w->name = primaryInput;//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
				w->is_primaryIn = true;
				wirePar.push_back(w);
            }*/
        } while (valid) ;

        cout << "Finish parse primaryInput" << endl ;

        do {
            string primaryOutput ;
            valid = vp.read_primary_output (primaryOutput) ;
			/*if (valid) 
            {  
				wire_par* w = new wire_par;
				w->name = primaryOutput;
				w->is_primaryOut = true;
				wirePar.push_back(w);
            }*/
		
        } while (valid) ;

        cout << "Finish parse primaryOutput" << endl ;


		do {
				string net ;
				valid = vp.read_wire (net) ;
		//		if (valid)
		//		cout << "Net: " << net << endl ;

			} while (valid) ;
		 cout << "Finish parse net, do nothing.." << endl ;


        cout << endl ;
        cout << "Cell insts: " << endl ;
		int num_v = 0;
		int num_e_up = 0; //upbound of edges no.

		do {
            string cellType, cellInst ;
            vector<pair<string, string> > pinNetPairs ;
			

            valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

            if (valid) {
				string  footprint;
                footprint.assign(cellType, 0, MAX_FOOTPRINT_LEN);//just get the first 4 char as the footprint

                map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
                H102_ASSERT(gate_lut_iter != gGate_Index.end());
                type = gate_lut_iter->second;
				if(FFMS != type){
					wire_par* w = new wire_par();
					w->name = pinNetPairs.back().second;
					w->out.push_back(num_v);
					wirePar[pinNetPairs.back().second] = w;
					num_v++;
					num_e_up += pinNetPairs.size();
				}
				}

        } while (valid) ;
		cout<<"We need to create "<<num_v<<" gates for this testbench."<<endl;
    
		//initialize Gates and Edges when parsing the cells info
		

	//go through the file again...
	{
        VerilogParser vp (filename) ;
        string moduleName ;
        bool valid = vp.read_module (moduleName) ;
        assert (valid) ;
		//primaryInput
        do {
            string primaryInput ;
            valid = vp.read_primary_input (primaryInput) ;
		} while (valid) ;
		//primaryOutput
        do {
            string primaryOutput ;
            valid = vp.read_primary_output (primaryOutput) ;
        } while (valid) ;
		//net
        do {
            string net ;
            valid = vp.read_wire (net) ;
        } while (valid) ;
	

		gates = new Gates[num_v];
	
		
		v = igraph_vector_init (&edges, num_e_up*2); // <----may use vector_resize_min to deallco, OR we just push back edges one by one//
		assert(v==0);
		
		int count = 0; 
        do {
            string cellType, cellInst ;
            vector<pair<string, string> > pinNetPairs ;
			

            valid = vp.read_cell_inst (cellType, cellInst, pinNetPairs) ;

            if (valid) {
                // set the type ID for the cell
                string  footprint;
                footprint.assign(cellType, 0, MAX_FOOTPRINT_LEN);//just get the first 4 char as the footprint

                map<string,COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(footprint);
                H102_ASSERT(gate_lut_iter != gGate_Index.end());
                type = gate_lut_iter->second;
				
	
			
				/*The gates may become PI or PO after removing the FF, so check back later == */
                if (FFMS != type)
				{
						gates[count].id = count;
						gates[count].type = type; //gates type and gates index saved!
						gates[count].name = pinNetPairs.back().second;
						
						//map<string,wire_par* >::iterator itr;
						for (int i=0; i<pinNetPairs.size()-1; i++)  //*.*
						{	
							if (wirePar.find(pinNetPairs[i].second)!= wirePar.end()){
								wirePar[pinNetPairs[i].second]->in.push_back(count);
							}
							else{
								wire_par* w = new wire_par();
								w->name = pinNetPairs[i].second;
								w->in.push_back(count);
								wirePar[pinNetPairs[i].second] = w;
								//cout<<"successfully read in vertex: "<<count<<endl;  //I only need to go through this once, just set num_e_up as a very large one
							}
						}

						// this is to help use the output pin name as the gate's name :)
						assert(pinNetPairs.back().first == "o");
						Gate_Id_Name_Map[count] = pinNetPairs.back().second;
					count++;
				}
            }
			
        } while (valid) ;
		cout<<count<<" "<<num_v<<endl;
		assert(count==num_v);

    }

    cout << endl ;
    cout << "Success! Topology file "<<filename<<" parsed!"<<endl;

	valid = createEdges(wirePar,gates);
	assert(valid);
	cout << "Success to create edges, now generate graph..."<<endl;
	//cout<<igraph_vector_size(&edges)<<endl;
	//cout<<igraph_vector_capacity(&edges)<<endl;
	igraph_create(&circuitGraph, &edges, 0, 1); // check if valid
	cout<<"Graph generates successfully!"<<endl;
	}


	int no_of_gates;
	no_of_gates = igraph_vcount(&circuitGraph);
	this->no_of_gates = no_of_gates;
	


	/*generate adj list, init level -1, PI,PO assign*/
	int no_of_edges;
	no_of_edges = igraph_vector_size(&edges);
	this->no_of_edges = no_of_edges;
	adj_edges = (int*)malloc(no_of_edges*sizeof(int));

	/*
		Use topological sort
	*/
	int tmp = 1;
	igraph_vector_init(&tmp_nodes, 0);
	tmp = igraph_topological_sorting(&circuitGraph, &tmp_nodes, IGRAPH_OUT);
	assert(!tmp);

	///////////////////////////////////////
	cout<<"Generate adjlist......"<<endl;
	adjlist_generate(gates,no_of_gates);
	cout<<"adjlist generation success!"<<endl;

	/*PI,PO assignment, set level 0, else -1*/
	for(int i=0; i<no_of_gates; i++)
	{
		/*In fact, this gate has a PI, but not level 0*/
		gates[i].gate_capacitance = 1.0;
		gates[i].gate_resistance = 1.0;
		gates[i].level = -1;
		//assert(gates[i].no_of_in||gates[i].no_of_out);
		if(gates[i].no_of_in == 0)
		{
			gates[i].is_input = true;
			gates[i].level = 0;
		}
		else{
			gates[i].is_input = false;
		}

		if(gates[i].no_of_out == 0)
		{
			gates[i].is_output = true;
		}
		else{
			gates[i].is_output = false;
		}
	}
	
	cout<<"PI,PO,LEVEL 0 Success!"<<endl;

	// well, map gate and gate's name here ~~~~
	// be aware that the value is the pointer to a gate
	for (int i = 0; i < no_of_gates; i++)
	{
		Gates_Map[Gate_Id_Name_Map[gates[i].id]] = &gates[i];
	}

	return true;
}

/*levelize the gates and quicksort them by level and .... update the adjlist accordingly*/
int ocl_opt::gates_levelize(Gates* gates, int no_of_gates/*, int* adj_edges*/){
	//gates_levelized = (Gates*)malloc(no_of_gates*sizeof(Gates));

	/*levelize*/
	cout<<"Levelizing each gate......"<<endl;
	for(int i=0; i<no_of_gates; i++)
	{
		if(gates[i].level!=-1)
			continue;
		//cout<<"This is the "<<i<<"'s gate."<<endl;
		gates[i].level = find_level(gates[i],gates/*,adj_edges*/);
		
	}	


	/*Quicksort by level*/
	cout<<"Sort gates by level........."<<endl;
	QuickSortbyLevel(gates,no_of_gates);


	

	/* use tmp index to find the index of older gate's index in new gates_levelized, O(n^2)....*/
	cout<<"Generate tmp index........."<<endl;
	int* index;
	index = (int*)malloc(no_of_gates*sizeof(int));
	for (int i=0;i<no_of_gates;i++)
		for(int j=0;j<no_of_gates;j++)
		{
			if(gates[j].id == i){
				//cout<<i<<endl<<j<<endl;
				//getchar();
				index[i] = j;
				break;
			}
		}


	/*reorder the adjlist*/
	cout<<"Reorder the adjlist........."<<endl;
	adjlist_regenerate(gates,no_of_gates,index);

	cout<<"Levelization completes & adjlist reordered!"<<endl;
	free(index);

	for (int i = 0; i < no_of_gates; i++)
	{
		if (i>0)
			assert((gates[i - 1].level == gates[i].level) || (gates[i - 1].level == (gates[i].level-1)));
	}

	return 0;
}

int ocl_opt::find_level(Gates gate, Gates* gates/*,int* adj_edges*/)
{
	assert(gate.level==-1);
	int level = -1;
	assert(gate.no_of_in!=0);
	
	for(int i=gate.start_in;i<gate.start_in+gate.no_of_in;i++)
	{
		if(gates[adj_edges[i]].level==-1){
			gates[adj_edges[i]].level = find_level(gates[adj_edges[i]],gates/*,adj_edges*/);
			if(gates[adj_edges[i]].level>level)
				level = gates[adj_edges[i]].level;
		}
		else
			if(gates[adj_edges[i]].level>level)
				level = gates[adj_edges[i]].level;
	}
	return level+1;
}


void ocl_opt::QuickSortbyLevel(Gates* gates, int no_of_gates){

	Gates gate_tmp;
	int i=0, j=no_of_gates-1;
	gate_tmp = gates[0];
	int key = gates[0].level*1000 + gates[0].no_of_in;
	if(no_of_gates>1){
	while(i<j)
	{
		for(;j>i;j--)
			if(gates[j].level*1000 + gates[j].no_of_in<key)
			{
				gates[i++]=gates[j];
				break;
			}
		for(;i<j;i++)
			if(gates[i].level*1000 + gates[i].no_of_in>key)
			{
				gates[j--] = gates[i];
				break;
			}
	}
	gates[i] = gate_tmp;
	//cout<<"key: "<<key<<endl;
	//cout<<"i: "<<i<<endl;
	QuickSortbyLevel(gates,i);
	QuickSortbyLevel(gates+i+1,no_of_gates-i-1);
	}
}

void ocl_opt::adjlist_regenerate(Gates* gate, int no_of_gates, int* index){

	igraph_vector_t nei_gates_in;
	igraph_vector_init(&nei_gates_in,0);
	igraph_vector_t nei_gates_out;
	igraph_vector_init(&nei_gates_out,0);
	int counter=0;
	for(int i=0; i<no_of_gates; i++)
	{	
		/*Generate adjlist*/
		int v = igraph_neighbors(&circuitGraph, &nei_gates_in, gates[i].id, IGRAPH_IN);
		assert(!v);
		v = igraph_neighbors(&circuitGraph, &nei_gates_out, gates[i].id, IGRAPH_OUT);
		assert(!v);
		int n1 = igraph_vector_size(&nei_gates_in);
		int n2 = igraph_vector_size(&nei_gates_out);
		gates[i].no_of_in = n1;
		gates[i].no_of_out = n2;
		gates[i].start_in = counter;
		
		for(int j=0; j<n1;j++)
		{
			int m = VECTOR(nei_gates_in)[j];
			adj_edges[counter++] = index[m];
		}
		gates[i].start_out = counter;
		for(int k=0; k<n2;k++)
		{
			int n = VECTOR(nei_gates_out)[k];
			adj_edges[counter++] = index[n];
		}
	}
	igraph_vector_destroy (&nei_gates_in);
	igraph_vector_destroy (&nei_gates_out);

}

void ocl_opt::adjlist_generate(Gates* gate, int no_of_gates){

	igraph_vector_t nei_gates_in;
	igraph_vector_init(&nei_gates_in,0);
	igraph_vector_t nei_gates_out;
	igraph_vector_init(&nei_gates_out,0);
	int counter=0;
	for(int ii=0; ii<no_of_gates; ii++)
	{	
		/*Unnecessasry step?*/
		int i;
		i = VECTOR(tmp_nodes)[ii];
		/*Generate adjlist*/
		int v = igraph_neighbors(&circuitGraph, &nei_gates_in, gates[i].id, IGRAPH_IN);
		assert(!v);
		v = igraph_neighbors(&circuitGraph, &nei_gates_out, gates[i].id, IGRAPH_OUT);
		assert(!v);
		int n1 = igraph_vector_size(&nei_gates_in);
		int n2 = igraph_vector_size(&nei_gates_out);
		gates[i].no_of_in = n1;
		gates[i].no_of_out = n2;
		gates[i].start_in = counter;
		
		for(int j=0; j<n1;j++)
		{

			adj_edges[counter++] = VECTOR(nei_gates_in)[j];
		}
		gates[i].start_out = counter;
		for(int k=0; k<n2;k++)
		{
			adj_edges[counter++] = VECTOR(nei_gates_out)[k];
		}
	}
	igraph_vector_destroy (&nei_gates_in);
	igraph_vector_destroy (&nei_gates_out);

}

void ocl_opt::init_timing(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, string place_bench, double time_constraint, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool>& adpt_block_exist, int& no_of_adpblk, bool use_coff){

	t = new TimingAnalysis;
	t->Init_Timing_Analysis(mu_sigma_local, arrival_time_slack_local, tmp_nodes, place_bench, time_constraint, this->Gates_Map, gate_sizes, gates, no_of_gates, adj_edges, no_of_edges, circuitGraph, gate_parameters, adpt_block_exist, no_of_adpblk, use_coff);
	//t.update_gates_delay(gates, no_of_gates, adj_edges);
	//t.SSTA(gates, no_of_gates, adj_edges, t.mat_class);

}

void ocl_opt::SSTA(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool> adpt_block_exist){
	t->SSTA(mu_sigma_local, arrival_time_slack_local, gate_parameters, adpt_block_exist);
	//t->cu_SSTA();
}

void ocl_opt::pass_pca_plac_info(map<string, vector<double> >& pca_param_local, map<string, vector<float> >& plac_info_local){
	for (int i = 0; i < no_of_gates; ++i){
		for (int j = 0; j < t->no_of_pc; ++j)
			pca_param_local[gates[i].name].push_back(gates[i].k_param[j]*3/0.15/gates[i].gate_resistance/gates[i].gate_capacitance);
		for (int j = 0; j < t->no_of_pc; ++j)
			pca_param_local[gates[i].name].push_back(gates[i].k_param[t->no_of_pc + j] * 3 / 0.05 / gates[i].gate_resistance / gates[i].gate_capacitance);
		plac_info_local[gates[i].name].push_back(gates[i].x / GRID_SIZE);
		plac_info_local[gates[i].name].push_back(gates[i].y / GRID_SIZE);
		plac_info_local[gates[i].name].push_back(gates[i].adapt_blocks_id);
	}

}


bool ocl_opt::load_celllib(char *filename){

	FLT32           bb_tmp;
	FLT32           delay_eql;
	string          old_name;//old footprint
	UINT8D           cell_pin;
	UINT8D           apt_index;
	UINT8D           sol_Cnt = 0;
	INT8D            tune_shift;
	UINT8D           size1;
	LibParserLUT    *pTimingarc;

	LibParser lp(filename);
	double maxTransition = 0.0;
	bool valid = lp.read_default_max_transition(maxTransition);
	H102_ASSERT(valid);
	cout << "The default max transition defined is " << maxTransition << endl;
	cell_lib.resize(MAX_CP_TYPE);

	int readCnt = 0;
	do
	{
		LibParserCellInfo cell;
		valid = lp.read_cell_info(cell);
		if (valid)
		{
			++readCnt;
			//cout << cell << endl ;//ostream& operator<< (ostream& os, LibParserCellInfo& cell) {

			map<string, COMPONENTTYPE>::iterator   gate_lut_iter = gGate_Index.find(cell.footprint);
			if (gate_lut_iter == gGate_Index.end())
			{
				//cout<<"Unknown type: "<<cell.footprint << " skip" << endl;
				continue;//ignore unknown type in contest.lib, e.g. ms00, vss, vcc
			}

			//cell already in local data lib
			if (old_name.compare(cell.footprint))
			{
				old_name = cell.footprint; sol_Cnt = 0;
			}
			else if ((++sol_Cnt % 10) >= 7)//(++sol_Cnt >= MAX_SOL_VALID_NUM)
			{
				continue;
			}

			// simply choose the first pin assuming that all pin have the same input cap, pins[0] is o
			cell_pin = (FFMS == gate_lut_iter->second) ? 2 : 1;
			{
				GATEOPT *pGate_opt = gateopt_cast(new GATEOPT);
				pGate_opt->area = cell.area;
				//config the Vt, gate size(width and length) to the option
				pGate_opt->set_solution(sol_Cnt);

				/************************************************************************/
				/* set input cap, output resist and delay with no adaptivity               */
				/************************************************************************/
				pGate_opt->cap = cell.pins[cell_pin].capacitance;
				H102_ASSERT("o" != cell.pins[cell_pin].name);
				// simply choose the first timing arc thinking that the timing arc are the same and 
				// simply choose the fallDelay LUT without other computations. e.g cell+transition
				pTimingarc = &(cell.timingArcs[0].fallDelay);
				// simply choose the third column's first entry as the offset of the arc delay
				size1 = pTimingarc->tableVals.size() - 1;
				pGate_opt->offset = 0;//pTimingarc->tableVals[0][INPUT_SLEW_INDEX];//offset[0]
				delay_eql = pTimingarc->tableVals[size1][INPUT_SLEW_INDEX] - pGate_opt->offset;

				H102_ASSERT(0 == pTimingarc->loadIndices[0]);
				// simply treat the delay and slew LUT is linear so we get equivalent resistance 
				pGate_opt->resist.push_back(delay_eql / pTimingarc->loadIndices[size1]);//resist[0]
				pGate_opt->off_power.push_back(cell.leakagePower);//off_power[0]
				H102_ASSERT((FFMS == gate_lut_iter->second) || cell.leakagePower);

				/************************************************************************/
				/*  generate the cap, resist and offset based on existing configuration  */
				/************************************************************************/
				for (apt_index = FBB_GRADE_1; apt_index < MAX_BB_GRADE; ++apt_index)
				{
					GET_SHIFT_UNIT(apt_index, tune_shift);
					bb_tmp = pGate_opt->get_leakpower_bb(tune_shift);
					pGate_opt->off_power.push_back(bb_tmp);
					bb_tmp = pGate_opt->get_resist_bb(tune_shift);
					pGate_opt->resist.push_back(bb_tmp);
				}

				cell_lib[gate_lut_iter->second].push_back(pGate_opt);
			}
		}
	} while (valid);

	cout << "Read " << readCnt << " number of library cells" << endl;

	return true;
}

void ocl_opt::load_at_raq(){

	GATEOPT         *pGate_opt;
	for (int i = 0; i < no_of_gates; i++){
		if (gates[i].is_input)
		{
			//assert(!gates[i].is_output);
			if (gates[i].is_output)
				continue;
			pGate_opt = gateopt_cast(new GATEOPT);
			pGate_opt->set_solution(0);//for input driver the solution is the first
			pGate_opt->cap = 0;
			pGate_opt->resist.resize(ADAPTIVITY_1 + 1);//how many resistance are there???????????
			pGate_opt->resist[ADAPTIVITY_0] = 0;//set the input driver.push_back(
			pGate_opt->resist[ADAPTIVITY_1] = 0;
			pGate_opt->offset = 0;
			pGate_opt->area = 0;
			pGate_opt->off_power.resize(ADAPTIVITY_1 + 1);
			pGate_opt->off_power[ADAPTIVITY_0] = 0;//set the off_power of input driver
			pGate_opt->off_power[ADAPTIVITY_1] = 0;
		}
		else if (gates[i].is_output)
		{
			pGate_opt = gateopt_cast(new GATEOPT);

			pGate_opt->set_solution(0);//for output load, the solution is the first
			pGate_opt->cap = 3;//set the output load
			pGate_opt->resist.push_back(0);
			pGate_opt->offset = 0;
			pGate_opt->area = 0;
		}
		else{

		}

	}
}

void ocl_opt::update_gates_info(){
	GATEOPT     *pGate_opt;
	float area_sum = 0;
	vector<float > gate_para;
	gate_para.resize(5); //w, resist, cap, vth, off_power
	for (int i = 0; i < no_of_gates; i++)
	{
		/*
		if (gates[i].is_input || gates[i].is_output)
		{
			continue;
		}
		else
		{
			//for initial timing, use default size
			pGate_opt = gateopt_cast(cell_lib[gates[i].type][0]);
		}
		*/
		pGate_opt = gateopt_cast(cell_lib[gates[i].type][0]);
		gate_para[0] = pGate_opt->w;
		//always give the resistance of the Low Grade ?
		gate_para[1] = pGate_opt->resist[ZERO_GRADE];
		gate_para[2] = pGate_opt->cap;
		gate_para[3] = pGate_opt->vth;
		gate_para[4] = pGate_opt->off_power[ZERO_GRADE];


		gate_sizes[gates[i].id] = gate_para;
		area_sum += pGate_opt->area;
	}
	//add gate sizes info.....! now just for implementing a-SSTA
	cout << "Total Area: " << area_sum << endl;
}

void ocl_opt::calculate_prob(map<int, vector<float> >&   block_prob_local){ t->calculate_prob(block_prob_local);}