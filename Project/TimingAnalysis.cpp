//.....
#define _CRT_SECURE_NO_WARNINGS
#include "TimingAnalysis.h"
#include <oclUtils.h>
#include <shrQATest.h>
#include <math.h>
#include "pdf.h"


#define PROFILING
#define MAX_TILE 1024
#define CL 0
#define N_O_I 5000
#define PARALLEL 1
bool coff = true;

extern "C" Gates_cu* cuSSTA(Gates_cu* gates, int no_of_gates, int* sort, int* edges, int no_of_edges, int no_of_pc, float *eigen_values, \
	float **eigen_vectors, int x_grid, float* k_param, int* l_count, int* l_start, int max_level, Gates* gates_t);

extern "C" void cuMC(Gates_cu2* gates, int no_of_gates, int* edges, int no_of_edges, int no_of_iterations, float** rad, float** rad_w, \
	int grids, int x_grid, int* l_count, int* l_start, int max_level);

bool KValueSort(Gates* i_begin, Gates* i_end){
	return i_begin->k_v < i_end->k_v;
}

bool LabelSort(Gates* i_begin, Gates* i_end){
	return i_begin->label > i_end->label;
}

bool myfunction(int i, int j) { return (i>j); }

//#define DEBUG
void parse_placement_file(string placement_file_name, float temp_max_min[], map<string, Gates* > Gates_Map){

	ifstream placement_file;
	char* word1;
	char* word12;
	char* c_string;

	string a_line_from_file;

	word1 = NULL;
	c_string = NULL;
	//temp_max_min = NULL;

	//temp_max_min    = (float *)malloc(4*sizeof(float));
	temp_max_min[0] = 0;
	temp_max_min[1] = 0;
	temp_max_min[2] = 0;
	temp_max_min[3] = 0;

	c_string = (char*)malloc(1000 * sizeof(char));

	if (!c_string){
		cout << "Memory cannot be assigned " << endl;
		exit(1);
	}
	//word1    = (char*)malloc(sizeof(100));
	//word12 = word1;

	placement_file.open(placement_file_name.c_str());

	if (!placement_file.is_open()){
		printf("The placement file can't be opened!! \n");
		return;
	}
	else{
		printf("The placement file is opened!! \n");
	}

	int x_grid, y_grid;
	while (getline(placement_file, a_line_from_file)){
		strcpy(c_string, a_line_from_file.c_str());

		//cout << c_string << endl;

		word1 = strtok(c_string, " :\t");

#ifdef DEBUG
		//printf("The word read from the placement file is %s \n", word1);
#endif
		if (word1 == NULL){
			continue;
		}

		//cout << word1 << "   " << ckt->Gates_Map.size() << endl;
		/*
		if (ckt->Gates_Map.find(string(word1)) == ckt->Gates_Map.end()){
			continue;
		}
		gate_name = (ckt->Gates_Map.find(word1))->second;
		gate_name->x_grid = atof(strtok(NULL, " :\t"));
		gate_name->y_grid = atof(strtok(NULL, " :\t"));
		*/
		x_grid = atof(strtok(NULL, " :\t"));
		y_grid = atof(strtok(NULL, " :\t"));
		if (Gates_Map.find(word1) == Gates_Map.end())
		{
			//cout << "Gate " << word1 << "is not found!" << endl;
			continue;
		}
		Gates_Map[word1]->x = x_grid;
		Gates_Map[word1]->y = y_grid;

		if (temp_max_min[0] < x_grid){
			temp_max_min[0] = x_grid;
		}

		if (temp_max_min[1] < y_grid){
			temp_max_min[1] = y_grid;
		}

		if (temp_max_min[2] > x_grid){
			temp_max_min[2] = x_grid;
		}

		if (temp_max_min[3] > y_grid){
			temp_max_min[3] = y_grid;
		}


	}

	free(c_string);
	//free(word1);
	placement_file.close();
}

void covmat_class::get_eigen_values(){
	
	this->eigen_values = new float[this->x_grids*this->y_grids];
	this->eigen_vectors = new float*[this->x_grids*this->y_grids];


	for (int i = 0; i < this->x_grids*this->y_grids; i++)
	{	
		this->eigen_vectors[i] = new float[this->x_grids*this->y_grids];
	}

	EigenSolver<MatrixXd> eg;
	eg.compute(this->cov_matrix, true);
	for (int i = 0; i<this->x_grids*this->y_grids; i++){
		this->eigen_values[i] = eg.eigenvalues()(i).real();
		for (int j = 0; j<x_grids*y_grids; j++){
			this->eigen_vectors[i][j] = eg.eigenvectors()(i, j).real();
		}
	}

//#ifdef DEBUG
	cout << "The eigenvalues are: " << eg.eigenvalues().transpose() << endl;
	//getchar();
//#endif

}

void TimingAnalysis::Init_Timing_Analysis(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, igraph_vector_t tmp_nodes, string place_bench, float time_constraint, map<string, Gates* > Gates_Map, map<int, vector<float> > gate_sizes_input, Gates* gates, \
	int no_of_gates, int* adj_edges, int no_of_edges, igraph_t circuitGraph, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool>& adpt_block_exist, int& no_of_adpblk, bool use_coff){

	//set param
	no_of_iterations = N_O_I;
	count = 0;
	gates_t = gates; //may change gates to gates_t within this function
	edges_t = adj_edges;
	this->no_of_gates = no_of_gates;
	this->no_of_edges = no_of_edges;
	this->req_arr_time = time_constraint;
	this->tmp_nodes = tmp_nodes;
	coff = use_coff;
	this->circuitGraph = circuitGraph;
	/*
	int l_count = 1;
	for (int i = 1; i < no_of_gates; ++i)
	{
		l_count++;
		if (gates_t[i].level != gates_t[i - 1].level)
		{
			cout << "Level: " << gates_t[i-1].level << " has " << l_count << " gates " << endl;
			l_count = 1;
		}	
	}
	getchar();
	*/
	naive_schedule();
	cout << "naive schedule finished!" << endl;
	//getchar();

	/*----parameter------*/
	gate_sizes = gate_sizes_input;
	int no_of_level = 4;
	mat_class = new covmat_class;
	const float correlation_level[] = { 0.3, 0.2, 0.15, 0.01 };
	parse_placement_file(place_bench, size_of_die, Gates_Map); // get the xmax, ymax, xmin, ymin

	// i need to remove this stupid iteration =-=!
	for (map<string, Gates*>::iterator it = Gates_Map.begin(); it != Gates_Map.end(); ++it)
	{	
		it->second->x -= size_of_die[2];
		it->second->y -= size_of_die[3];
		assert(it->second->x >= 0);
		assert(it->second->y >= 0);
	}

	mat_class->x_grids = ceil( (size_of_die[0] - size_of_die[2]) /GRID_SIZE);
	mat_class->y_grids = ceil((size_of_die[1] - size_of_die[3]) / GRID_SIZE);
#if DEBUG
	for (int i = 0; i<4; i++)
		cout << size_of_die[i] << endl;
	cout << "X grids :" << mat_class->x_grids << endl;
	cout << "Y grids :" << mat_class->y_grids << endl;
	getchar();
#endif
	int x_grids = mat_class->x_grids;
	int y_grids = mat_class->y_grids;
	mat_class->no_of_pc = x_grids*y_grids;
	no_of_pc = mat_class->no_of_pc;
	MatrixXd tmp(x_grids*y_grids, x_grids*y_grids);
	mat_class->cov_matrix = tmp;
	float** itrtr;
	itrtr = (float **)malloc(x_grids*y_grids*sizeof(float*));
	for (int i = 0; i<x_grids*y_grids; i++){
		*(itrtr + i) = (float *)malloc(x_grids*y_grids*sizeof(float));
	}
	mat_class->MC_matrix = itrtr;
	
/*-----------------Start Generate the covariance matrix--------------------*/
	/*---now we have sized matrix & correlation level, init matrix---*/
	for (int n1 = 0; n1 < x_grids*y_grids; n1++){
		for (int n2 = 0; n2 < x_grids*y_grids; n2++){
			*(*(itrtr + n1) + n2) = 0.0;
			mat_class->cov_matrix(n1, n2) = 0.0;
		}
	}

	/*now begin iteration*/
	int rows = -1; 
	int columns = -1;
	int rows_1 = -1; 
	int columns_1 = -1; 
	int row_diff;
	int column_diff;

	for (int n1 = 0; n1 < x_grids*y_grids; n1++){
		columns_1++;
		if ((n1) % x_grids == 0){
			rows_1++;
			columns_1 = 0;

		}
		rows = -1;
		columns = -1;

		for (int n2 = 0; n2 < x_grids*y_grids; n2++){
			columns++;
			if ((n2) % x_grids == 0){
				rows++;
				columns = 0;
			}

			row_diff = rows - rows_1;
			column_diff = columns - columns_1;

			if (row_diff == 0 && column_diff == 0){
				mat_class->cov_matrix(n1, n2) = 1;
				*(*(itrtr + n1) + n2) = 1;
			}
			else if (abs(row_diff) <= no_of_level && abs(column_diff) <= no_of_level){
				int a;
				a = abs(row_diff) > abs(column_diff) ? abs(row_diff) : abs(column_diff);
				mat_class->cov_matrix(n1, n2) = *(correlation_level + a - 1);
				*(*(itrtr + n1) + n2) = *(correlation_level + a - 1);
			}
			else {
				mat_class->cov_matrix(n1, n2) = 0.0;
				*(*(itrtr + n1) + n2) = 0;
			}
		} // End of n2 For Loop 
	}


	/*---now get the eigens---*/
	mat_class->get_eigen_values();
	cout << "Eigen Done!" << endl;

	int x_adapt_blocks = 4;
	int y_adapt_blocks = 4;
	if (x_adapt_blocks > mat_class->x_grids){
		x_adapt_blocks = mat_class->x_grids;
	}

	if (y_adapt_blocks > mat_class->y_grids){
		y_adapt_blocks = mat_class->y_grids;
	}

	
	adpt_blk = (bool*)malloc(sizeof(bool)*x_adapt_blocks*y_adapt_blocks);
	no_of_adpt_blk = x_adapt_blocks*y_adapt_blocks;
	no_of_adpblk = no_of_adpt_blk;

	for (int i = 0; i < no_of_adpt_blk; i++){
		Block_inst* block;
		block = new Block_inst;
		block->id = i;
		block->adaptive = false;
		Blocks.push_back(block);
		adpt_blk[i] = false;
	}

	int x, y;

	x = ceil((float)mat_class->x_grids / x_adapt_blocks);
	y = ceil((float)mat_class->y_grids / y_adapt_blocks);
	for (int i = 0; i < no_of_gates; i++){
		gates[i].adapt_blocks_id = (int)(gates[i].x / GRID_SIZE) / x + (int)(gates[i].y / GRID_SIZE) / y * (x_adapt_blocks);
	}
	
	

	find_blk_PO();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	//this->prede_test();

#ifdef	PROFILING
	shrDeltaT(0);
#endif
	this->SSTA(mu_sigma_local, arrival_time_slack_local, gate_parameters, adpt_block_exist);

#ifdef	PROFILING
	double kernel_time = shrDeltaT(0);
#endif

	
	
	
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
#ifdef	PROFILING
	shrDeltaT(0);
#endif
	//Monte_Carlo_Gen();
	//Monte_Carlo(mat_class, true, 1);
#ifdef	PROFILING
	kernel_time = shrDeltaT(0);
	//cout << "CPU MC Running time(s):" << kernel_time << endl;
#endif
	//getchar();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	backup();
#if PARALLEL
	//backup();

	//cu_MC();
	cu_SSTA();
	
#endif
	get_max_delay_prob();


}


void TimingAnalysis::update_gates_delay(Gates* gates, int no_of_gates, int* edges, std::map<std::string, std::vector<float> > gate_parameters){

	//import R,C value here, elmore delay just use the two para, for now
	for (int i = 0; i < no_of_gates; i++){
		// I ignore the PI PO here, caz in my model the PI PO are real first/last level gates!
		/*
		if (gates[i].is_input)
		{
			gates[i].gate_resistance = 1;
			gates[i].gate_capacitance = 0;
		}
		else if (gates[i].is_output)
		{
			gates[i].gate_resistance = 0;
			gates[i].gate_capacitance = 1;
		}
		*/
		//else{
		gates[i].gate_resistance = gate_parameters[gates[i].name][1];
		gates[i].gate_capacitance = gate_parameters[gates[i].name][2];
		//}
	}

	/*-------- just use Elmore------------*/
	for (int i = 0; i < no_of_gates; i++)
	{
		gates[i].load_capacitance = 0.0;
		gates[i].delay = 0.0;
		gates[i].delay_mu = 0;
		gates[i].delay_sigma = 0;
		if (gates[i].is_output)
			gates[i].load_capacitance = 3;
		else {
			for (int j = gates[i].start_out; j < gates[i].start_out + gates[i].no_of_out; j++){
				//gates' output gate cap sum
				gates[i].load_capacitance += gates[edges[j]].gate_capacitance;
			}
		}
			//comput delay using elmore model
			gates[i].delay = gates[i].gate_resistance * gates[i].load_capacitance;
	}
}

void TimingAnalysis::update_gates_delay2(Gates2* gates, int no_of_gates, int* edges){

	//import R,C value here, elmore delay just use the two para, for now
	for (int i = 0; i < no_of_gates; i++){
		// I ignore the PI PO here, caz in my model the PI PO are real first/last level gates!
		/*
		if (gates[i].is_input)
		{
		gates[i].gate_resistance = 1;
		gates[i].gate_capacitance = 0;
		}
		else if (gates[i].is_output)
		{
		gates[i].gate_resistance = 0;
		gates[i].gate_capacitance = 1;
		}
		*/
		//else{
		gates[i].gate_resistance = gate_sizes[gates[i].id][1];
		gates[i].gate_capacitance = gate_sizes[gates[i].id][2];
		//}
	}

	/*-------- just use Elmore------------*/
	for (int i = 0; i < no_of_gates; i++)
	{
		gates[i].load_capacitance = 0.0;
		gates[i].delay = 0.0;
		gates[i].delay_mu = 0;
		gates[i].delay_sigma = 0;
		if (gates[i].is_output)
			gates[i].load_capacitance = 3;
		else {
			for (int j = gates[i].start_out; j < gates[i].start_out + gates[i].no_of_out; j++){
				//gates' output gate cap sum
				gates[i].load_capacitance += gates[edges[j]].gate_capacitance;
			}
		}
		//comput delay using elmore model
		gates[i].delay = gates[i].gate_resistance * gates[i].load_capacitance;
	}
}

//#ifdef aa
void TimingAnalysis::get_req_slack(Gates* gates, int no_of_gates, int* edges, double time_constraint){
	/*--------, get the require time----------*/
	for (int i = no_of_gates - 1; i >= 0; i--){
		if (gates[i].is_output){
			gates[i].req_arrive_time = time_constraint - gates[i].gate_mu- 3*gates[i].gate_sigma;
			continue;
		}
		for (int j = gates[i].start_in; j < gates[i].start_in + gates[i].no_of_in; j++)
		{	
			gates[edges[j]].req_arrive_time = gates[i].req_arrive_time - gates[i].gate_mu - 3 * gates[i].gate_sigma;
			// the frontier's req_arrive_time should be positive!!!
			assert(gates[edges[j]].req_arrive_time >= 0);
		}
	}
}


void TimingAnalysis::SSTA(map<string, pair<double, double> >& mu_sigma_local, map<string, pair<double, double> >& arrival_time_slack_local, std::map<std::string, std::vector<float> > gate_parameters, map<int, bool> adpt_block_exist){
	/*---this is hell, let me do it step by step , currently I don't consider the adaptive - 。- ---*/
	/*------------------------Cov_matrix, gates_parameters OK, first, let's get the k values(L,W) for each gates-----------------------*/

#ifdef	PROFILING
	shrDeltaT(0);
#endif

	this->update_gates_delay(gates_t, no_of_gates, edges_t, gate_parameters);

	for (int i = 0; i < no_of_adpt_blk; ++i)
		adpt_blk[i] = adpt_block_exist[i];

	k_para_matrix = (float*)malloc(sizeof(float) * no_of_gates * no_of_pc * 2);
	int offset = 2 * no_of_pc;
	for (int ii = 0; ii < no_of_gates; ii++)
	{
		int i = tmp_nodes2[ii];
		float adpt = 0.0f;
		if (adpt_blk[gates_t[i].adapt_blocks_id])
			adpt = -0.05f;

		double sigma_of_delay = 0.0; // this is only the sigma of the gate's delay, just for debugging purpose
		gates_t[i].gate_mu = 0.0;  // This is to init mu
		gates_t[i].gate_sigma = 0.0; // This is to init sigma
		gates_t[i].gate_mu = gates_t[i].delay / pow(1.2 - (0.2 + adpt), 2);  // nominal delay
		gates_t[i].k_param = (float*)malloc(sizeof(float) * 2 * mat_class->no_of_pc);

		/*------------ start k-------------------*/

		int row = floor(gates_t[i].y / GRID_SIZE);
		int column = floor(gates_t[i].x / GRID_SIZE);
		assert(row >= 0 && row < mat_class->no_of_pc);
		assert(column >= 0 && column < mat_class->no_of_pc);

		int i_of_j = row*(mat_class->x_grids) + column;


		//These two need to get from .... i don't know now
		double sigma_of_L = 1.0;
		double sigma_of_W = 1.0;

		double k_tmp = 0.0;

		//gates[i].k_para = new double[2 * cov->no_of_pc];
		//gates[i].no_of_k = 2 * cov->no_of_pc;
		//double* p = gates[i].k_para;


		for (int j = 0; j < mat_class->no_of_pc; j++)
		{
			// I'll figure out the minus situation
			if (mat_class->eigen_values[j] < 0){
				mat_class->eigen_values[j] = 0;
			}

			//L = sqr eg_v * eg_vec * sigma, the dRdL is the constain for a specific size of gate, i'll update it
			k_tmp = (0.15 / 3)* gates_t[i].delay * sqrt(mat_class->eigen_values[j]) * mat_class->eigen_vectors[i_of_j][j] / pow(1.2 - (0.2 + adpt), 2);// *sigma_of_L / pow(1.2 - (0.2 + adpt), 2);
			//*(p++) = k_tmp;
			gates_t[i].k_param[j] = k_tmp;
			k_para_matrix[i*offset + j] = k_tmp;
			sigma_of_delay += k_tmp * k_tmp;
		//	cout << k_tmp << endl;
		}

		for (int j = 0; j < no_of_pc; j++)
		{
			// I'll figure out the minus situation
			if (mat_class->eigen_values[j] < 0){
				mat_class->eigen_values[j] = 0;
			}
			//W = sqr eg_v * eg_vec * sigma, the dRdW is the constant for a specific size of gate, i'll update it
			k_tmp = -(0.08 / 3)* gates_t[i].delay * sqrt(mat_class->eigen_values[j]) * mat_class->eigen_vectors[i_of_j][j] / pow(1.2 - (0.2 + adpt), 2);// *sigma_of_W / pow(1.2 - (0.2 + adpt), 2);
			//*(p++) = k_tmp;
			gates_t[i].k_param[j + no_of_pc] = k_tmp;
			k_para_matrix[i*offset + no_of_pc + j] = k_tmp;
			sigma_of_delay += k_tmp * k_tmp;
		//	cout << k_tmp << endl;
		}
	
		/*--------------------Done with K values------------------------*/

		/*----get gate[i]'s sigma----*/
		assert(sigma_of_delay >= 0);
		gates_t[i].gate_sigma = (float)sqrt(sigma_of_delay);
		if (gates_t[i].gate_sigma > 10000000 || gates_t[i].gate_sigma < -10000000){
			cout << "Wrong sigma in gate " << i << endl;
			getchar();
		}


	}



	/*----------- Do Max and Sum function --------------*/
	/*------------ I decide not to use recursive function, because the gates have been sorted by level -------------------*/
	int current_level = 0;
	int i;
	for (int ii = 0; ii < no_of_gates; ii++)
	{
		i = tmp_nodes2[ii];
		if (gates_t[i].is_input){
			//PI use the gate's delay as accumulated delay
			gates_t[i].delay_mu = gates_t[i].gate_mu;
			gates_t[i].delay_sigma = gates_t[i].gate_sigma;
		}
		else{
				// need to add a flag here, to check if all the inputs have been visited
				//--max--//
			mu_sigma_struct1 max_strc;
			max_strc = max_function(&gates_t[i]);  //no change supposed to be on gates[i]
				//--sum--//
			sum_function(max_strc, &gates_t[i]); // gates[i] param should changes
		}

		mu_sigma_local[gates_t[i].name].first = gates_t[i].delay_mu;
		mu_sigma_local[gates_t[i].name].second = gates_t[i].delay_sigma;
		//cout << gates[i].delay_mu << " " << gates[i].delay_sigma << endl;
		//getchar();
	}
	
	this->update_req_arr_time_and_slack(arrival_time_slack_local);
}
//#endif

void TimingAnalysis::cu_SSTA(){
	//this->update_gates_delay(gates_cu, no_of_gates, edges_t);

	gates_cu_p = cuSSTA(gates_cu, no_of_gates, tmp_nodes2, edges_cu, no_of_edges, mat_class->no_of_pc, mat_class->eigen_values, mat_class->eigen_vectors,\
		mat_class->x_grids, k_para_matrix, l_count, l_start, max_level, gates_t);

}

void TimingAnalysis::cu_MC(){
	//this->update_gates_delay(gates_cu, no_of_gates, edges_t);

	int k;
	k = mat_class->x_grids*mat_class->y_grids;
	cuMC(gates_cu2, no_of_gates, edges_cu, no_of_edges, no_of_iterations, mat_class->random_numbers, mat_class->random_numbers_width, \
		k, mat_class->x_grids, l_count, l_start, max_level);

}


void TimingAnalysis::cl_SSTA(Gates2* gates, int no_of_gates, int* edges, int no_of_edges, covmat_class* cov){
	/*---this is hell, let me do it step by step , currently I don't consider the adaptive - 。- ---*/
	/*------------------------Cov_matrix, gates_parameters OK, first, let's get the k values(L,W) for each gates-----------------------*/
	this->update_gates_delay2(gates, no_of_gates, edges);


	/*-------------Preparing cl data-----------*/
	int cl_no_of_pc;
	cl_no_of_pc = cov->no_of_pc;
	float* cl_eigen_vec;
	cl_eigen_vec = (float*)malloc(sizeof(float)*(cl_no_of_pc*cl_no_of_pc));

	int idx = 0;
	for (int i = 0; i < cl_no_of_pc; i++)
		for (int j = 0; j < cl_no_of_pc; j++){
		cl_eigen_vec[idx++] = cov->eigen_vectors[i][j];
	}

	/*--------------OpenCL starts here!---------------*/
		int G = no_of_gates;  //global size
		int L; //local size
		int no_of_work_groups;

		try{

			vector<cl::Device> devices;
			//unsigned numDevices = getDeviceList(devices);
			cl::Context context(CL_DEVICE_TYPE_GPU);
			cl::CommandQueue queue(context);


			cl::Buffer d_gates(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(Gates2)*no_of_gates, gates);
			cl::Buffer d_edges(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(int)*no_of_edges, edges);
			cl::Buffer d_eigen_v(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*cl_no_of_pc, cov->eigen_values);
			cl::Buffer d_eigen_vec(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(float)*cl_no_of_pc*cl_no_of_pc, cl_eigen_vec);
			cl::Buffer d_flag(context, CL_MEM_WRITE_ONLY, sizeof(float));
			cl::Buffer d_matrix(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(float)*no_of_gates*no_of_pc * 2, k_para_matrix);
			cl::Buffer d_max_k(context, CL_MEM_READ_WRITE, sizeof(float)*no_of_pc * 2);

			queue.enqueueWriteBuffer(d_gates, CL_TRUE, 0, sizeof(Gates2)*no_of_gates, gates, 0, 0);
			queue.enqueueWriteBuffer(d_edges, CL_TRUE, 0, sizeof(int)*no_of_edges, edges, 0, 0);
			queue.enqueueWriteBuffer(d_eigen_v, CL_TRUE, 0, sizeof(float)*cl_no_of_pc, cov->eigen_values, 0, 0);
			queue.enqueueWriteBuffer(d_eigen_vec, CL_TRUE, 0, sizeof(float)*cl_no_of_pc*cl_no_of_pc, cl_eigen_vec, 0, 0);
			queue.enqueueWriteBuffer(d_matrix, CL_TRUE, 0, sizeof(float)*no_of_gates*no_of_pc * 2, k_para_matrix, 0, 0);

			cl::Program program(context, util::loadProgram("Kernels.cl"), true);
//			cl::Kernel ko_ssta(program, "ssta");

			// Get the work group size, so-called local size
			//cl::Device device = devices[0];

			L = 1024;
			if (G < L*ssta_iter){
				L = G / ssta_iter;
				no_of_work_groups = 1;
			}
			else{
				if (!(G % (L*ssta_iter)))
					no_of_work_groups = G / (L * ssta_iter);
				else
				{
					while ((G % (L*ssta_iter)))
						G--;
					no_of_work_groups = G / (L * ssta_iter);
				}

				if (no_of_work_groups < 1) {
					no_of_work_groups = 2;
					assert(!(G / (no_of_work_groups*ssta_iter)));
					L = G / (no_of_work_groups*ssta_iter);
				}
			}
			cl::make_kernel<cl::Buffer/*gates*/, cl::Buffer/*eigen_v*/, cl::Buffer/*eigen_vec*/, int, int, int, cl::Buffer> ssta1(program, "ssta1");

			printf("\n========== Now starting OpenCL program, global %d, local %d ===========\n", G, L);

			

			//cl::NDRange global(G/ssta_iter);
			//cl::NDRange local(L);
			cl::NDRange global(no_of_gates);

#ifdef	PROFILING
			shrDeltaT(0);
#endif
			int xx = cov->x_grids;
			ssta1(cl::EnqueueArgs(queue, global), d_gates, d_eigen_v, d_eigen_vec, cl_no_of_pc, ssta_iter, xx, d_flag);
			queue.finish();
#ifdef	PROFILING
			double total;
			double kernel_time = shrDeltaT(0);
			cout << "1st GPU Running time(s):" << kernel_time << endl;
			total = kernel_time;
#endif
			//float flag = 0.0f;
			//queue.enqueueReadBuffer(d_flag, CL_TRUE, 0, sizeof(float), &flag, 0, 0);

			//while ((G % (ssta_iter2)))
				//G--;
			//cl::NDRange global1(G);
			//cl::Buffer d_gates1(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(Gates)*no_of_gates, gates);
			//queue.enqueueWriteBuffer(d_gates1, CL_TRUE, 0, sizeof(Gates)*no_of_gates, gates, 0, 0);
			cl::make_kernel<cl::Buffer/*gates*/, cl::Buffer,/*edges*/ int, int, int, cl::Buffer, cl::Buffer, cl::Buffer> ssta2(program, "ssta2");

			/*
			int flag_level = 0;
			gates[0].nn = 0;
			for (int i = 1; i < no_of_gates; ++i){
				if (gates[i - 1].level == gates[i].level)
					if (abs(gates[i - 1].no_of_in - gates[i].no_of_in) == 0){
						gates[i].nn = flag_level;
					}
					else
					{
						gates[i].nn = ++flag_level;
					}
				else{
					gates[i].nn = ++flag_level;
				}
			}
			for (int i = 0; i < no_of_gates; ++i){
				gates[i].level = gates[i].nn;
			}
			*/




#ifdef	PROFILING
			shrDeltaT(0);
#endif
			//cl::NDRange local2(L);
			//cl::NDRange global2(550);
			//cout << "Gates2: " << sizeof(Gates2) << endl;
			//cl::LocalSpaceArg localmem = cl::Local(sizeof(Gates2) * 550);

			//int current_level = 0;
			//while (current_level <= gates[no_of_gates-1].level){
			ssta2(cl::EnqueueArgs(queue, global), d_gates, d_edges, gates[no_of_gates - 1].level, ssta_iter, no_of_pc, d_matrix, d_max_k, d_flag);
			//}
			/*
			int begin = 0;
			int end = 39;
			while (end < no_of_gates){
				ssta2(cl::EnqueueArgs(queue, global2), d_gates, d_edges, begin, end, ssta_iter, no_of_pc, d_matrix, d_max_k, d_flag, localmem);
				begin = end + 1;
				end += 40;
			}

			if (begin < no_of_gates)
			{
				end = no_of_gates - 1;
				ssta2(cl::EnqueueArgs(queue, global2), d_gates, d_edges, begin, end, ssta_iter, no_of_pc, d_matrix, d_max_k, d_flag, localmem);
			}
			*/
			//float flag = 3.0f;
			//queue.enqueueReadBuffer(d_flag, CL_TRUE, 0, sizeof(float), &flag, 0, 0);
			queue.finish();
			//cout << "flag: " << flag << endl;
#ifdef	PROFILING
			kernel_time = shrDeltaT(0);
			cout << "2nd GPU Running time(s):" << kernel_time << endl;
			total += kernel_time;
			cout << "Total: " << total << endl;
			getchar();
#endif
			queue.enqueueReadBuffer(d_gates, CL_TRUE, 0, sizeof(Gates2)*no_of_gates, gates, 0, 0);
			
			int count1 = 0;
			for (int i = 0; i < no_of_gates; i++){
				if ((gates[i].delay_mu * 0.2) < gates[i].delay_sigma)
					count1++;
			}

			cout <<"Unnormal sigma count: " << count1 << endl;
			/*
			for (int i = 0; i < no_of_gates; i++){
				cout << "gate " << i << " level: "<< gates[i].level <<" delay: " << gates[i].delay_mu << " sigma: " << gates[i].delay_sigma << endl;
				getchar();
			}
			*/
		}
		catch (cl::Error err){
			std::cout << "Exception\n";
			std::cerr << "ERROR: "
				<< err.what()
				<< "("
				//<< err_code(err.err())
				<< ")"
				<< err.err()
				<< std::endl;
			getchar();
		}

}

//#ifdef aa
void TimingAnalysis::sum_function(mu_sigma_struct1 max_strc, Gates* currentG){
	
	currentG->delay_mu = currentG->gate_mu + max_strc.mu;
	float tmp1 = 0.0f;
	float tmp2 = 0.0f;
	float sigma = 0.0f;
	for (int i = 0; i < 2 * no_of_pc; i++)
	{
		tmp1 += pow(currentG->k_param[i],2);
		tmp2 += pow(max_strc.k_parameters[i],2);
		//sigma += pow(tmp, 2);
	}
	//if (tmp2>200)
	//cout << tmp1 << " " << tmp2 << endl;
	currentG->delay_sigma = 2*sqrt(tmp1) + sqrt(tmp2);
}

mu_sigma_struct1 TimingAnalysis::max_function(Gates* currentG){
	float mu_1, mu_2, sigma_1, sigma_2;
	mu_sigma_struct1 max_strc;
	max_strc.mu = gates_t[edges_t[currentG->start_in]].delay_mu;
	max_strc.sigma = gates_t[edges_t[currentG->start_in]].delay_sigma;
	for (int j = 0; j < 2 * no_of_pc; j++){
		max_strc.k_parameters.push_back(gates_t[edges_t[currentG->start_in]].k_param[j]);
		assert(k_para_matrix[gates_t[edges_t[currentG->start_in]].id * 2 * no_of_pc + j] == gates_t[edges_t[currentG->start_in]].k_param[j]);
	}


	
	for (int i = currentG->start_in+1; i < currentG->start_in + currentG->no_of_in; i++){
		/*
		if (currentG->level != (gates_t[edges_t[i]].level + 1))
		{	
			cout << currentG->level << " " << gates_t[edges_t[i]].level << endl;
			getchar();
		}
		*/
		mu_1 = max_strc.mu;
		sigma_1 = max_strc.sigma;

		mu_2 = gates_t[edges_t[i]].delay_mu;
		sigma_2 = gates_t[edges_t[i]].delay_sigma;

		
		if (mu_1 - 3 * sigma_1 > mu_2 + 3 * sigma_2)
		{
			continue;
		}
	
		if (mu_1 + 3 * sigma_1 < mu_2 - 3 * sigma_2)
		{
			max_strc.mu = mu_2;
			max_strc.sigma = sigma_2;
			for (int j = 0; j < 2 * no_of_pc; j++){
				max_strc.k_parameters.push_back(gates_t[edges_t[i]].k_param[j]);
			}
			continue;
		}
		
		
		//step 2
		float co_variance = 0.0f;
		float correlation = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			co_variance += max_strc.k_parameters[j] * gates_t[edges_t[i]].k_param[j];
		}
		
		//assert((sigma_1 * sigma_2));
		if ((sigma_1 * sigma_2))
			correlation = co_variance / (sigma_1 * sigma_2);
		else
			correlation = co_variance / 0.001;
		if (correlation > 0.99 && abs(sigma_1 - sigma_2) < 0.1){
			if (mu_1 > mu_2){
				continue;
			}
			else{
				max_strc.mu = mu_2;
				max_strc.sigma = sigma_2;
				for (int j = 0; j < 2 * no_of_pc; j++){
					max_strc.k_parameters.push_back(gates_t[edges_t[i]].k_param[j]);
				}
				continue;
			}
		}
		


		//step 3
		float alpha = sqrt(abs(pow(sigma_1, 2) + pow(sigma_2, 2) - 2 * co_variance));
		
		float beta = (mu_1 - mu_2) / alpha;

		float phi = pow(2.718281828f, -beta*beta / 2) / sqrt(2 * 3.141592654f);

		float phi_intg = integrate1(beta);
		//float phi_intg1 = integrate2(beta);
	//	if (abs(phi_intg - phi_intg1) > 0.1){
	//		cout << beta << " , " << phi_intg << " , " << phi_intg1 << endl;
	//		getchar();
	//	}

		float phi_intg_m = integrate1(-beta);
	//	float phi_intg_m1 = integrate2(-beta);

		float sigma_3, mu_3;
	//	float sigma_31, mu_31;

		mu_3 = mu_1 * phi_intg + mu_2 * phi_intg_m +alpha * phi;
	//	mu_31 = mu_1 * phi_intg1 + mu_2 * phi_intg_m1 + alpha * phi;

		float sigma_tmp = (pow(mu_1, 2) + pow(sigma_1, 2)) * phi_intg + (pow(mu_2, 2) + pow(sigma_2, 2)) * phi_intg_m + (mu_1 + mu_2) * alpha * phi - mu_3*mu_3;
	//	float sigma_tmp2 = (pow(mu_1, 2) + pow(sigma_1, 2)) * phi_intg1 + (pow(mu_2, 2) + pow(sigma_2, 2)) * phi_intg_m1 + (mu_1 + mu_2) * alpha * phi - mu_3*mu_3;

		sigma_3 = sqrt(abs(sigma_tmp));
	//	sigma_31 = sqrt(abs(sigma_tmp2));

//		if (abs(sigma_3 - sigma_31) > 10 || abs(mu_3 - mu_31) > 10){
	//		cout << mu_3 << " " << mu_31 << endl;
		//	cout << sigma_3 << " " << sigma_31 << endl;
			//getchar();
	//	}
		
		//step 4
		vector<float> a_r;
		float S0 = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			float r_1 = max_strc.k_parameters[j];
			float r_2 = gates_t[edges_t[i]].k_param[j];
			//assert(sigma_3>0);
			if (!sigma_3)
				sigma_3 = 0.1;
			float ar_tmp = (sigma_1*r_1*phi_intg + sigma_2*r_2*phi_intg_m) / sigma_3;
			a_r.push_back(ar_tmp);
			S0 += ar_tmp * ar_tmp;
		}

		for (int j = 0; j < 2 * no_of_pc; j++){
			a_r[j] = a_r[j] * sigma_3 / sqrt(abs(S0));
		}
		
		
		max_strc.mu = mu_3;
			max_strc.sigma = sigma_3;

		max_strc.k_parameters = a_r;
		if (mu_3*0.2 < sigma_3)
			count++;
			
	}
	
	return max_strc;
}

float TimingAnalysis::integrate(float beta){
	
	float step = 0.01f;
	float down_idx = -5;
	int up_idx = (int)(beta - down_idx) / step;
	float micro = down_idx;

	float sum = 0.0f;
	float tmp;

	for (int ii = 0; ii < up_idx; ii++){
		tmp = pow(micro, 2) / 2;
		sum += pow(2.718281828f, -tmp) * step;
		micro += step;
	}
	
	return sum / (sqrt(2 * 3.141592654f));

}

double TimingAnalysis::integrate2(double x){
	float step = 0.01;
	int iterations;
	double index;
	double lli;
	double mean = 0.0;
	double sigma = 1.0;
	double sum = 0.0;

	lli = -2.0 + mean - 3 * sigma;

	if (x < lli){
		return 0.00;
	}

	iterations = (int)(ceil((x - lli) / step));

	for (int i = 0; i < iterations; i++){
		index = pow(((lli + i*step - step / 2) - mean) / sigma, 2) / 2;
		sum += step*pow(2.718281828, -index);
	}

	return sum / (sigma * sqrt(2 * 3.141592654));
}

float TimingAnalysis::integrate1(float beta){

	bool positive;
	float b = abs(beta);
	float result;
	int flag = 0;

	int i, j;
	for (i = 0; i < 36; ++i){
		for (j = 0; j < 10; ++j){
			if (b >(row[i] + col[j]))
				continue;
			else{
				flag = 1;
				break;
			}
		}
		if (flag)
			break;
	}

	if ((i - 1) < 0 || (j - 1) < 0)
		result = pdf_m[i][j];
	else if (i > 35 || j > 35)
		result = pdf_m[i - 1][j - 1];
	else
		result = 0.5*(pdf_m[i][j] + pdf_m[i - 1][j - 1]);


	if (beta >= 0){
		return result + 0.5;
	}
	else{
		return 0.5 - result;
	}
}
//#endif

void TimingAnalysis::update_req_arr_time_and_slack(map<string, pair<double, double> >& arrival_time_slack_local){

	
	for (int ii = this->no_of_gates - 1; ii >= 0; ii--){

		int i = tmp_nodes2[ii];
		if (gates_t[i].is_output)
		{
			gates_t[i].req_arrive_time = this->req_arr_time;
			for (int j = gates_t[i].start_in; j < gates_t[i].start_in + gates_t[i].no_of_in; j++){
				gates_t[edges_t[j]].req_arrive_time = gates_t[i].req_arrive_time - gates_t[edges_t[j]].gate_mu - 3 * gates_t[edges_t[j]].gate_sigma;
			}
		}
		else{
			for (int j = gates_t[i].start_in; j < gates_t[i].start_in + gates_t[i].no_of_in; j++){
				gates_t[edges_t[j]].req_arrive_time = gates_t[i].req_arrive_time - gates_t[edges_t[j]].gate_mu - 3 * gates_t[edges_t[j]].gate_sigma;
			}
		}
	}
	
	for (int i = 0; i < this->no_of_gates; i++){
		this->req_arr_time;
		gates_t[i].slack = gates_t[i].req_arrive_time - gates_t[i].delay_mu - gates_t[i].delay_sigma * 3;
		arrival_time_slack_local[gates_t[i].name].first = gates_t[i].req_arrive_time;
		arrival_time_slack_local[gates_t[i].name].second = gates_t[i].slack;
	}
}

void TimingAnalysis::find_blk_PO(){

	for (int i = 0; i < no_of_gates; i++){
		for (int j = gates_t[i].start_out; j < gates_t[i].start_out + gates_t[i].no_of_out; j++){
			if (gates_t[i].adapt_blocks_id != gates_t[edges_t[j]].adapt_blocks_id)
			{
				int blockID = gates_t[i].adapt_blocks_id;
				Blocks[blockID]->gates_idx.push_back(i);
			}
		}
		if (gates_t[i].is_output){
			int blockID = gates_t[i].adapt_blocks_id;
			Blocks[blockID]->gates_idx.push_back(i);
		}

	}

}


// may change...
void TimingAnalysis::calculate_prob(map<int, vector<float> > &probability){

	//this->update_req_arr_time_and_slack();

	vector<Block_inst* >::iterator blk_iter;
	for (blk_iter = Blocks.begin(); blk_iter != Blocks.end(); ++blk_iter){
		float min_slack = 5000.0f;
		float _mu = 0.0f, _sigma = 0.0f, _rat = 0.0f;
		vector<int>::iterator it = (*blk_iter)->gates_idx.begin();
		for (; it != (*blk_iter)->gates_idx.end(); ++it){
			if (gates_t[(*it)].slack < min_slack){
				min_slack = gates_t[(*it)].slack;
				_mu = gates_t[(*it)].delay_mu;
				_sigma = gates_t[(*it)].delay_sigma;
				_rat = gates_t[(*it)].req_arrive_time;
			}
		}
		(*blk_iter)->w_slack = min_slack;
	//	cout << "Min Slack: "<< min_slack << endl;

		float phi = 0.0f;
		if (min_slack >= 0){
			probability[(*blk_iter)->id].push_back(1); //Low
			probability[(*blk_iter)->id].push_back(0); //High
			probability[(*blk_iter)->id].push_back(0); //useless
		}
		else{
			phi = (_rat - _mu) / _sigma;
			probability[(*blk_iter)->id].push_back(integrate1(phi)); //Low
			probability[(*blk_iter)->id].push_back(1 - integrate1(phi)); //High
			probability[(*blk_iter)->id].push_back(0); //useless
		}
	//	cout << "Prob: "<< probability[(*blk_iter)->id][1] << endl;
	//////////////////////////////////
	}
}

void TimingAnalysis::adaptive_SSTA(){
	//--------------------update slack--------------------//
	//update_req_arr_time_and_slack();
	//calculate probability after the FIRST SSTA!!!!
	vector<Block_inst* >::iterator blk_iter;
	for (blk_iter = Blocks.begin(); blk_iter != Blocks.end(); ++blk_iter){
		cout << (*blk_iter)->id << "Prob: " << (*blk_iter)->prob << endl;
	}
	
}



void TimingAnalysis::calculate_yield1(Gates* gates){

	vector<Gates> PO_gates;
	int offset = 2*no_of_pc;
	for (int i = 0; i < no_of_gates; i++){
		if (gates[i].is_output)
			PO_gates.push_back(gates[i]);
	}

	float mu_1, mu_2, sigma_1, sigma_2;
	mu_sigma_struct1 max_strc;
	max_strc.mu = PO_gates[0].delay_mu;
	max_strc.sigma = PO_gates[0].delay_sigma;
	for (int j = 0; j < 2 * no_of_pc; j++){
		max_strc.k_parameters.push_back(k_para_matrix[index_tmp[PO_gates[0].id]*offset + j]);
		//max_strc.k_parameters.push_back(PO_gates[0].k_param[j]);
	}

	/*
	for (int i = 0; i < PO_gates.size(); ++i){
		cout << PO_gates[i].delay_mu << " " << PO_gates[i].delay_sigma << endl;
		
	}
	*/
	vector<Gates>::iterator itr = PO_gates.begin();

	for (itr++; itr != PO_gates.end(); itr++){

		mu_1 = max_strc.mu;
		sigma_1 = max_strc.sigma;

		mu_2 = (*itr).delay_mu;
		sigma_2 = (*itr).delay_sigma;


		if (mu_1 - 3 * sigma_1 > mu_2 + 3 * sigma_2)
		{
			continue;
		}

		if (mu_1 + 3 * sigma_1 < mu_2 - 3 * sigma_2)
		{
			max_strc.mu = mu_2;
			max_strc.sigma = sigma_2;
			for (int j = 0; j < 2 * no_of_pc; j++){
				max_strc.k_parameters.push_back(k_para_matrix[index_tmp[(*itr).id] * offset + j]);
			}
			continue;
		}


		//step 2
		float co_variance = 0.0f;
		float correlation = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			co_variance += max_strc.k_parameters[j] * k_para_matrix[index_tmp[(*itr).id] * offset + j];
		}
		//cout << sigma_1 << " " << sigma_2 << " " << sigma_1*sigma_2 << endl;
		assert((sigma_1 * sigma_2));
		correlation = co_variance / (sigma_1 * sigma_2);
		if (correlation > 0.99 && abs(sigma_1 - sigma_2) < 0.1){
			if (mu_1 > mu_2){
				continue;
			}
			else{
				max_strc.mu = mu_2;
				max_strc.sigma = sigma_2;
				for (int j = 0; j < 2 * no_of_pc; j++){
					max_strc.k_parameters.push_back(k_para_matrix[index_tmp[(*itr).id] * offset + j]);
				}
				continue;
			}

		}

		//step 3
		float alpha = sqrt(abs(pow(sigma_1, 2) + pow(sigma_2, 2) - 2 * co_variance));
		if (!(alpha>0))
			cout << alpha << endl;
		float beta = (mu_1 - mu_2) / alpha;

		float phi = pow(2.718281828, -beta*beta / 2) / sqrt(2 * 3.141592654);
		float phi_intg = integrate1(beta);
		float phi_intg_m = integrate1(-beta);

		float sigma_3, mu_3;

		mu_3 = mu_1 * phi_intg + mu_2 * phi_intg_m +alpha * phi;
		float sigma_tmp = (pow(mu_1, 2) + pow(sigma_1, 2)) * phi_intg + (pow(mu_2, 2) + pow(sigma_2, 2)) * phi_intg_m + (mu_1 + mu_2) * alpha * phi - mu_3*mu_3;
		sigma_3 = sqrt(abs(sigma_tmp));

		//step 4
		vector<float> a_r;
		float S0 = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			float r_1 = max_strc.k_parameters[j];
			float r_2 = k_para_matrix[tmp_nodes2[(*itr).id]*offset + j];
			assert(sigma_3);
			float ar_tmp = (sigma_1*r_1*phi_intg + sigma_2*r_2*phi_intg_m) / sigma_3;
			a_r.push_back(ar_tmp);
			S0 += ar_tmp * ar_tmp;
		}

		for (int j = 0; j < 2 * no_of_pc; j++){
			a_r[j] = a_r[j] * sigma_3 / sqrt(abs(S0));
		}

		max_strc.mu = mu_3;
		max_strc.sigma = sigma_3;
		max_strc.k_parameters = a_r;
	}

	cout << max_strc.mu << " " << max_strc.sigma << endl;
	float phi = 0.0f;
	float yield = 0.0f;
	phi = (req_arr_time - max_strc.mu) / max_strc.sigma;
	cout << phi << endl;
	yield = integrate1(phi);
	cout << "Yield: " << yield << endl;

}

void TimingAnalysis::calculate_yield1(Gates_cu* gates){

	vector<Gates_cu> PO_gates;
	int offset = 2 * no_of_pc;
	for (int i = 0; i < no_of_gates; i++){
		if (gates[i].is_output)
			PO_gates.push_back(gates[i]);
	}

	float mu_1, mu_2, sigma_1, sigma_2;
	mu_sigma_struct1 max_strc;
	max_strc.mu = PO_gates[0].delay_mu;
	max_strc.sigma = PO_gates[0].delay_sigma;
	for (int j = 0; j < 2 * no_of_pc; j++){
		max_strc.k_parameters.push_back(k_para_matrix[PO_gates[0].id*offset + j]);
		//max_strc.k_parameters.push_back(PO_gates[0].k_param[j]);
	}

	vector<Gates_cu>::iterator itr = PO_gates.begin();

	for (itr++; itr != PO_gates.end(); itr++){

		mu_1 = max_strc.mu;
		sigma_1 = max_strc.sigma;

		mu_2 = (*itr).delay_mu;
		sigma_2 = (*itr).delay_sigma;


		if (mu_1 - 3 * sigma_1 > mu_2 + 3 * sigma_2)
		{
			continue;
		}

		if (mu_1 + 3 * sigma_1 < mu_2 - 3 * sigma_2)
		{
			max_strc.mu = mu_2;
			max_strc.sigma = sigma_2;
			for (int j = 0; j < 2 * no_of_pc; j++){
				max_strc.k_parameters.push_back(k_para_matrix[(*itr).id*offset + j]);
			}
			continue;
		}


		//step 2
		float co_variance = 0.0f;
		float correlation = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			co_variance += max_strc.k_parameters[j] * k_para_matrix[(*itr).id*offset + j];
		}
		//cout << sigma_1 << " " << sigma_2 << " " << sigma_1*sigma_2 << endl;
		//assert((sigma_1 * sigma_2));
		if (!(sigma_1 * sigma_2)){
			correlation = co_variance/0.5;
		}
		else
			correlation = co_variance / (sigma_1 * sigma_2);
		if (correlation > 0.99 && abs(sigma_1 - sigma_2) < 0.1){
			if (mu_1 > mu_2){
				continue;
			}
			else{
				max_strc.mu = mu_2;
				max_strc.sigma = sigma_2;
				for (int j = 0; j < 2 * no_of_pc; j++){
					max_strc.k_parameters.push_back(k_para_matrix[(*itr).id*offset + j]);
				}
				continue;
			}

		}

		//step 3
		float alpha = sqrt(abs(pow(sigma_1, 2) + pow(sigma_2, 2) - 2 * co_variance));
		if (!(alpha>0))
			cout << alpha << endl;
		float beta = (mu_1 - mu_2) / alpha;

		float phi = pow(2.718281828, -beta*beta / 2) / sqrt(2 * 3.141592654);
		float phi_intg = integrate1(beta);
		float phi_intg_m = integrate1(-beta);

		float sigma_3, mu_3;

		mu_3 = mu_1 * phi_intg + mu_2 * phi_intg_m + alpha * phi;
		float sigma_tmp = (pow(mu_1, 2) + pow(sigma_1, 2)) * phi_intg + (pow(mu_2, 2) + pow(sigma_2, 2)) * phi_intg_m + (mu_1 + mu_2) * alpha * phi - mu_3*mu_3;
		sigma_3 = sqrt(abs(sigma_tmp));

		//step 4
		vector<float> a_r;
		float S0 = 0.0f;
		for (int j = 0; j < 2 * no_of_pc; j++){
			float r_1 = max_strc.k_parameters[j];
			float r_2 = k_para_matrix[(*itr).id*offset + j];
			assert(sigma_3);
			float ar_tmp = (sigma_1*r_1*phi_intg + sigma_2*r_2*phi_intg_m) / sigma_3;
			a_r.push_back(ar_tmp);
			S0 += ar_tmp * ar_tmp;
		}

		for (int j = 0; j < 2 * no_of_pc; j++){
			a_r[j] = a_r[j] * sigma_3 / sqrt(abs(S0));
		}

		max_strc.mu = mu_3;
		max_strc.sigma = sigma_3;
		max_strc.k_parameters = a_r;
	}

	cout << max_strc.mu << " " << max_strc.sigma << endl;
	float phi = 0.0f;
	float yield = 0.0f;
	phi = (req_arr_time - max_strc.mu) / max_strc.sigma;
	cout << phi << endl;
	yield = integrate1(phi);
	cout << "Yield_cu: " << yield << endl;
	
}

void TimingAnalysis::get_max_delay_prob(){

	float max_mu = 0, max_sigma = 0, delay = 0;

	for (int i = 0; i < no_of_gates; ++i){
		if (gates_t[i].is_output){
			if (delay < gates_t[i].delay_mu + 3 * gates_t[i].delay_sigma){
				delay = gates_t[i].delay_mu + 3 * gates_t[i].delay_sigma;
				max_mu = gates_t[i].delay_mu;
				max_sigma = gates_t[i].delay_sigma;
			}
		}

	}
	float phi = (req_arr_time - max_mu) / max_sigma;
	float yield = integrate1(phi);
	cout << "Yield: " << yield << endl;
}

void TimingAnalysis::naive_schedule(){
	vector<Gates* > old_graph;
	vector<Gates* > ready_gates;
	vector<Gates* > new_graph;
	int level = 0;
	int max_fanin = 0;
	

	int i;

	for (int ii = 0; ii < no_of_gates; ++ii){
		i = VECTOR(tmp_nodes)[ii];
		gates_t[i].ready_cnt = 0;
		gates_t[i].in_degree_v = 0;
		gates_t[i].k_v = 0.0f;
		gates_t[i].label = 0;
		gates_t[i].finish = false;
		old_graph.push_back(&gates_t[i]);
	}
	
	if (coff){
		//coff-gra priorty calculate
		vector<int> tmp1;
		vector<int> tmp2; //store the lexi of optimal indx
		int label_cnt = 1;
		for (int i = 0; i < old_graph.size(); ++i)
		{
			if (old_graph[i]->ready_cnt == old_graph[i]->no_of_out){
				ready_gates.push_back(old_graph[i]);
				old_graph.erase(old_graph.begin() + i);
				i--;
			}
		}
		while (!ready_gates.empty()){

			int index = 0;
			tmp2.clear();
			for (int j = ready_gates[0]->start_out; j < ready_gates[0]->start_out + ready_gates[0]->no_of_out; ++j){
				tmp2.push_back(gates_t[edges_t[j]].label);
			}
			if (tmp2.empty())
				tmp2.push_back(0);
			else
				sort(tmp2.begin(), tmp2.end(), myfunction);

			for (int i = 1; i < ready_gates.size(); ++i){
				tmp1.clear();
				for (int j = ready_gates[i]->start_out; j < ready_gates[i]->start_out + ready_gates[i]->no_of_out; ++j){
					tmp1.push_back(gates_t[edges_t[j]].label);
				}
				if (tmp1.empty())
					tmp1.push_back(0);
				else
					sort(tmp1.begin(), tmp1.end(), myfunction);
				int size;
				if (tmp1.size() < tmp2.size())
					size = tmp1.size();
				else
					size = tmp2.size();
				int flag = 1;
				for (int k = 0; k < size; ++k){
					if (tmp1[k] == tmp2[k]){
						if ((k == size - 1) && (tmp1.size() < tmp2.size()))
							flag = 0;
						continue;
					}
					if (tmp1[k] < tmp2[k]){
						flag = 0;
						break;
					}
				}
				if (!flag){
					index = i;
					tmp2 = tmp1;
				}
			}

			ready_gates[index]->label = label_cnt++;
			for (int k = ready_gates[index]->start_in; k < ready_gates[index]->start_in + ready_gates[index]->no_of_in; ++k){
				if (!(gates_t[edges_t[k]].finish)){
					gates_t[edges_t[k]].ready_cnt++;
					if (gates_t[edges_t[k]].ready_cnt == gates_t[edges_t[k]].no_of_out){
						ready_gates.push_back(&gates_t[edges_t[k]]);
						gates_t[edges_t[k]].finish = true;
					}
				}	
			}
			ready_gates.erase(ready_gates.begin() + index);
			/*
			for (int i = 0; i < old_graph.size(); ++i)
			{
				if (old_graph[i]->ready_cnt == old_graph[i]->no_of_out){
					ready_gates.push_back(old_graph[i]);
					old_graph.erase(old_graph.begin() + i);
					i--;
				}
			}
			*/
			
			//cout << "Size: " << old_graph.size()<<endl;
		}


		old_graph.clear();
		ready_gates.clear();
		for (int ii = 0; ii < no_of_gates; ++ii){
			i = VECTOR(tmp_nodes)[ii];
			gates_t[i].ready_cnt = 0;
			gates_t[i].in_degree_v = 0;
			gates_t[i].k_v = 0.0f;
			old_graph.push_back(&gates_t[i]);
		}
	}
	
	assert(old_graph.size() == no_of_gates);
	int cc = 0;
	for (int i = 0; i < old_graph.size(); ++i)
	{
		if (old_graph[i]->ready_cnt == old_graph[i]->no_of_in){
			cc++;
			ready_gates.push_back(old_graph[i]);
			old_graph.erase(old_graph.begin() + i);
			i--;
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	while (!ready_gates.empty()){

		if (ready_gates.size() < MAX_TILE){
			for (int i = 0; i < ready_gates.size(); i++){
				ready_gates[i]->level = level;
				for (int j = ready_gates[i]->start_out; j < ready_gates[i]->start_out + ready_gates[i]->no_of_out; ++j)
					gates_t[edges_t[j]].ready_cnt++;
			}
			new_graph.insert(new_graph.end(), ready_gates.begin(), ready_gates.end());
			ready_gates.clear();
			level++;
		}
		else{
			if (coff){
				sort(ready_gates.begin(), ready_gates.end(), LabelSort);
				for (int i = 0; i < ready_gates.size(); ++i){
				}
				int i_begin = 0, i_end = MAX_TILE;
				for (int i = i_begin; i < i_end; ++i){
					ready_gates[i]->level = level;
					//cout << ready_gates[i]->label << endl;
					for (int j = ready_gates[i]->start_out; j < ready_gates[i]->start_out + ready_gates[i]->no_of_out; ++j)
						gates_t[edges_t[j]].ready_cnt++;
				}
				new_graph.insert(new_graph.end(), ready_gates.begin() + i_begin, ready_gates.begin() + i_end);
				ready_gates.erase(ready_gates.begin() + i_begin, ready_gates.begin() + i_end);
				int t2 = ready_gates.size();
				level++;
			}
			else {
				max_fanin = 0;
				for (int i = 0; i < ready_gates.size(); ++i){
					if (ready_gates[i]->no_of_in > max_fanin)
						max_fanin = ready_gates[i]->no_of_in;
				}

				if (max_fanin)
					for (int i = 0; i < ready_gates.size(); ++i)
						ready_gates[i]->k_v = ready_gates[i]->no_of_in / max_fanin;
				sort(ready_gates.begin(), ready_gates.end(), KValueSort);
				int i_begin, i_end;
				float min_k_diff = 10000.0f;
				for (int i = 0, j = MAX_TILE; j < ready_gates.size(); ++i, ++j){
					if (min_k_diff >(ready_gates[j]->k_v - ready_gates[i]->k_v))
					{
						min_k_diff = ready_gates[j]->k_v - ready_gates[i]->k_v;
						i_begin = i; i_end = j;
					}
				}

				for (int i = i_begin; i < i_end; ++i){
					ready_gates[i]->level = level;
					for (int j = ready_gates[i]->start_out; j < ready_gates[i]->start_out + ready_gates[i]->no_of_out; ++j)
						gates_t[edges_t[j]].ready_cnt++;
				}
				new_graph.insert(new_graph.end(), ready_gates.begin() + i_begin, ready_gates.begin() + i_end);
				ready_gates.erase(ready_gates.begin() + i_begin, ready_gates.begin() + i_end);
				level++;
			}
		}

		for (int i = 0; i < old_graph.size(); ++i)
		{
			if (old_graph[i]->ready_cnt == old_graph[i]->no_of_in){
				cc++;
				ready_gates.push_back(old_graph[i]);
				old_graph.erase(old_graph.begin() + i);
				i--;
			}
		}
	}
	
	

	assert(new_graph.size() == no_of_gates);
	tmp_nodes2 = (int*)malloc(no_of_gates*sizeof(int));
	for (int i = 0; i < new_graph.size(); ++i){
		tmp_nodes2[i] = new_graph[i]->id;
	}



	/*
	cout << new_graph.size() << " " << no_of_gates << " " << cc << endl;
	int l_count = 1;
	int insuff = 0;
	for (int i = 1; i < no_of_gates; ++i)
	{
		l_count++;
		assert(new_graph[i]->level >= new_graph[i - 1]->level);
		if (new_graph[i]->level != new_graph[i - 1]->level)
		{
			cout << "Level: " << new_graph[i-1]->level << " has " << l_count << " gates " << endl;
			if (l_count < 500)
				insuff++;
			l_count = 1;
		}
	}
	cout << "Level: " << new_graph[no_of_gates - 1]->level << " has " << l_count << " gates " << endl;
	cout << "Insufficient group: " << insuff << endl;

	/*
	for (int i = 0; i < new_graph.size(); ++i){
		for (int j = new_graph[i]->start_out; j < new_graph[i]->start_out + new_graph[i]->no_of_out; ++j){
			if (gates_t[j].level >)
		}
	}
	*/
}

void TimingAnalysis::backup(){

#if CL
	gates_cl = (Gates2*)malloc(no_of_gates * sizeof(Gates2));
	for (int i = 0; i < no_of_gates; i++){
		gates_cl[i].adapt_blocks_id = gates_t[i].adapt_blocks_id;
		gates_cl[i].delay = gates_t[i].delay;
		gates_cl[i].delay_mu = gates_t[i].delay_mu;
		gates_cl[i].delay_sigma = gates_t[i].delay_sigma;
		gates_cl[i].gate_capacitance = gates_t[i].gate_capacitance;
		gates_cl[i].gate_mu = gates_t[i].gate_mu;
		gates_cl[i].gate_resistance = gates_t[i].gate_resistance;
		gates_cl[i].gate_sigma = gates_t[i].gate_sigma;
		gates_cl[i].id = gates_t[i].id;
		gates_cl[i].is_input = gates_t[i].is_input;
		gates_cl[i].is_output = gates_t[i].is_output;
		gates_cl[i].level = gates_t[i].level;
		gates_cl[i].load_capacitance = gates_t[i].load_capacitance;
		gates_cl[i].name = gates_t[i].name;
		gates_cl[i].no_of_in = gates_t[i].no_of_in;
		gates_cl[i].no_of_k = gates_t[i].no_of_k;
		gates_cl[i].no_of_out = gates_t[i].no_of_out;
		gates_cl[i].start_in = gates_t[i].start_in;
		gates_cl[i].start_out = gates_t[i].start_out;
		gates_cl[i].type = gates_t[i].type;
		gates_cl[i].x = gates_t[i].y;
		gates_cl[i].y = gates_t[i].y;

	}
#endif
	for (int i = 0; i < no_of_gates; i++){
		gates_t[i].selfV = 0;
		gates_t[i].getV = 0;
	}
	

	for (int i = 0; i < no_of_gates; ++i){
		int ii = tmp_nodes2[i];
		//gates_t[ii].id = i;
		for (int j = 0; j < 2 * no_of_pc; ++j){
			k_para_matrix[i * 2 * no_of_pc + j] = gates_t[ii].k_param[j];	
		}
	}

	gates_cu = (Gates_cu*)malloc(no_of_gates * sizeof(Gates_cu));
	gates_cu2 = (Gates_cu2*)malloc(no_of_gates * sizeof(Gates_cu2));
	for (int ii = 0; ii < no_of_gates; ii++){
		int i = tmp_nodes2[ii];
		gates_cu[ii].id = gates_t[i].id;
		gates_cu[ii].is_input = gates_t[i].is_input;
		gates_cu[ii].is_output = gates_t[i].is_output;
		gates_cu[ii].level = gates_t[i].level;
		gates_cu[ii].delay = gates_t[i].delay;
		gates_cu[ii].gate_mu = gates_t[i].gate_mu;
		gates_cu[ii].gate_sigma = gates_t[i].gate_sigma;
		gates_cu[ii].adapt_blocks_id = gates_t[i].adapt_blocks_id;
		gates_cu[ii].x = gates_t[i].x;
		gates_cu[ii].y = gates_t[i].y;
		gates_cu[ii].output_time_of_gate = gates_t[i].output_time_of_gate;
	}


	cout << "Generate tmp index........." << endl;
	index_tmp = (int*)malloc(no_of_gates*sizeof(int));
	for (int i = 0; i<no_of_gates; i++)
		for (int j = 0; j<no_of_gates; j++)
		{
		if (gates_cu[j].id == i){
			index_tmp[i] = j;
			break;
		}
	}

	igraph_vector_t nei_gates_in;
	igraph_vector_init(&nei_gates_in, 0);
	igraph_vector_t nei_gates_out;
	igraph_vector_init(&nei_gates_out, 0);
	int counter = 0;
	edges_cu = (int*)malloc(no_of_edges*sizeof(int));
	for (int i = 0; i<no_of_gates; i++)
	{
		/*Generate adjlist*/
		int v = igraph_neighbors(&circuitGraph, &nei_gates_in, gates_cu[i].id, IGRAPH_IN);
		assert(!v);
		v = igraph_neighbors(&circuitGraph, &nei_gates_out, gates_cu[i].id, IGRAPH_OUT);
		assert(!v);
		int n1 = igraph_vector_size(&nei_gates_in);
		int n2 = igraph_vector_size(&nei_gates_out);
		gates_cu[i].no_of_in = n1;
		gates_cu[i].no_of_out = n2;
		gates_cu[i].start_in = counter;

		for (int j = 0; j<n1; j++)
		{
			int m = VECTOR(nei_gates_in)[j];
			edges_cu[counter++] = index_tmp[m];
		}
		gates_cu[i].start_out = counter;
		for (int k = 0; k<n2; k++)
		{
			int n = VECTOR(nei_gates_out)[k];
			edges_cu[counter++] = index_tmp[n];
		}

		//////////////////
		gates_cu[i].id = i;
	}
	igraph_vector_destroy(&nei_gates_in);
	igraph_vector_destroy(&nei_gates_out);


	max_level = gates_cu[no_of_gates - 1].level + 1;
	l_count = (int*)malloc(max_level*sizeof(int));
	l_start = (int*)malloc(max_level*sizeof(int));
	l_start[0] = 0;
	int lcount_tmp = 1;
	for (int ii = 1; ii < no_of_gates; ++ii){
		lcount_tmp++;
		if (gates_cu[ii].level != gates_cu[ii - 1].level){
			l_start[gates_cu[ii].level] = ii;
			l_count[gates_cu[ii].level - 1] = lcount_tmp - 1;
			lcount_tmp = 1;
		}
	}
	l_count[max_level - 1] = lcount_tmp;
	cout << max_level << endl;

	for (int ii = 0; ii < no_of_gates; ii++){
		gates_cu2[ii].is_input = gates_cu[ii].is_input;
		gates_cu2[ii].start_in = gates_cu[ii].start_in;
		gates_cu2[ii].no_of_in = gates_cu[ii].no_of_in;
		gates_cu2[ii].delay = gates_cu[ii].delay;
		gates_cu2[ii].x = gates_cu[ii].x;
		gates_cu2[ii].y = gates_cu[ii].y;
		gates_cu2[ii].output_time_of_gate = gates_cu[ii].output_time_of_gate;
	}

}


void TimingAnalysis::prede_test(){

	int ii;
	for (int i = 0; i < no_of_gates; i++){
		ii = tmp_nodes2[i];
		gates_t[ii].selfV = gates_t[ii].no_of_in;
		for (int j = gates_t[ii].start_in; j < gates_t[ii].start_in + gates_t[ii].no_of_in; j++)
			gates_t[ii].getV += gates_t[edges_t[j]].selfV;
	}

}

void TimingAnalysis::Monte_Carlo_Gen(){

	static int init = 0;
	if (init == 0){

		mat_class->random_numbers = (float**)malloc(mat_class->x_grids*mat_class->y_grids*sizeof(float*));
		mat_class->random_numbers_width = (float**)malloc(mat_class->x_grids*mat_class->y_grids*sizeof(float*));

		for (int allocate = 0; allocate < mat_class->x_grids*mat_class->y_grids; allocate++){
			mat_class->random_numbers[allocate] = gaussian_random_numbers1(no_of_iterations);
			mat_class->random_numbers_width[allocate] = gaussian_random_numbers1(no_of_iterations);
		}

		float **chol;
		chol = (float**)malloc(mat_class->x_grids*mat_class->y_grids*sizeof(float*));


		for (int allocate = 0; allocate < mat_class->x_grids*mat_class->y_grids; allocate++){
			//*(chol+allocate) = (double*)malloc(corr->x_grids*corr->y_grids*sizeof(double));
			chol[allocate] = (float*)malloc(mat_class->x_grids*mat_class->y_grids*sizeof(float));
			for (int init = 0; init < mat_class->x_grids*mat_class->y_grids; init++){
				chol[allocate][init] = 0.0;
			}
		}

		// This one for length
		cholesky(mat_class->MC_matrix, mat_class->x_grids*mat_class->y_grids, NULL, 0, chol, NULL, 0);
		matrix_multiplication(mat_class->random_numbers, chol, mat_class->x_grids*mat_class->y_grids, no_of_iterations);

		// This one for width
		cholesky(mat_class->MC_matrix, mat_class->x_grids*mat_class->y_grids, NULL, 0, chol, NULL, 0);
		matrix_multiplication(mat_class->random_numbers_width, chol, mat_class->x_grids*mat_class->y_grids, no_of_iterations);

		init = 1;
	}
}

void TimingAnalysis::Monte_Carlo(covmat_class *mat_class, bool adaptive, int blocks){
	
	Gates* temp_gate;
	float bin_low = 0;
	float bin_high = 0;
	int   row, column;
	int   random_number;
	int	  idx;
	/*
	for (int ii = 0; ii < no_of_gates; ii++)
	{
		idx = tmp_nodes2[ii];
		cout << ii << " " <<gates_t[idx].is_input <<" "<<gates_t[idx].is_output<< endl;
	}
	getchar();
	*/
	double max_arr_avg = 0;
	for (int i = 0; i < no_of_iterations; i++){

		for (int j = 0; j < no_of_gates; ++j)
			gates_t[j].output_time_of_gate = 0.0;



		for (int ii = 0; ii < no_of_gates; ii++)
		{
			idx = tmp_nodes2[ii];
			row = floor(gates_t[idx].y / GRID_SIZE);
			column = floor(gates_t[idx].x / GRID_SIZE);
			random_number = row*(mat_class->x_grids) + column;

			float temp_max = -10000;
			if (gates_t[idx].is_input){
				temp_max = 0;
			}

			for (int in = gates_t[idx].start_in; in < gates_t[idx].start_in + gates_t[idx].no_of_in; ++in){
				if (gates_t[edges_t[in]].output_time_of_gate > temp_max){
					temp_max = gates_t[edges_t[in]].output_time_of_gate;
				}
			}

			gates_t[idx].output_time_of_gate = temp_max + gates_t[idx].delay + (0.15 / 3)*gates_t[idx].delay*(mat_class->random_numbers[random_number][i]) \
												- (0.08/3)*(mat_class->random_numbers_width[random_number][i])*gates_t[idx].delay;
		}

		float max_arr = 0;
		for (int i = 0; i < no_of_gates; ++i){
			if (max_arr < gates_t[i].output_time_of_gate)
				max_arr = gates_t[i].output_time_of_gate;
		}


		if (max_arr > 100000 || max_arr < -100000){
			cout << max_arr << endl;
			continue;
		
		}
		max_arr_avg = max_arr_avg + (double)max_arr;
		max_arr_avg = max_arr_avg / 2;
		if (max_arr < req_arr_time)
			bin_low++;
		else
			bin_high++;
	}

	//max_arr_avg = max_arr_avg / no_of_iterations;
	cout << max_arr_avg << endl;
	cout << "Yield of MC: " << bin_low / no_of_iterations << endl;


}