#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Gates.h"
#include <stdio.h>
// Gates_cu struct contains all the computation parameters
#define GRID_SIZE 900
#define MAX_TILE 1024
#define UNROLL 64
#define CASES 1  
#ifdef __cplusplus
extern "C" {
#endif
	extern __device__ float pdf_m[36][10] = {
			{ 0.0000, 0.0040, 0.0080, 0.0120, 0.0160, 0.0199, 0.0239, 0.0279, 0.0319, 0.0359 },
			{ 0.0398, 0.0438, 0.0478, 0.0517, 0.0557, 0.0596, 0.0636, 0.0675, 0.0714, 0.0753 },
			{ 0.0793, 0.0832, 0.0871, 0.0910, 0.0948, 0.0987, 0.1026, 0.1064, 0.1103, 0.1141 },
			{ 0.1179, 0.1217, 0.1255, 0.1293, 0.1331, 0.1368, 0.1406, 0.1443, 0.1480, 0.1517 },
			{ 0.1554, 0.1591, 0.1628, 0.1664, 0.1700, 0.1736, 0.1772, 0.1808, 0.1844, 0.1879 },
			{ 0.1915, 0.1950, 0.1985, 0.2019, 0.2054, 0.2088, 0.2123, 0.2157, 0.2190, 0.2224 },
			{ 0.2257, 0.2291, 0.2324, 0.2357, 0.2389, 0.2422, 0.2454, 0.2486, 0.2517, 0.2549 },
			{ 0.2580, 0.2611, 0.2642, 0.2673, 0.2704, 0.2734, 0.2764, 0.2794, 0.2823, 0.2852 },
			{ 0.2881, 0.2910, 0.2939, 0.2967, 0.2995, 0.3023, 0.3051, 0.3078, 0.3106, 0.3133 },
			{ 0.3159, 0.3186, 0.3212, 0.3238, 0.3264, 0.3289, 0.3315, 0.3340, 0.3365, 0.3389 },
			{ 0.3413, 0.3438, 0.3461, 0.3485, 0.3508, 0.3531, 0.3554, 0.3577, 0.3599, 0.3621 },
			{ 0.3643, 0.3665, 0.3686, 0.3708, 0.3729, 0.3749, 0.3770, 0.3790, 0.3810, 0.3830 },
			{ 0.3849, 0.3869, 0.3888, 0.3907, 0.3925, 0.3944, 0.3962, 0.3980, 0.3997, 0.4015 },
			{ 0.4032, 0.4049, 0.4066, 0.4082, 0.4099, 0.4115, 0.4131, 0.4147, 0.4162, 0.4177 },
			{ 0.4192, 0.4207, 0.4222, 0.4236, 0.4251, 0.4265, 0.4279, 0.4292, 0.4306, 0.4319 },
			{ 0.4332, 0.4345, 0.4357, 0.4370, 0.4382, 0.4394, 0.4406, 0.4418, 0.4429, 0.4441 },
			{ 0.4452, 0.4463, 0.4474, 0.4484, 0.4495, 0.4505, 0.4515, 0.4525, 0.4535, 0.4545 },
			{ 0.4554, 0.4564, 0.4573, 0.4582, 0.4591, 0.4599, 0.4608, 0.4616, 0.4625, 0.4633 },
			{ 0.4641, 0.4649, 0.4656, 0.4664, 0.4671, 0.4678, 0.4686, 0.4693, 0.4699, 0.4706 },
			{ 0.4713, 0.4719, 0.4726, 0.4732, 0.4738, 0.4744, 0.4750, 0.4756, 0.4761, 0.4767 },
			{ 0.4772, 0.4778, 0.4783, 0.4788, 0.4793, 0.4798, 0.4803, 0.4808, 0.4812, 0.4817 },
			{ 0.4821, 0.4826, 0.4830, 0.4834, 0.4838, 0.4842, 0.4846, 0.4850, 0.4854, 0.4857 },
			{ 0.4861, 0.4864, 0.4868, 0.4871, 0.4875, 0.4878, 0.4881, 0.4884, 0.4887, 0.4890 },
			{ 0.4893, 0.4896, 0.4898, 0.4901, 0.4904, 0.4906, 0.4909, 0.4911, 0.4913, 0.4916 },
			{ 0.4918, 0.4920, 0.4922, 0.4925, 0.4927, 0.4929, 0.4931, 0.4932, 0.4934, 0.4936 },
			{ 0.4938, 0.4940, 0.4941, 0.4943, 0.4945, 0.4946, 0.4948, 0.4949, 0.4951, 0.4952 },
			{ 0.4953, 0.4955, 0.4956, 0.4957, 0.4959, 0.4960, 0.4961, 0.4962, 0.4963, 0.4964 },
			{ 0.4965, 0.4966, 0.4967, 0.4968, 0.4969, 0.4970, 0.4971, 0.4972, 0.4973, 0.4974 },
			{ 0.4974, 0.4975, 0.4976, 0.4977, 0.4977, 0.4978, 0.4979, 0.4979, 0.4980, 0.4981 },
			{ 0.4981, 0.4982, 0.4982, 0.4983, 0.4984, 0.4984, 0.4985, 0.4985, 0.4986, 0.4986 },
			{ 0.4987, 0.4987, 0.4987, 0.4988, 0.4988, 0.4989, 0.4989, 0.4989, 0.4990, 0.4990 },
			{ 0.4990, 0.4991, 0.4991, 0.4991, 0.4992, 0.4992, 0.4992, 0.4992, 0.4993, 0.4993 },
			{ 0.4993, 0.4993, 0.4994, 0.4994, 0.4994, 0.4994, 0.4994, 0.4995, 0.4995, 0.4995 },
			{ 0.4995, 0.4995, 0.4995, 0.4996, 0.4996, 0.4996, 0.4996, 0.4996, 0.4996, 0.4997 },
			{ 0.4997, 0.4997, 0.4997, 0.4997, 0.4997, 0.4997, 0.4997, 0.4997, 0.4997, 0.4998 },
			{ 0.4998, 0.4998, 0.4998, 0.4998, 0.4998, 0.4998, 0.4998, 0.4998, 0.4998, 0.4998 },
	};

	extern __device__ float row[36] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5 };
	extern __device__ float col[10] = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09 };
#ifdef __cplusplus
}
#endif


texture<int, 1, cudaReadModeElementType> count_ter;
texture<int, 1, cudaReadModeElementType> start_ter;
//texture<int, 1, cudaReadModeElementType> edges_ter;
//texture<int, 1, cudaReadModeElementType> sort_ter;
texture<float, 1, cudaReadModeElementType> k_m_ter;

typedef struct {
	float mu;
	float sigma;
}mu_sigma_struct;

__device__ float integrate1(float beta){

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

__device__ float integrate(float beta){

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

	return sum / sqrt(2 * 3.141592654f);

}


__device__ void sum_function(mu_sigma_struct max_strc, Gates_cu* currentG, float* k_para_matrix, int no_of_pc, \
	int max_k_idx){

	currentG->delay_mu = currentG->gate_mu + max_strc.mu;
	int offset = 2 * no_of_pc;
	int g_id_tmp = currentG->id;
	float tmp = 0.0f;
	float tmp2 = 0.0f;
	float sigma = 0.0f;

	for (int i = 0; i < offset; i++)
	{
		//tmp = tex1Dfetch(k_m_ter, (g_id_tmp * 2 * no_of_pc + i)) + tex1Dfetch(k_m_ter, (max_k_idx * 2 * no_of_pc + i));
		tmp = k_para_matrix[g_id_tmp * 2 * no_of_pc + i] + k_para_matrix[max_k_idx * 2 * no_of_pc + i];
		//tmp =  l_k_m[idx * 40 + i] + l_k_m[idx * 40 + i];
		sigma += pow(tmp, 2);
	}

	currentG->delay_sigma = sqrt(sigma);

}


__device__ void max_function(Gates_cu* gates_t, Gates_cu* currentG, int* edges_t, int no_of_pc, float* k_para_matrix,  \
	int bidx, int current_level){

	int e_idx = edges_t[currentG->start_in];
	float mu_1, mu_2, sigma_1, sigma_2;
	int offset = 2 * no_of_pc;
	int max_k_idx;
	mu_sigma_struct max_strc;
	int current_k_idx;
	max_strc.mu = gates_t[e_idx].delay_mu;
	max_strc.sigma = gates_t[e_idx].delay_sigma;
	max_k_idx = gates_t[e_idx].id;

	for (int i = currentG->start_in + 1; i < currentG->start_in + currentG->no_of_in; i++){
		
		e_idx = edges_t[i];
		mu_1 = max_strc.mu;
		sigma_1 = max_strc.sigma;
		mu_2 = gates_t[e_idx].delay_mu;
		sigma_2 = gates_t[e_idx].delay_sigma;
		current_k_idx = gates_t[e_idx].id;

		
		if (mu_1 - 3 * sigma_1 > mu_2 + 3 * sigma_2)
		{
			continue;
		}

		if (mu_1 + 3 * sigma_1 < mu_2 - 3 * sigma_2)
		{
			max_strc.mu = mu_2;
			max_strc.sigma = sigma_2;
			max_k_idx = current_k_idx;
			continue;
		}
		
		//step 2
		float co_variance = 0.0f;
		float correlation = 0.0f;

		for (int j = 0; j < offset; j++){
			//co_variance += tex1Dfetch(k_m_ter, (max_k_idx * offset + j)) * tex1Dfetch(k_m_ter, (current_k_idx * offset + j));
			co_variance += k_para_matrix[max_k_idx * offset + j] * k_para_matrix[current_k_idx * offset + j];
			//co_variance += k_l[fake + j] * k_l[fake + j];
		}

		correlation = co_variance / (sigma_1 * sigma_2);
		if (correlation > 0.99 && abs(sigma_1 - sigma_2) < 0.1){
			if (mu_1 > mu_2){
				continue;
			}
			else{
				max_strc.mu = mu_2;
				max_strc.sigma = sigma_2;
				max_k_idx = current_k_idx;
				continue;
			}

		}

		//step 3
		
		float alpha = sqrt(abs(pow(sigma_1, 2) + pow(sigma_2, 2) - 2 * co_variance));
		float beta = (mu_1 - mu_2) / alpha;

		float phi = pow(2.718281828f, -beta*beta / 2) / sqrt(2 * 3.141592654f);
		float phi_intg = integrate1(beta);
		float phi_intg_m = integrate1(-beta);

		float sigma_3, mu_3;

		mu_3 = mu_1 * phi_intg + mu_2 * phi_intg_m + alpha * phi;
		float sigma_tmp = (pow(mu_1, 2) + pow(sigma_1, 2)) * phi_intg + (pow(mu_2, 2) + pow(sigma_2, 2)) * phi_intg_m + (mu_1 + mu_2) * alpha * phi - mu_3*mu_3;
		sigma_3 = sqrt(abs(sigma_tmp));
		/*
		step 4
		float S0 = 0.0f;
		for (int j = 0; j < offset; j++){
		float r_1 = tex2D(k_m_ter, max_k_idx, j);
		float r_2 = tex2D(k_m_ter, current_k_idx, j);
		float ar_tmp = (sigma_1*r_1*phi_intg + sigma_2*r_2*phi_intg_m) / sigma_3;
		max_k[j] = ar_tmp;
		S0 += ar_tmp * ar_tmp;
		}


		for (int j = 0; j < offset; j++){
		max_k[j] = max_k[j] * sigma_3 / sqrt(abs(S0));
		}
		*/
		max_strc.mu = mu_3;
		max_strc.sigma = sigma_3;
		

	}

	sum_function(max_strc, currentG, k_para_matrix, no_of_pc, max_k_idx);


}


__global__ void PCA(Gates_cu* gates, int no_of_gates, int* sort, int* edges, float* k_param, \
	int* l_count, int* l_start, int current_level, int no_of_pc, int max_level){

	Gates_cu tmp;
	//current_level=0;
	int i = threadIdx.x;// +blockIdx.x*blockDim.x;
	int idx = threadIdx.x;
	int b = blockIdx.x;

	//while (current_level < max_level){

		int count; int start;

		count = l_count[current_level];
		start = l_start[current_level];


		if (i < count){

			tmp = gates[start + i];

			if (tmp.is_input)
			{
				tmp.delay_mu = tmp.gate_mu;
				tmp.delay_sigma = tmp.gate_sigma;
			}
			else{
				//--max--//
				max_function(gates, &tmp, edges, no_of_pc, k_param, b, current_level); //no change supposed to be on gates[i]
				//printf("%f\n", l_gates[idx].delay_sigma);
				//--sum--//
				//sum_function(max_strc, &l_gates[i], k_param, no_of_pc, gates, edges); // gates[i] param should changes
			}
			//l_gates[tmp.id%MAX_TILE] = tmp;
			gates[start + i] = tmp;
			//test(gates, &l_gates[i], edges);
			//gates[b][sort[start + i]] = l_gates[idx];
			//gates[sort[start + i]] = tmp;
		}
		//current_level++;
		//__syncthreads();
	//}
}

__global__ void timing(Gates_cu* gates, int no_of_pc, float* eigen_v, float** eigen_vec, int x_grids, int no_of_gates)
{
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	if (k < no_of_gates){
		float sigma_of_delay = 0.0f;
		gates[k].gate_mu = 0.0f;  // This is to init mu
		gates[k].gate_sigma = 0.0f; // This is to init sigma
		gates[k].gate_mu = gates[k].delay;
		int row = (int)floor((float)(gates[k].y / GRID_SIZE));
		int column = (int)floor((float)(gates[k].x / GRID_SIZE));

		int i_of_j = row*x_grids + column;

		float k_tmp = 0.0f;


		for (int j = 0; j < no_of_pc; j++)
		{
			if (eigen_v[j] < 0){
				eigen_v[j] = 0;
			}

			//L = sqr eg_v * eg_vec * sigma, the dRdL is the constain for a specific size of gate, i'll update it
			k_tmp = (0.15f / 3)* gates[k].delay * sqrt(eigen_v[j]) * eigen_vec[i_of_j][j];// *sigma_of_L;
			sigma_of_delay += k_tmp * k_tmp;
			//gates[k].k_param[j] = k_tmp;

		}


		for (int j = 0; j < no_of_pc; j++)
		{
			//W = sqr eg_v * eg_vec * sigma, the dRdW is the constain for a specific size of gate, i'll update it
			k_tmp = -(0.08f / 3)* gates[k].delay * sqrt(eigen_v[j]) * eigen_vec[i_of_j][j];// *sigma_of_W;
			sigma_of_delay += k_tmp * k_tmp;
			//gates[k].k_param[j + no_of_pc] = k_tmp;
		}

		/*----get gate[i]'s sigma----*/

		gates[k].gate_sigma = sqrt(sigma_of_delay);
	}
}

extern "C" Gates_cu* cuSSTA(Gates_cu* gates, int no_of_gates, int* sort, int* edges, int no_of_edges, int no_of_pc, float *eigen_values, \
	float **eigen_vectors, int x_grid, float* k_param, int* l_count, int* l_start, int max_level, Gates* gates_t){

	cudaError_t cudaStatus;

	int* c_sort;
	cudaMallocManaged(&c_sort, sizeof(int) * no_of_gates);
	for (int i = 0; i < no_of_gates; ++i){
		c_sort[i] = sort[i];
	}


	float* k_parameters;
	cudaMallocManaged(&k_parameters, sizeof(float) * no_of_gates * no_of_pc * 2);
	for (int i = 0; i < no_of_gates * no_of_pc * 2; ++i){
		k_parameters[i] = k_param[i];
	}


	int* c_l_count;
	cudaMallocManaged(&c_l_count, sizeof(int) * max_level);
	for (int i = 0; i < max_level; ++i)
		c_l_count[i] = l_count[i];

	int* c_l_start;
	cudaMallocManaged(&c_l_start, sizeof(int) * max_level);
	for (int i = 0; i < max_level; ++i)
		c_l_start[i] = l_start[i];

	Gates_cu* gates_cu1;
	Gates_cu** gates_cu;
	Gates_cu* g_str;
	if (CASES == 1){

		cudaMallocManaged(&gates_cu1, no_of_gates*sizeof(Gates_cu));
		cudaMallocManaged(&g_str, no_of_gates*sizeof(Gates_cu));
		for (int j = 0; j < no_of_gates; ++j){
			gates_cu1[j] = gates[j];
		}
	}
	else{

		cudaMallocManaged(&gates_cu, CASES * sizeof(Gates_cu*));
		for (int i = 0; i < CASES; ++i){
			cudaMallocManaged(&gates_cu[i], no_of_gates * sizeof(Gates_cu));
		}


		for (int i = 0; i < CASES; ++i){
			for (int j = 0; j < no_of_gates; ++j){
				gates_cu[i][j] = gates[j];
			}
		}
	}

	float* eigen_v;
	float** eigen_vec;
	int* edges_cu;

	cudaMallocManaged(&edges_cu, no_of_edges * sizeof(int));
	cudaMallocManaged(&eigen_v, no_of_pc * sizeof(float));
	cudaMallocManaged(&eigen_vec, no_of_pc * sizeof(float*));

	for (int i = 0; i < no_of_edges; ++i){
		edges_cu[i] = edges[i];
	}
	for (int i = 0; i < no_of_pc; ++i){
		eigen_v[i] = eigen_values[i];
	}
	for (int i = 0; i < no_of_pc; ++i){
		cudaMallocManaged(&eigen_vec[i], no_of_pc * sizeof(float));
		for (int j = 0; j < no_of_pc; ++j){
			eigen_vec[i][j] = eigen_vectors[i][j];
		}
	}

	for (int i = 0; i < no_of_pc; ++i){
		assert(eigen_v[i] == eigen_values[i]);
		for (int j = 0; j < no_of_pc; ++j)
			assert(eigen_vectors[i][j] == eigen_vec[i][j]);
	}



	/////////////////////////////////////////////////////////////////////////////////////////
	dim3 blockn = (no_of_gates - 1) / 1024 + 1;
	dim3 threadn = 1024;

	//cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

	if (CASES == 1){
		timing << < blockn, threadn >> > (gates_cu1, no_of_pc, eigen_v, eigen_vec, x_grid, no_of_gates);
		cudaDeviceSynchronize();
	}
	else{
		for (int i = 0; i < CASES; ++i){
			timing << < blockn, threadn >> > (gates_cu[i], no_of_pc, eigen_v, eigen_vec, x_grid, no_of_gates);
			cudaDeviceSynchronize();
		}
	}


	float time_elapsed = 0;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	int current_level = 0;
	while (current_level < max_level){
		PCA << <CASES, MAX_TILE >> > (gates_cu1, no_of_gates, c_sort, edges_cu, k_parameters, \
			c_l_count, c_l_start, current_level, no_of_pc, max_level);
		cudaDeviceSynchronize();
		current_level++;
	}

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_elapsed, start, stop);
	printf("GPU running：%f(ms)\n", time_elapsed);


	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "!!!! GPU program execution error in 2: cuda Error=%d,(%s)\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}


	return gates_cu1;

}
