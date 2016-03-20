#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "Gates.h"
#include <stdio.h>
// Gates_cu struct contains all the computation parameters
#define GRID_SIZE 900
#define MAX_TILE 1024
#define UNROLL 64
#define CASES 768  


__global__ void MC(Gates_cu2** gates, int no_of_gates, int* edges, float* random_numbers, \
	float* random_numbers_width, int* l_count, int* l_start, int current_level, int no_of_iterations, int grids, int x_grid){


	Gates_cu2 tmp;

	int i = threadIdx.x;// +blockIdx.x*blockDim.x;
	int idx = threadIdx.x;
	int b = blockIdx.x;
	float temp_max = -10000;

	int count; int start;

	count = l_count[current_level]; //tex1Dfetch(count_ter, current_level);
	start = l_start[current_level];// tex1Dfetch(start_ter, current_level);


	if (i < count){

		tmp = gates[b][start + i];
		int row = (int)floor((float)(tmp.y / GRID_SIZE));
		int column = (int)floor((float)(tmp.x / GRID_SIZE));
		int random_number = row*x_grid + column;

		if (tmp.is_input)
		{
			temp_max = 0;
		}
		else{

			for (int in = tmp.start_in; in < tmp.start_in + tmp.no_of_in; ++in){
				if (gates[b][edges[in]].output_time_of_gate > temp_max){
					temp_max = gates[b][edges[in]].output_time_of_gate;
				}
			}

			tmp.output_time_of_gate = temp_max + (0.15 / 3)*tmp.delay + tmp.delay*(random_numbers[random_number*no_of_iterations+b]) \
				- (0.08 / 3)*(random_numbers_width[random_number*no_of_iterations + b])*tmp.delay;
		}
		gates[b][start + i] = tmp;

	}

}

extern "C" void cuMC(Gates_cu2* gates, int no_of_gates, int* edges, int no_of_edges, int no_of_iterations, float** rad, float** rad_w, \
	int grids, int x_grid, int* l_count, int* l_start, int max_level){

	cudaError_t cudaStatus;


	Gates_cu2* gates_cu1;
	Gates_cu2** gates_cu;
	
	if (CASES == 1){

		cudaMallocManaged(&gates_cu1, no_of_gates*sizeof(Gates_cu2));
		for (int j = 0; j < no_of_gates; ++j){
			gates_cu1[j] = gates[j];
		}
	}
	else{

		cudaMallocManaged(&gates_cu, CASES * sizeof(Gates_cu2*));
		for (int i = 0; i < CASES; ++i){
			cudaMallocManaged(&gates_cu[i], no_of_gates * sizeof(Gates_cu2));
		}


		for (int i = 0; i < CASES; ++i){
			for (int j = 0; j < no_of_gates; ++j){
				gates_cu[i][j] = gates[j];
			}
		}
	}

	int* edges_cu;

	cudaMallocManaged(&edges_cu, no_of_edges * sizeof(int));
	for (int i = 0; i < no_of_edges; ++i){
		edges_cu[i] = edges[i];
	}



	float* random_numbers;
	cudaMallocManaged(&random_numbers, sizeof(float) * grids * no_of_iterations);
	int count = 0;
	for (int i = 0; i < grids; ++i){
		for (int j = 0; j < no_of_iterations; ++j){
			random_numbers[count++] = rad[i][j];
		}
	}


	float* random_numbers_width;
	cudaMallocManaged(&random_numbers_width, sizeof(float) * grids * no_of_iterations);
	count = 0;
	for (int i = 0; i < grids; ++i){
		for (int j = 0; j < no_of_iterations; ++j){
			random_numbers_width[count++] = rad_w[i][j];
		}
	}

	int* c_l_count;
	cudaMallocManaged(&c_l_count, sizeof(int) * max_level);
	for (int i = 0; i < max_level; ++i)
		c_l_count[i] = l_count[i];

	int* c_l_start;
	cudaMallocManaged(&c_l_start, sizeof(int) * max_level);
	for (int i = 0; i < max_level; ++i)
		c_l_start[i] = l_start[i];



	float time_elapsed = 0;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	int current_level = 0;
	while (current_level < max_level){
		MC << <CASES, MAX_TILE >> > (gates_cu, no_of_gates, edges_cu, random_numbers, random_numbers_width,\
			c_l_count, c_l_start, current_level, no_of_iterations, grids，x_grid);
		cudaDeviceSynchronize();

		current_level++;
	}

	cudaEventRecord(stop, 0);

	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time_elapsed, start, stop);
	printf("GPU MC running：%f(ms)\n", time_elapsed);

	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
	{
		fprintf(stderr, "!!!! GPU program execution error in 2: cuda Error=%d,(%s)\n", cudaStatus, cudaGetErrorString(cudaStatus));
	}





}