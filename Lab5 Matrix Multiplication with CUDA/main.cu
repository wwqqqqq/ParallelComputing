#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEFAULT_M 1000
#define DEFAULT_N 1000
#define DEFAULT_L 1000

#define TILE_WIDTH 16

//Total amount of shared memory per block:			49152 bytes
//Maximum number of threads per multiprocessor:		2048
//Maximum number of threads per block :				1024
//Max dimension size of a thread block(x, y, z) :	(1024, 1024, 64)
//Max dimension size of a grid size(x, y, z) :		(2147483647, 65535, 65535)

void multMatrixNoCUDA(const float* a, const float* b, float* c, int m, int n, int l);
__global__  void multMatrixCUDA(const float* a, const float* b, float* c, int m, int n, int l);
__global__  void multMatrixCUDA_tiled(const float* a, const float* b, float* c, int m, int n, int l);
bool InitCUDA();
void generate_matrix(float* mat, int m, int n);

void printMat(const float* mat, int m, int n) {
	printf("\n***********\n");
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++)
			printf("%f\t", mat[i*n + j]);
		printf("\n");
	}
	printf("***********\n");
}



int main(void) {
	int m, n, l;
	m = DEFAULT_M;
	n = DEFAULT_N;
	l = DEFAULT_L;
	//input m, n, and l
	printf("Please input m, n, and l. Input 0 to use default setting.\n");
	int temp = 0;
	scanf("%d", &temp);
	if (temp > 0) {
		m = temp;
		scanf("%d", &temp);
		if (temp > 0) n = temp;
		scanf("%d", &temp);
		if (temp > 0) l = temp;
	}


	if (!InitCUDA()) return 0;

	float *a = (float*)malloc(sizeof(float)*m*n);
	float *b = (float*)malloc(sizeof(float)*n*l);
	float *c = (float*)malloc(sizeof(float)*m*l);
	float *d = (float*)malloc(sizeof(float)*m*l);
	srand((unsigned int)time(NULL));
	generate_matrix(a, m, n);
	generate_matrix(b, n, l);
	clock_t st, ed;

	//CPU
	st = clock();
	multMatrixNoCUDA(a, b, c, m, n, l);
	ed = clock();
	printf("Using CPU time = %lfms\n", (double)(ed - st) / CLOCKS_PER_SEC * 1000);

	//GPU
	float *cuda_a, *cuda_b, *cuda_c, *cuda_d;
	cudaMalloc((void**)&cuda_a, sizeof(float)*m*n);
	cudaMalloc((void**)&cuda_b, sizeof(float)*n*l);
	cudaMalloc((void**)&cuda_c, sizeof(float)*m*l);
	cudaMalloc((void**)&cuda_d, sizeof(float)*m*l);
	cudaMemcpy(cuda_a, a, sizeof(float)*m*n, cudaMemcpyHostToDevice);
	cudaMemcpy(cuda_b, b, sizeof(float)*n*l, cudaMemcpyHostToDevice);
	int nblocks = (m * l + 255) / 256;
	//calculating...
	st = clock();
	multMatrixCUDA <<< nblocks, 256 >>> (cuda_a, cuda_b, cuda_c, m, n, l);
	cudaDeviceSynchronize(); //force CPU to wait until CUDA kernel executing finished
	ed = clock();
	//calculation end

	cudaMemcpy(d, cuda_c, sizeof(float)*m*l, cudaMemcpyDeviceToHost);
	printf("Using CUDA time = %lfms\n", (double)(ed - st) / CLOCKS_PER_SEC * 1000);
	//test correctness of result
	float error = 0;
	float maxerror = 0;
	int counterror = 0;
	for (int i = 0; i < m * l; i++) {
		float temp = (c[i] > d[i]) ? (c[i] - d[i]) : (d[i] - c[i]);
		if (temp > maxerror) maxerror = temp;
		error += temp;
		if (temp > 1e-5) counterror++;
	}
	if (counterror == 0) {
		printf("Result correct!\n");
	}
	printf("average error = %f, max error = %f\n", error / m / l, maxerror);

	//GPU - tiled algorithm
	dim3 gridSize, blockSize;
	blockSize.x = TILE_WIDTH;
	blockSize.y = TILE_WIDTH; 
	blockSize.z = 1;
	gridSize.x = (m + blockSize.x - 1) / blockSize.x; 
	gridSize.y = (l + blockSize.y - 1) / blockSize.y;
	gridSize.z = 1;
	//calculating...
	st = clock();
	multMatrixCUDA_tiled <<<gridSize, blockSize>>>(cuda_a, cuda_b, cuda_d, m, n, l);
	cudaDeviceSynchronize();
	ed = clock();
	//calculation end

	cudaMemcpy(d, cuda_d, sizeof(float)*m*l, cudaMemcpyDeviceToHost);
	printf("Using CUDA time = %lfms (tiled)\n", (double)(ed - st) / CLOCKS_PER_SEC * 1000);
	//test correctness of result
	error = 0;
	maxerror = 0;
	counterror = 0;
	for (int i = 0; i < m * l; i++) {
		float temp = (c[i] > d[i]) ? (c[i] - d[i]) : (d[i] - c[i]);
		if (temp > maxerror) maxerror = temp;
		error += temp;
		if (temp > 1e-5) counterror++;
	}
	if (counterror == 0) {
		printf("Result correct!\n");
	}
	printf("average error = %f, max error = %f\n", error / m / l, maxerror);

	//end

	/*//print matrix
	printMat(a, m, n);
	printMat(b, n, l);
	printMat(c, m, l);
	printMat(d, m, l);
	*/
	free(a);
	free(b);
	free(c);
	free(d);
	cudaFree(cuda_a);
	cudaFree(cuda_b);
	cudaFree(cuda_c);
	cudaFree(cuda_d);

	getchar();
	getchar();
	return 0;
}


bool InitCUDA() {
	int count;
	cudaGetDeviceCount(&count);

	if (count == 0) {
		fprintf(stderr, "There is no device.\n");
		return false;
	}
	int i;
	for (i = 0; i<count; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		//printDeviceProp(prop);
		if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if (prop.major >= 1) break;
		}
	}
	if (i == count) {
		fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
		return false;
	}
	cudaSetDevice(i);
	return true;
}

void generate_matrix(float* mat, int m, int n) {
	//srand((unsigned int)time(NULL));
	for (int i = 0; i<m; i++)
		for (int j = 0; j<n; j++) {
			//float temp = (float)rand() / RAND_MAX;
			float temp = 0;
			if (rand() % 2 == 1) temp = -temp;
			float temp2 = rand() % 100;
			if (rand() % 2 == 1) temp2 = -temp2;
			mat[i*n + j] = temp2 + temp;
		}
}

__global__  void multMatrixCUDA(const float* a, const float* b, float* c, int m, int n, int l) {
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int row = idx / l;
	int column = idx % l;

	if (row < m && column < l) {
		float sum = 0;
		for (int i = 0; i < n; i++) {
			sum += a[row * n + i] * b[i * l + column];
		}
		c[idx] = sum;

	}
}

void multMatrixNoCUDA(const float* a, const float* b, float* c, int m, int n, int l) {
	for (int i = 0; i<m; i++) {
		for (int j = 0; j<l; j++) {
			double sum = 0;
			for (int k = 0; k<n; k++) {
				sum += a[i*n + k] * b[k*l + j];
			}
			c[i*l + j] = sum;
		}
	}
}

__global__  void multMatrixCUDA_tiled(const float* a, const float* b, float* c, int m, int n, int l) {
	__shared__ float shared_a[TILE_WIDTH][TILE_WIDTH];
	__shared__ float shared_b[TILE_WIDTH][TILE_WIDTH];
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int row = bx * TILE_WIDTH + tx;
	int col = by * TILE_WIDTH + ty;
	float sum = 0.0;
	for (int i = 0; i*TILE_WIDTH <= n; i++) {
		if (i * TILE_WIDTH + ty < n && row < m)
			shared_a[tx][ty] = a[row * n + i * TILE_WIDTH + ty];
		else shared_a[tx][ty] = 0;
		if (i * TILE_WIDTH + tx < n && col < l)
			shared_b[tx][ty] = b[(i * TILE_WIDTH + tx)*l + col];
		else shared_b[tx][ty] = 0;
		__syncthreads();
		for (int k = 0; k < TILE_WIDTH; k++)
			sum += shared_a[tx][k] * shared_b[k][ty];
		__syncthreads();
	}
	if (row < m && col < l)
		c[row*l + col] = sum;
}
