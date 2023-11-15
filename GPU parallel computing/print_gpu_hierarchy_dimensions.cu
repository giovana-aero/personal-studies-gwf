#include<stdio.h>

__global__
void print_sizes();

int main(){

	print_sizes<<<2,2>>>();
	
	return 0;
}

__global__
void print_sizes(){
	printf("block indexes: %d %d %d\n",blockIdx.x,blockIdx.y,blockIdx.z);
	printf("thread indexes: %d %d %d\n",threadIdx.x,threadIdx.y,threadIdx.z);
	printf("grid dimensions: %d, %d, %d\n",gridDim.x,gridDim.y,gridDim.z);
	printf("block dimensions: %d, %d, %d\n",blockDim.x,blockDim.y,blockDim.z);
}