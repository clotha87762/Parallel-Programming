#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>


int* d;
int* graph;



__constant__ int cuda_bf;
__constant__ int cuda_total_vertex;
__constant__ int cuda_tempVertex;

#define INF 1e9
#define H2D cudaMemcpyHostToDevice
#define D2H cudaMemcpyDeviceToHost

using namespace std;

int
init_device ()
{	
	cudaSetDevice(0);
	return 0;
}

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
} while (0)

//extern __shared__ int D[];
__global__ void floyd_warshall_1(int* dist,int k ,int kbf){

	int idx,idy;


	idx = k ;
	idy = k ;

	int i = cuda_bf * idx + threadIdx.y;
	int j = cuda_bf * idy + threadIdx.x;
	if(i>=cuda_total_vertex||j>=cuda_total_vertex)
		return ;

	__shared__ int D[32*32];
	D[threadIdx.y*cuda_bf + threadIdx.x] = dist[i*cuda_tempVertex + j];
	__syncthreads();
	// Put to shared memory???
	int x = 0;

	//int dij = dist[i*total_vertex + j];
	//int dik = dist[i*total_vertex + k];
	//int dkj = dist[k*total_vertex + j];
	int dij ,dik,dkj;
	int a = threadIdx.y * cuda_bf + threadIdx.x;
	int b = threadIdx.y * cuda_bf;
	while( x < cuda_bf ){
		dij = D[a];
		dik = D[b + x];
		dkj = D[x*cuda_bf + threadIdx.x];
		if(dij>dik+dkj){
			D[a] = dik + dkj;
		}
		__syncthreads();
		x++;
	}

	dist[i*cuda_tempVertex + j] = D[threadIdx.y*cuda_bf + threadIdx.x];

	return ;
}

__global__ void floyd_warshall_2(int* dist,int k , int kbf ){

	int idx,idy;
	if(blockIdx.x % 2 == 0 ){
		idx = (blockIdx.x/2) >= k ? (blockIdx.x/2+1):(blockIdx.x/2);
		idy = k;
	}
	else {
		idx = k;
		idy = (blockIdx.x/2) >= k ? (blockIdx.x/2+1):(blockIdx.x/2);
	}

	int i = cuda_bf * idx + threadIdx.y;
	int j = cuda_bf * idy + threadIdx.x;
	//bool flag = 0;
	//if(i>=cuda_total_vertex||j>=cuda_total_vertex)
	//	return;

	__shared__ int D2[32*32*2];
	D2[threadIdx.y * cuda_bf + threadIdx.x] = dist[i*cuda_tempVertex + j];
	D2[(cuda_bf*cuda_bf) + (threadIdx.y *cuda_bf ) + (threadIdx.x)] = dist[ (kbf+threadIdx.y) * cuda_tempVertex + (kbf +threadIdx.x)];
	__syncthreads();

	// Put to shared memory???
	int x = 0;

	int dij ,dik,dkj;
	int a = (threadIdx.y * cuda_bf + threadIdx.x);
	int b;
	if(blockIdx.x%2==0){
		b = cuda_bf*cuda_bf + threadIdx.x;
	}
	else{
		b = cuda_bf*cuda_bf + cuda_bf*threadIdx.y;
	}

	dij = D2[a];

	while(x<cuda_bf){

		if(blockIdx.x%2==0){
			dik = D2[cuda_bf*threadIdx.y + x];
			dkj = D2[b + (x*cuda_bf)];
		}
		else{
			dik = D2[b + x];
			dkj = D2[x*cuda_bf + threadIdx.x];
		}
		if(dij>dik+dkj){
			dij = dik + dkj;
		}
		__syncthreads();
		x++;
	}
	dist[i*cuda_tempVertex + j] = dij;

	return ;
}

__global__ void floyd_warshall_3(int* dist, int k ,int kbf){

	int idx,idy;
	idy = blockIdx.y >= k? blockIdx.y + 1 : blockIdx.y;
	idx = blockIdx.x >= k? blockIdx.x + 1 : blockIdx.x;

	int i = cuda_bf * idx + threadIdx.y;
	int j = cuda_bf * idy + threadIdx.x;
	//if(i>=cuda_total_vertex||j>=cuda_total_vertex)
	//	return ;

	__shared__ int D3[32*32*3];
	D3[threadIdx.y * cuda_bf + threadIdx.x] = dist[i*cuda_tempVertex + j];
	D3[(cuda_bf*cuda_bf) + (threadIdx.y*cuda_bf) + threadIdx.x] = dist[(cuda_bf*idx+threadIdx.y)*cuda_tempVertex + (kbf + threadIdx.x)];
	D3[(2*cuda_bf*cuda_bf) + (threadIdx.y*cuda_bf) + threadIdx.x] = dist[(kbf+threadIdx.y)*cuda_tempVertex + (idy*cuda_bf+threadIdx.x)];
	__syncthreads();

	// Put to shared memory???
	int x = 0;
	int dij ,dik,dkj;

	int a =threadIdx.y * cuda_bf + threadIdx.x;
	int b = cuda_bf*cuda_bf  + threadIdx.y*cuda_bf;
	int c = 2*cuda_bf*cuda_bf + threadIdx.x; 
	
	dij = D3[a];

	while(x<cuda_bf){
		dik = D3[b + x];
		dkj = D3[x*cuda_bf + c];
		if(dij>dik+dkj){
			dij = dik + dkj;
		}
		x++;
	}
	dist[i*cuda_tempVertex + j] = dij;

	return ;
}



__global__ void floyd_warshall_beta_1(int* dist, int k , int kbf  ){

	int idx,idy;
	idx = k;
	idy = k;
	int i = cuda_bf * idx + (blockIdx.x%cuda_bf);
	int j = cuda_bf * idy + threadIdx.x;
		if(i>=cuda_total_vertex||j>=cuda_total_vertex)
		return ;	

	// Put to shared memory???

	int dij = dist[i*cuda_tempVertex + j];
	int dik = dist[i*cuda_tempVertex + kbf];
	int dkj = dist[kbf*cuda_tempVertex + j];

	if(dij>dik+dkj){
		dist[i*cuda_tempVertex+j] = dik + dkj;
	}

	return ;
}

__global__ void floyd_warshall_beta_2(int* dist, int k , int kbf  ){

	int idx,idy;
	int temp = blockIdx.x / cuda_bf;
	if( (temp) % 2 == 0 ){
		idx = (temp/2) >= k ? (temp/2+1):(temp/2);
		idy = k;
	}
	else {
		idx = k;
		idy = (temp/2) >= k ? (temp/2+1):(temp/2);
	}

	int i = cuda_bf * idx + (blockIdx.x%cuda_bf);
	int j = cuda_bf * idy + threadIdx.x;
		if(i>=cuda_total_vertex||j>=cuda_total_vertex)
		return ;	

	// Put to shared memory???

	int dij = dist[i*cuda_tempVertex + j];
	int dik = dist[i*cuda_tempVertex + kbf];
	int dkj = dist[kbf*cuda_tempVertex + j];

	if(dij>dik+dkj){
		dist[i*cuda_tempVertex+j] = dik + dkj;
	}

	return ;
}

__global__ void floyd_warshall_beta_3(int* dist, int k , int kbf ,int grid_size ){

	int idx,idy;
	int temp = ((blockIdx.y*gridDim.x) + blockIdx.x) / cuda_bf;
	idx = temp/grid_size  >= k?  temp/grid_size + 1 : temp/grid_size;
	idy = temp % grid_size >= k? temp%grid_size + 1 : temp % grid_size;

	int i = cuda_bf * idx + (blockIdx.x%cuda_bf);
	int j = cuda_bf * idy + threadIdx.x;
		if(i>=cuda_total_vertex||j>=cuda_total_vertex)
		return ;	

	// Put to shared memory???

	int x = kbf + cuda_bf;
	int dij ,dik,dkj;

	while(kbf<x){
		dij = dist[i*cuda_tempVertex + j];
		dik = dist[i*cuda_tempVertex + kbf];
		dkj = dist[kbf*cuda_tempVertex + j];
		if(dij>dik+dkj){
			dist[i*cuda_tempVertex + j] = dik + dkj;
		}
		//__syncthreads();
		kbf++;
	}
	return;
}




int main(int argc,char* argv[]){

	cudaEvent_t total_start, total_stop;
    cudaEvent_t com_start, com_stop;
    cudaEvent_t mem_start, mem_stop;
	cudaEvent_t io_start, io_stop;

	float total_temp=0,total_total=0,io_temp =0 , io_total=0 , com_temp =0,com_total=0 , mem_temp=0 , mem_total=0;

	cudaEventCreate(&total_start);
	cudaEventCreate(&total_stop);
	cudaEventCreate(&com_start);
	cudaEventCreate(&com_stop);
	cudaEventCreate(&mem_start);
	cudaEventCreate(&mem_stop);
	cudaEventCreate(&io_start);
	cudaEventCreate(&io_stop);


	cudaEventRecord(total_start); 

	init_device();
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop,0);
	fprintf(stderr,"clock rate %lf\n",prop.clockRate);

	int bf = atoi(argv[3]);
	int total_vertex;
	int edge_num;
	ifstream input;
	ofstream output;
	input.open(argv[1]);

	input >> total_vertex;
	input >> edge_num;

	int tempVertex = total_vertex % bf ?  (total_vertex + (bf - (total_vertex%bf) )): total_vertex;
	

	//fprintf(stderr,"tempVertex:%d\n",tempVertex);

	d = new int[tempVertex*tempVertex];
	graph = new int[tempVertex*tempVertex];

	for(int i=0;i<tempVertex;i++){
		for(int j=0;j<tempVertex;j++){
			graph[i*tempVertex+j] = INF;
		}
		graph[i*tempVertex + i ]=0;
	}


	cudaEventRecord(io_start);
	for(int i=0;i<edge_num;i++){
		int a,b;
		input >> a;
		input >> b;
		input >> graph[(a-1)*tempVertex + (b-1) ];
		//fprintf(stderr,"graph %d %d :%d\n",a,b,graph[a*tempVertex+b]);
	}
	cudaEventRecord(io_stop);
	cudaEventSynchronize(io_stop);
	cudaEventElapsedTime(&io_temp,io_start,io_stop);
	io_total += io_temp;

	int* cuda_graph;

	cudaMalloc((void**)&cuda_graph,sizeof(int)*tempVertex*tempVertex);
	cudaCheckErrors("malloc gpu");

	cudaEventRecord(mem_start);	
	cudaMemcpy(cuda_graph,graph,sizeof(int)*tempVertex*tempVertex,H2D);
	cudaCheckErrors("memcpy gpu");
	
	cudaMemcpyToSymbol(cuda_bf,&bf,sizeof(int));
	cudaMemcpyToSymbol(cuda_total_vertex,&total_vertex,sizeof(int));
	cudaMemcpyToSymbol(cuda_tempVertex,&tempVertex,sizeof(int));


	cudaEventRecord(mem_stop);
	cudaEventSynchronize(mem_stop);
	cudaEventElapsedTime(&mem_temp,mem_start,mem_stop);
	mem_total += mem_temp;

	//int FWblockDim = total_vertex%bf ? (total_vertex/bf + 1) : total_vertex/bf;
	//int remainBF = total_vertex%bf? total_vertex%bf : bf ;
	int FWblockDim = tempVertex / bf ;

	dim3 threadStr(bf,bf);
	dim3 blockStr(FWblockDim-1,FWblockDim-1);
	dim3 blockStr2((FWblockDim-1)*bf,FWblockDim-1);
	cudaEventRecord(com_start);


	if( bf ==20 && edge_num/total_vertex <= 6){	
		for(int K=0;K<FWblockDim;K++){
				printf("K=%d phase1\n",K);
			// Phase 1
			floyd_warshall_1<<< 1,threadStr>>>( cuda_graph,K,K*bf);
			cudaCheckErrors("phase 1");

					//cudaDeviceSynchronize();

			// Phase 2
			//printf("K=%d phase2\n",K);
			
			if(FWblockDim>1){
				floyd_warshall_2<<< (FWblockDim-1)*2 ,threadStr>>>( cuda_graph,K,K*bf);
				cudaCheckErrors("phase 2 col");
				
					//cudaDeviceSynchronize();

				// Phase 3
				//printf("K=%d phase3\n",K);

			
				floyd_warshall_3<<<blockStr,threadStr>>>(cuda_graph,K,K*bf);
				
				cudaCheckErrors("phase 3");

					//cudaDeviceSynchronize();

			}
		}
	}
	else{
		for(int K=0;K<FWblockDim;K++){
			// Phase 1
			//printf("K=%d phase1\n",K);
			for(int i=0;i<bf;i++){
				floyd_warshall_beta_1<<<bf,bf>>>(cuda_graph,K,K*bf + i);
				cudaCheckErrors("phase 1");
			}


			//printf("K=%d phase2\n",K);
			//Phase 2
			

			if(FWblockDim>1){
				for(int i=0;i<bf;i++){
					floyd_warshall_beta_2<<<(FWblockDim-1)*2*bf,bf>>>(cuda_graph,K,K*bf + i );
					cudaCheckErrors("phase 2 col");
				}
				
				

				//printf("K=%d phase3\n",K);
				//Phase 3
				
				//fprintf(stderr,"qqq %d\n",(FWblockDim-1)*(FWblockDim-1)*bf);
				floyd_warshall_beta_3<<<blockStr2,bf>>>(cuda_graph,K,K*bf,FWblockDim-1);
				cudaCheckErrors("phase 3");

			}
		}
	}

	cudaDeviceSynchronize();
	// 時間計算是否要擺到前面??
	cudaEventRecord(com_stop);
	cudaEventSynchronize(com_stop);
	cudaEventElapsedTime(&com_temp,com_start,com_stop);
	com_total+=com_temp;

	cudaEventRecord(mem_start);
	cudaMemcpy(graph,cuda_graph,sizeof(int)*tempVertex*tempVertex,D2H);
	cudaCheckErrors("copy back error");
	cudaEventRecord(mem_stop);
	cudaEventSynchronize(mem_stop);
	cudaEventElapsedTime(&mem_temp,mem_start,mem_stop);
	mem_total += mem_temp;



	cudaEventRecord(io_start);
	output.open(argv[2]);
	// 每行最後面到底要不要加SPACE!!!!!!!

	for(int i=0;i<total_vertex;i++){

		for(int j=0;j<total_vertex;j++){
			if(graph[i*tempVertex+j]==INF){
				output<<"INF";
			}
			else{
				output<<graph[i*tempVertex+j];
			}
			output<<" ";
		}

	
		output<<endl;
	}
	cudaEventRecord(io_stop);
	cudaEventSynchronize(io_stop);
	cudaEventElapsedTime(&io_temp,io_start,io_stop);
	io_total += io_temp;

	cudaEventRecord(total_stop);
    cudaEventSynchronize(total_stop);
	cudaEventElapsedTime(&total_temp, total_start, total_stop);



    fprintf(stderr, "\n\n");
    fprintf(stderr, "TOTAL = %f\n", total_temp);
    fprintf(stderr, "COMPUTE = %f\n", com_total);
    fprintf(stderr, "MEMORY = %f\n", mem_total);
	fprintf(stderr, "IO = %f\n", io_total);

	return 0;
}