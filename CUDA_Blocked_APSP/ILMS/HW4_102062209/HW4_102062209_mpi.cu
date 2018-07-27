#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <mpi.h>
int* d;




__constant__ int cuda_bf;
__constant__ int cuda_total_vertex;
__constant__ int cuda_tempVertex;
__constant__ int cuda_device_num;
__constant__ int cuda_FW_block;

#define INF 1e9
#define H2D cudaMemcpyHostToDevice
#define D2H cudaMemcpyDeviceToHost
#define D2D cudaMemcpyDeviceToDevice

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

__global__ void floyd_warshall_2(int* dist,int k , int kbf  ){

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

__global__ void floyd_warshall_3(int* dist, int k ,int kbf,int ID){

	int idx,idy;

	int blockIdx_x = ((cuda_FW_block-1)/cuda_device_num)*ID + blockIdx.x;


	idy = blockIdx.y >= k? blockIdx.y + 1 : blockIdx.y;
	idx = blockIdx_x >= k? blockIdx_x + 1 : blockIdx_x;

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

__global__ void floyd_warshall_beta_3(int* dist, int k , int kbf ,int grid_size,int ID ){

	int idx,idy;

	int blockIdx_y = ((cuda_FW_block-1)/cuda_device_num)*ID + blockIdx.y;
	int temp = ((blockIdx_y*gridDim.x) + blockIdx.x) / cuda_bf;
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
	MPI_Status status;
	MPI_Request req;

	float total_temp=0,total_total=0,io_temp =0 , io_total=0 , com_temp =0,com_total=0 , mem_temp=0 , mem_total=0;

	int rc = MPI_Init(&argc,&argv);
	int rank , process_num;
	if(rc!= MPI_SUCCESS){
		printf("Error when initializing mpi \n");
  	}
  	MPI_Comm_size(MPI_COMM_WORLD,&process_num);
  	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  	/*
  	if(rank==0){
		cudaEventCreate(&total_start);
		cudaEventCreate(&total_stop);
		cudaEventCreate(&com_start);
		cudaEventCreate(&com_stop);
		cudaEventCreate(&mem_start);
		cudaEventCreate(&mem_stop);
		cudaEventCreate(&io_start);
		cudaEventCreate(&io_stop);
	}
*/

	cudaSetDevice(rank);
	cudaCheckErrors("???");

	//if(rank==0)
	//cudaEventRecord(total_start); 

//
	//struct cudaDeviceProp prop;
	//cudaGetDeviceProperties(&prop,0);
	//fprintf(stderr,"clock rate %lf\n",prop.clockRate);

	int bf = atoi(argv[3]);
	int total_vertex;
	int edge_num;

	int DEVICE_NUM = 2;
	int tempVertex;
	int * graph;// = new int[(tempVertex)*(tempVertex)];
	ifstream input;
	ofstream output;
	fprintf(stderr,"IM here\n");	
	if(rank==0){
		
		input.open(argv[1]);

		input >> total_vertex;
		input >> edge_num;

		tempVertex = total_vertex % bf ?  (total_vertex + (bf - (total_vertex%bf) )): total_vertex;
		graph = new int[(tempVertex)*(tempVertex)];

		for(int i=0;i<tempVertex;i++){
			for(int j=0;j<tempVertex;j++){
				graph[i*tempVertex+j] = INF;
			}
			graph[i*tempVertex + i ]=0;
		}

		
		//cudaEventRecord(io_start);
		//cudaCheckErrors("4");
		for(int i=0;i<edge_num;i++){
			int a,b;
			input >> a;
			input >> b;
			input >> graph[(a-1)*tempVertex + (b-1) ];
			//fprintf(stderr,"graph %d %d :%d\n",a,b,graph[a*tempVertex+b]);
		}
		MPI_Send(&total_vertex,1,MPI_INT,1,0,MPI_COMM_WORLD);
		MPI_Send(&edge_num,1,MPI_INT,1,0,MPI_COMM_WORLD);
		MPI_Send(&tempVertex,1,MPI_INT,1,0,MPI_COMM_WORLD);

		MPI_Send(graph,tempVertex*tempVertex,MPI_INT,1,0,MPI_COMM_WORLD);

	}
	else{
		MPI_Recv(&total_vertex,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
		MPI_Recv(&edge_num,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);

		MPI_Recv(&tempVertex,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
		graph = new int[(tempVertex)*(tempVertex)];
		MPI_Recv(graph,tempVertex*tempVertex,MPI_INT,0,0,MPI_COMM_WORLD,&status);

	}
	fprintf(stderr,"IM there\n");
	//fprintf(stderr,"tempVertex:%d\n",tempVertex);

	//d = new int[tempVertex*tempVertex];
	
	//cudaMallocHost((void**)&graph,sizeof(int)*tempVertex*tempVertex);

	/*
	graph = new int[(tempVertex)*(tempVertex)];

	for(int i=0;i<tempVertex;i++){
		for(int j=0;j<tempVertex;j++){
			graph[i*tempVertex+j] = INF;
		}
		graph[i*tempVertex + i ]=0;
	}

	if(rank==0)
	cudaEventRecord(io_start);
	for(int i=0;i<edge_num;i++){
		int a,b;
		input >> a;
		input >> b;
		input >> graph[(a-1)*tempVertex + (b-1) ];
		//fprintf(stderr,"graph %d %d :%d\n",a,b,graph[a*tempVertex+b]);
	}
	*/
	/*
	if(rank==0){
		cudaEventRecord(io_stop);
		cudaCheckErrors("1");
		cudaEventSynchronize(io_stop);
		cudaCheckErrors("2");
		cudaEventElapsedTime(&io_temp,io_start,io_stop);
		cudaCheckErrors("3");
		io_total += io_temp;
	}
*/
	int* cuda_graph;
	

	fprintf(stderr,"1111\n");
	cudaMalloc((void**)&cuda_graph,sizeof(int)*(tempVertex)*(tempVertex));
	cudaCheckErrors("malloc gpu");
	fprintf(stderr,"2222\n");
	

	int FWblockDim = tempVertex / bf ;

	//cudaSetDevice(0);
	if(rank==0){
		//cudaEventRecord(mem_start);	
	}
	cudaCheckErrors("oao");
	cudaSetDevice(rank);
		
	cudaMemcpy(cuda_graph,graph,sizeof(int)*tempVertex*tempVertex ,H2D);
	cudaCheckErrors("memcpy gpu");
			
	cudaMemcpyToSymbol(cuda_bf,&bf,sizeof(int));
	cudaMemcpyToSymbol(cuda_total_vertex,&total_vertex,sizeof(int));
	cudaMemcpyToSymbol(cuda_tempVertex,&tempVertex,sizeof(int));
	cudaMemcpyToSymbol(cuda_device_num,&DEVICE_NUM,sizeof(int));
	cudaMemcpyToSymbol(cuda_FW_block,&FWblockDim,sizeof(int));
		
	/*
	if(rank==0){
		cudaSetDevice(rank);
		cudaEventRecord(mem_stop);
		cudaEventSynchronize(mem_stop);
		cudaEventElapsedTime(&mem_temp,mem_start,mem_stop);
		mem_total += mem_temp;
	}
	*/
	//int FWblockDim = total_vertex%bf ? (total_vertex/bf + 1) : total_vertex/bf;
	//int remainBF = total_vertex%bf? total_vertex%bf : bf ;
	

	dim3 threadStr(bf,bf);
	dim3 blockStr((FWblockDim-1)/DEVICE_NUM,(FWblockDim-1)/DEVICE_NUM);
	dim3 blockStr_mod((FWblockDim-1)%DEVICE_NUM,(FWblockDim-1)%DEVICE_NUM);
	dim3 blockStr2((FWblockDim-1)*bf,FWblockDim-1);
	//if(rank==0)
	//cudaEventRecord(com_start);

	int* copy = new int[tempVertex*tempVertex];
	fprintf(stderr,"IM HERE\n");
	if( bf ==20 && edge_num/total_vertex <= 6){	

			//int threadId = omp_get_thread_num();
			//cudaSetDevice(threadId);
			int* type = new int[DEVICE_NUM];
		

			for(int K=0;K<FWblockDim;K++){
				
				// Phase 1

				
					int threadId = rank;
					cudaSetDevice(threadId);
					printf("K=%d phase1 id=%d\n",K,threadId);
					floyd_warshall_1<<<1,threadStr>>>(cuda_graph,K,K*bf);
					
					cudaCheckErrors("phase 1");

					//cudaDeviceSynchronize()
					// Phase 2
					printf("K=%d phase2 id=%d\n",K,threadId);
					
					if(FWblockDim>1){

							floyd_warshall_2<<< ((FWblockDim-1))*2 ,threadStr>>>(cuda_graph,K,K*bf);
							cudaCheckErrors("phase 2 col");

							// Phase 3
							if(threadId!=DEVICE_NUM-1){
								if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId+1)<=K){
									type[threadId] = 0;
								}
								else if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId+1)>K){
									type[threadId] = 1;
								}
								else{
									type[threadId] = 2;
								}
							}
							else{
								if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId) + ((FWblockDim-1)%DEVICE_NUM + (FWblockDim-1)/DEVICE_NUM)<=K){
									type[threadId] = 0;
								}
								else if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId) + ((FWblockDim-1)%DEVICE_NUM+(FWblockDim-1)/DEVICE_NUM)>K){
									type[threadId] = 1;
								}
								else{
									type[threadId] = 2;
								}
							}

									
				
							dim3 Str_normal((FWblockDim-1)/DEVICE_NUM,FWblockDim-1);
							dim3 Str_last((FWblockDim-1)/DEVICE_NUM + ((FWblockDim-1)%DEVICE_NUM), FWblockDim-1);
							printf("K=%d phase3\n",K);
							
							if(threadId==(DEVICE_NUM-1)&&(((FWblockDim-1)%DEVICE_NUM)!=0)){

								floyd_warshall_3<<<Str_last,threadStr>>>(cuda_graph,K,K*bf,threadId);
								cudaCheckErrors("phase 3 last");	
							}
							else if((FWblockDim-1)/DEVICE_NUM!=0){
								floyd_warshall_3<<<Str_normal,threadStr>>>(cuda_graph,K,K*bf,threadId);
								cudaCheckErrors("phase 3 normal");	
							}
							

					}

				

				if(FWblockDim>1){

					if(rank==0){
						cudaSetDevice(rank);
						cudaEventRecord(com_stop);
						cudaEventSynchronize(com_stop);
						cudaEventElapsedTime(&com_temp,com_start,com_stop);
						com_total += com_temp;
						cudaEventRecord(mem_start);
					}
					int offset,count;
					
					int i = rank;
								if(type[i]==2){
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*i  + tempVertex*bf;
								}
								else{
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*i ;
								}

								if(i != DEVICE_NUM-1){
									count =  tempVertex*sizeof(int)*((FWblockDim-1)/DEVICE_NUM*bf) ;
								}
								else{
									count = tempVertex*sizeof(int)*(((FWblockDim-1)/DEVICE_NUM*bf)+((FWblockDim-1)%DEVICE_NUM*bf));
								}
								if(type[i]==1){
									count += tempVertex * bf * sizeof(int);
								}

								cudaMemcpy(graph+offset,cuda_graph+offset,count,D2H);
								cudaCheckErrors("memcpy");
							//	fprintf(stderr,"ori count %d : %d\n",i,count);


					for(int j=0;j<DEVICE_NUM;j++){
						if(i==j)
							continue;
						MPI_Isend(&type[i],1,MPI_INT,j,j,MPI_COMM_WORLD,&req);	
					}
					for(int j=0;j<DEVICE_NUM;j++){
						if(i==j)
							continue;
						MPI_Recv(&type[j],1,MPI_INT,j,i,MPI_COMM_WORLD,&status);	
						//fprintf(stderr,"rank %d type%d : %d\n",rank,j,type[j]);
					}
					if(count>0){
						for(int j=0;j<DEVICE_NUM;j++){
							if(i==j)
								continue;
							MPI_Isend(&graph[offset],count/sizeof(int),MPI_INT,j,j,MPI_COMM_WORLD,&req);
						}
					}
								
					for(int j=0;j<DEVICE_NUM;j++){
								
								if(i==j)
									continue;

								fprintf(stderr,"%d %d\n",i,j);
								if(type[j]==2){
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*j  + tempVertex*bf;
								}
								else{
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*j ;
								}

								//fprintf(stderr,"OAO\n",i,j);
								if(j != DEVICE_NUM-1){
									count =  tempVertex*sizeof(int)*((FWblockDim-1)/DEVICE_NUM*bf) ;
								}
								else{
									count = tempVertex*sizeof(int)*(((FWblockDim-1)/DEVICE_NUM*bf)+((FWblockDim-1)%DEVICE_NUM*bf));
								}
								if(type[j]==1){
									count += tempVertex * bf * sizeof(int);
								}
								
								//fprintf(stderr,"i:%d j:%d offset:%d count:%d  addi%d   addoff%d  typei:%d  typej:%d \n",i,j,offset,count,cuda_graph[i],cuda_graph[i]+offset,type[i],type[j]);
								if(count>0){
									MPI_Recv(&graph[offset],count/sizeof(int),MPI_INT,j,i,MPI_COMM_WORLD,&status);
								}
								//fprintf(stderr,"i:%d j:%d offset:%d count:%d  addi%d   addoff%d  typei:%d  typej:%d \n",i,j,offset,count,cuda_graph[i],cuda_graph[i]+offset,type[i],type[j]);
								cudaMemcpy(cuda_graph+offset,graph+offset,count,H2D);
								cudaCheckErrors("memcpy");	
								
					}
					//fprintf(stderr, "QQ %d\n",rank );
						
					if(rank==0){
			           
				        cudaEventRecord(mem_stop);
				        cudaEventSynchronize(mem_stop);
				        cudaEventElapsedTime(&mem_temp, mem_start, mem_stop);
						mem_total += mem_temp;
						cudaCheckErrors("mem end");
					}
				}
				if(rank==0){
				cudaEventRecord(com_start);
				cudaCheckErrors("com start");	
				}
			}
	}
	else{
			int* type = new int[DEVICE_NUM];
		
			
			for(int K=0;K<FWblockDim;K++){
				// Phase 1
				
					int threadId = rank;
					cudaSetDevice(threadId);

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
						
							if(threadId!=DEVICE_NUM-1){
									if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId+1)<=K){
										type[threadId] = 0;
									}
									else if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId+1)>K){
										type[threadId] = 1;
									}
									else{
										type[threadId] = 2;
									}
								}
								else{
									if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId) + ((FWblockDim-1)%DEVICE_NUM + (FWblockDim-1)/DEVICE_NUM)<=K){
										type[threadId] = 0;
									}
									else if(((FWblockDim-1)/DEVICE_NUM)*threadId<K&&((FWblockDim-1)/DEVICE_NUM)*(threadId) + ((FWblockDim-1)%DEVICE_NUM+(FWblockDim-1)/DEVICE_NUM)>K){
										type[threadId] = 1;
									}
									else{
										type[threadId] = 2;
									}
								}


						//printf("K=%d phase3\n",K);
						//Phase 3
						dim3 Str_normal((FWblockDim-1)*bf,(FWblockDim-1)/DEVICE_NUM);
						dim3 Str_last((FWblockDim-1)*bf, (FWblockDim-1)/DEVICE_NUM + ((FWblockDim-1)%DEVICE_NUM) );
						if(threadId==(DEVICE_NUM-1)&&(((FWblockDim-1)%DEVICE_NUM)!=0)){
							floyd_warshall_beta_3<<<Str_last,bf>>>(cuda_graph,K,K*bf,FWblockDim-1,threadId);
						}
						else if((FWblockDim-1)/DEVICE_NUM!=0){
							floyd_warshall_beta_3<<<Str_normal,bf>>>(cuda_graph,K,K*bf,FWblockDim-1,threadId);
						}
						
						cudaCheckErrors("phase 3");

					}
				
				
				if(FWblockDim>1){

					if(rank==0){
						/*
						cudaSetDevice(rank);
						cudaEventRecord(com_stop);
						cudaEventSynchronize(com_stop);
						cudaEventElapsedTime(&com_temp,com_start,com_stop);
						com_total += com_temp;
						cudaEventRecord(mem_start);
						*/
					}
					int offset,count;
					
					int i = rank;
								if(type[i]==2){
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*i  + tempVertex*bf;
								}
								else{
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*i ;
								}

								if(i != DEVICE_NUM-1){
									count =  tempVertex*sizeof(int)*((FWblockDim-1)/DEVICE_NUM*bf) ;
								}
								else{
									count = tempVertex*sizeof(int)*(((FWblockDim-1)/DEVICE_NUM*bf)+((FWblockDim-1)%DEVICE_NUM*bf));
								}
								if(type[i]==1){
									count += tempVertex * bf * sizeof(int);
								}

								cudaMemcpy(graph+offset,cuda_graph+offset,count,D2H);
								cudaCheckErrors("memcpy");
							//	fprintf(stderr,"ori count %d : %d\n",i,count);


					for(int j=0;j<DEVICE_NUM;j++){
						if(i==j)
							continue;
						MPI_Isend(&type[i],1,MPI_INT,j,j,MPI_COMM_WORLD,&req);	
					}
					for(int j=0;j<DEVICE_NUM;j++){
						if(i==j)
							continue;
						MPI_Recv(&type[j],1,MPI_INT,j,i,MPI_COMM_WORLD,&status);	
						//fprintf(stderr,"rank %d type%d : %d\n",rank,j,type[j]);
					}
					if(count>0){
						for(int j=0;j<DEVICE_NUM;j++){
							if(i==j)
								continue;
							MPI_Isend(&graph[offset],count/sizeof(int),MPI_INT,j,j,MPI_COMM_WORLD,&req);
						}
					}
								
					for(int j=0;j<DEVICE_NUM;j++){
								
								if(i==j)
									continue;

								//fprintf(stderr,"%d %d\n",i,j);
								if(type[j]==2){
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*j  + tempVertex*bf;
								}
								else{
									offset = tempVertex*((FWblockDim-1)/DEVICE_NUM*bf )*j ;
								}

								//fprintf(stderr,"OAO\n",i,j);
								if(j != DEVICE_NUM-1){
									count =  tempVertex*sizeof(int)*((FWblockDim-1)/DEVICE_NUM*bf) ;
								}
								else{
									count = tempVertex*sizeof(int)*(((FWblockDim-1)/DEVICE_NUM*bf)+((FWblockDim-1)%DEVICE_NUM*bf));
								}
								if(type[j]==1){
									count += tempVertex * bf * sizeof(int);
								}
								
								//fprintf(stderr,"i:%d j:%d offset:%d count:%d  addi%d   addoff%d  typei:%d  typej:%d \n",i,j,offset,count,cuda_graph[i],cuda_graph[i]+offset,type[i],type[j]);
								if(count>0){
									MPI_Recv(&graph[offset],count/sizeof(int),MPI_INT,j,i,MPI_COMM_WORLD,&status);
								}
								//fprintf(stderr,"i:%d j:%d offset:%d count:%d  addi%d   addoff%d  typei:%d  typej:%d \n",i,j,offset,count,cuda_graph[i],cuda_graph[i]+offset,type[i],type[j]);
								cudaMemcpy(cuda_graph+offset,graph+offset,count,H2D);
								cudaCheckErrors("memcpy");	
								
					}
					//fprintf(stderr, "QQ %d\n",rank );
						
					if(rank==0){
			           /*
				        cudaEventRecord(mem_stop);
				        cudaEventSynchronize(mem_stop);
				        cudaEventElapsedTime(&mem_temp, mem_start, mem_stop);
						mem_total += mem_temp;
						cudaCheckErrors("mem end");
						*/
					}
				}
				if(rank==0){
					/*
					cudaEventRecord(com_start);
					cudaCheckErrors("com start");	
					*/
				}
		}
	}

	MPI_Finalize();


	fprintf(stderr,"IM THERE %d\n",rank);
	//fprintf(stderr,"%d QAQ \n",rank);
//	cudaSetDevice(rank);
	//cudaDeviceSynchronize();

	// 時間計算是否要擺到前面??
	//if(rank==0){
	//cudaEventRecord(com_stop);
	//cudaEventSynchronize(com_stop);
	//cudaEventElapsedTime(&com_temp,com_start,com_stop);
	//com_total+=com_temp;
//	}
	
	//fprintf(stderr,"%d QAQ 2\n",rank);
	//fprintf(stderr,"qqq %s   %s  \n",typeid((&graph[0])+10),typeid(cuda_graph[1]+10));
	if(rank==0){
		//cudaEventRecord(mem_start);
		cudaMemcpy(graph,cuda_graph,sizeof(int)*tempVertex*tempVertex,D2H);

		//fprintf(stderr,"%d QAQ QAQ\n",rank);
		//cudaCheckErrors("copy back error");
		//cudaEventRecord(mem_stop);
		//cudaEventSynchronize(mem_stop);
		//cudaEventElapsedTime(&mem_temp,mem_start,mem_stop);
		mem_total += mem_temp;
	}
	cudaCheckErrors("QQQ");
		//fprintf(stderr,"%d QAQ3 \n",rank);
	/*
	#pragma omp parallel num_threads(DEVICE_NUM)
	{

		int* tempGraph = new int[tempVertex*tempVertex];
		int threadId = omp_get_thread_num();

		cudaSetDevice(threadId);
		
		int offset = tempVertex*(tempVertex/DEVICE_NUM) * threadId;
		int count = threadId==(DEVICE_NUM-1)? tempVertex*( (tempVertex/DEVICE_NUM) + (tempVertex%DEVICE_NUM)) : tempVertex*(tempVertex/DEVICE_NUM) ;
			


		cudaMemcpy(tempGraph,(cuda_graph[threadId]),sizeof(int)*tempVertex*tempVertex,D2H);
		for(int i=offset;i<offset+count;i++){
			graph[i] = tempGraph[i];
		}
		cudaCheckErrors("copy back error");
	

	}
	*/

	if(rank==0){
		//cudaEventRecord(io_start);

	
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
		
	//
		//cudaEventRecord(io_stop);
		//cudaEventSynchronize(io_stop);
		//cudaEventElapsedTime(&io_temp,io_start,io_stop);
		//io_total += io_temp;
	
	///cudaEventRecord(total_stop);
   // cudaEventSynchronize(total_stop);
	//cudaEventElapsedTime(&total_temp, total_start, total_stop);
	}

	if(rank==0){
	    fprintf(stderr, "\n\n");
	    fprintf(stderr, "TOTAL = %f\n", total_temp);
	    fprintf(stderr, "COMPUTE = %f\n", com_total);
	    fprintf(stderr, "MEMORY = %f\n", mem_total);
		fprintf(stderr, "IO = %f\n", io_total);
	}
	return 0;
}