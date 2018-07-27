#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>



using namespace std;




int main(int argc,char* argv[]){

	int m_dist;
int father;
int** weight;

int source_id;
	int rank,rc,process_num;
	int edge_num,total_vertex;
	int* outList;
	int* outWeight;
	int* recvBuffer;
	int* outBuffer;
	int outCount = 0;
	int inCount = 0;
	int nowDist;
	int updateFlag = 0;
	int* allFlags;
	MPI_Request* reqs;
	MPI_Status status;

	rc = MPI_Init(&argc,&argv);
	int i,j;
	 if(rc!= MPI_SUCCESS){
	printf("Error when initializing mpi \n");
  }



  	MPI_Comm_size(MPI_COMM_WORLD,&process_num);
  	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	source_id = atoi(argv[4]) - 1;

	weight = new int*[process_num];

	for(i=0;i<total_vertex;i++)
		weight[i] = new int[process_num];



	


	printf("rank%d process_num%d\n",rank,process_num);

	if(rank==0){


		ifstream input;
		printf("rank%d 01\n",rank);
		input.open(argv[2]);
		input >> total_vertex;
		input >> edge_num;

		
		printf("total_vertex%d  edge%d\n",total_vertex,edge_num);
		for(i=0;i<total_vertex;i++)
			for(j=0;j<total_vertex;j++)
				weight[i][j] = -1;

		for(i=0;i<total_vertex;i++){
			weight[i][i] = 0;
		}
		int from,to;
		printf("www \n");
		for(i=0;i<edge_num;i++){
			input>>from;
			input>>to;
			input>>weight[from-1][to-1];
			weight[to-1][from-1] = weight[from-1][to-1];
		}

		printf("rank%d 02\n",rank);
	}
	
	printf("000 %d\n",rank);
	

	MPI_Bcast(&total_vertex,1,MPI_INT,0,MPI_COMM_WORLD);
	printf("111 %d\n",rank);
	MPI_Bcast(&edge_num,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&weight[0][0],total_vertex*total_vertex,MPI_INT,0,MPI_COMM_WORLD);
	

	if(rank==source_id){
		m_dist = 0;
	}
	else{
		m_dist = 210000000;
	}
	father = rank;

	
	allFlags = new int[process_num];


	for(i=0;i<total_vertex;i++){
		if(weight[rank][i]>=0&&i!=rank)
			outCount++;
		if(weight[i][rank]>=0&&i!=rank&&rank!=source_id)
			inCount++;
	}
	printf("222 %d\n",rank);

	outList = new int[outCount];
	outWeight = new int[outCount];
	recvBuffer = new int[inCount*2];
	reqs = new MPI_Request[outCount];
	outBuffer = new int[outCount*2];
	outCount = 0;


	for(i=0;i<total_vertex;i++){
		if(weight[rank][i]>=0&&i!=rank){
			outList[outCount] = i;
			outWeight[outCount] = weight[rank][i];
			outBuffer[outCount*2] = 210000000;
			outBuffer[outCount*2+1] = rank;
			outCount++;
		}
	}
	printf("out %d  = %d\n",rank,outCount);
	printf("in %d  = %d\n",rank,inCount);
	bool termination = false;

	if(rank==source_id){
		for(i=0;i<outCount;i++){
			outBuffer[i*2 +0] = m_dist+outWeight[i];
		}
	}

	printf("aaa %d\n",rank);

	

	while(1){
		
		termination = true;
		updateFlag = 0;

		for(i=0;i<outCount;i++){
			printf("QQ %d\n",rank);
			MPI_Isend(&outBuffer[i*2],2,MPI_INT,outList[i],0,MPI_COMM_WORLD,&reqs[i]);
		}
		for(i=0;i<inCount;i++){
			MPI_Recv(&recvBuffer[i*2],2,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
			if(recvBuffer[i*2]<m_dist){
				m_dist = recvBuffer[i*2];
				father = recvBuffer[i*2+1];
				updateFlag = 1;
			}
		}


		// 不知道這邊要不要@@
		for(i=0;i<outCount;i++){
			MPI_Wait(&reqs[i],&status);
		}

	
		for(i=0;i<outCount;i++){
			outBuffer[i*2 +0] = m_dist+outWeight[i];
		}

		MPI_Allgather(&updateFlag,1,MPI_INT,allFlags,1,MPI_INT,MPI_COMM_WORLD);
		for(i=0;i<process_num;i++){
			if(allFlags[i]==1){
				termination = false;
				break;
			}
		}
		if(termination)
			break;

	}
	
	printf("bbb %d\n",rank);

	if(rank==0){

		int* father_array = new int[process_num]; 
		for(i=0;i<process_num;i++){
			MPI_Recv(&father_array[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
		}

			ofstream output;
			output.open(argv[3]);
			int* tempArray = new int[total_vertex];
			int x,count;
			for(i=1;i<total_vertex+1;i++){
				x = i-1;
				tempArray[0] = i;
				count = 1;
				while(father_array[x]!=source_id){
					tempArray[count++] = father_array[x]+1;
					x = father_array[x] ;
				}
				tempArray[count++]=source_id + 1;
				for(int j=count-1;j>=0;j--){
					output<<tempArray[j]<<" ";
				}
				output<<endl;
			}

	}
	else{
		MPI_Send(&father,1,MPI_INT,0,0,MPI_COMM_WORLD);
	}

	
	printf("ccc %d\n",rank);

	MPI_Finalize();

	return 0;
}