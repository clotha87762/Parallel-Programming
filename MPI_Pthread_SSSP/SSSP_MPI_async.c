#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int** weight;


int main(int argc,char* argv[]){

	int m_dist;
	int father;
	

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
	int reactive =0;
	int* startArray;
	int* endArray;
	int* valueArray;
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
	


	//printf("rank%d process_num%d\n",rank,process_num);

	if(rank==0){


		//ifstream input;
		FILE* fi;
		//printf("rank%d 01\n",rank);
		//input.open(argv[2]);
		//input >> total_vertex;
		//input >> edge_num;
		fi = fopen(argv[2],"r");
		fscanf(fi,"%d %d",&total_vertex,&edge_num);

		startArray = (int*)malloc(sizeof(int)*edge_num);
		endArray = (int*)malloc(sizeof(int)*edge_num);
		valueArray = (int*)malloc(sizeof(int)*edge_num);
	
		reqs = (MPI_Request*)malloc(sizeof(MPI_Request)*(total_vertex+3));
		
		weight = (int**) malloc (sizeof(int*) * total_vertex);

		for(i=0;i<total_vertex;i++)
			weight[i] = (int*) malloc(sizeof(int)*total_vertex);

		//printf("total_vertex%d  edge%d\n",total_vertex,edge_num);
		for(i=0;i<total_vertex;i++)
			for(j=0;j<total_vertex;j++)
				weight[i][j] = -1;

		for(i=0;i<total_vertex;i++){
			weight[i][i] = 0;
		}
		int from,to;
	//	printf("www \n");
		for(i=0;i<edge_num;i++){
			//input>>from;
			//input>>to;
			//input>>weight[from-1][to-1];
			fscanf(fi,"%d %d ",&from,&to);
			fscanf(fi,"%d ",&valueArray[i]);
			//weight[to-1][from-1] = weight[from-1][to-1];
			startArray[i] = from-1;
			endArray[i] = to -1;
			//valueArray[i] = weight[to-1][from-1];
			//printf("weight[%d][%d]: %d\n",from-1,to-1,weight[from-1][to-1]);
		}

		//printf("rank%d 02\n",rank);
	}
	
//	printf("000 %d\n",rank);
	

	MPI_Bcast(&total_vertex,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&edge_num,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank!=0){
		startArray = (int*)malloc(sizeof(int)*edge_num);
		endArray = (int*)malloc(sizeof(int)*edge_num);
		valueArray = (int*)malloc(sizeof(int)*edge_num);
		reqs = (MPI_Request*)malloc(sizeof(MPI_Request)*(total_vertex+3));
		
		weight = (int**) malloc (sizeof(int*) * total_vertex);

		for(i=0;i<total_vertex;i++){
			weight[i] = (int*) malloc(sizeof(int)*total_vertex);
			for(j=0;j<total_vertex;j++)
				weight[i][j] = -1;
		}
	}
//	printf("111 %d\n",rank);


	for(i=0;i<total_vertex;i++)
		weight[i][i] = 0;
	
	MPI_Bcast(&startArray[0],edge_num,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&endArray[0],edge_num,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&valueArray[0],edge_num,MPI_INT,0,MPI_COMM_WORLD);
	for(i=0;i<edge_num;i++){
		weight[startArray[i]][endArray[i]] = valueArray[i];
		weight[endArray[i]][startArray[i]] = valueArray[i];
	}

	//for(i=0;i<total_vertex;i++)
	//MPI_Bcast(weight[i],total_vertex,MPI_INT,0,MPI_COMM_WORLD);
	
//	if(rank==1){
//		for(i=0;i<total_vertex;i++)
//			for(j=0;j<total_vertex;j++)
//				printf("[%d][%d]:%d\n",i,j,weight[i][j]);
//	}

	if(rank==source_id){
		m_dist = 0;
	}
	else{
		m_dist = 210000000;
	}
	father = rank;

	
	allFlags = (int*) malloc(sizeof(int)*process_num);


	for(i=0;i<total_vertex;i++){
		if(weight[rank][i]>=0&&i!=rank)
			outCount++;
		if(weight[i][rank]>=0&&i!=rank&&rank!=source_id)
			inCount++;
	}
	//printf("222 %d\n",rank);

	outList =(int*) malloc(sizeof(int)*outCount);
	outWeight = (int*) malloc(sizeof(int)*outCount);
	recvBuffer = (int*) malloc(sizeof(int)*3);
	outBuffer = (int*) malloc(sizeof(int)*3);
	outCount = 0;


	for(i=0;i<total_vertex;i++){
		if(weight[rank][i]>=0&&i!=rank){
			outList[outCount] = i;
			outWeight[outCount] = weight[rank][i];
			//outBuffer[outCount*2] = 210000000;
			//outBuffer[outCount*2+1] = rank;
			outCount++;
		}
	}
//	printf("out %d  = %d\n",rank,outCount);
//	printf("in %d  = %d\n",rank,inCount);

	//if(rank==source_id){
	//	for(i=0;i<outCount;i++){
	//		outBuffer[i*2 +0] = m_dist+outWeight[i];
	//	}
	//}

	//printf("aaa %d\n",rank);

	
	for(i=0;i<outCount;i++){
		//printf("!! %d\n",rank);
		if(outList[i]<rank){
			reactive = 1;
			break;
		}
	}
	

//	printf("@@ %d\n",rank);
	// BUFFER : [0] = value  [1] = source [2] = what(0 = dist,1=token,2=terminate)   white = 0 , black = 1


	// 注意！！！！！！！！！   之後要回來完成
	//  process 數 = 1 得版本！！

	outBuffer[1] = rank;
	if(rank==source_id){
			//printf("kkk %d\n",rank);
		for(i=0;i<outCount;i++){
			outBuffer[0] = m_dist+outWeight[i];
			outBuffer[1] = rank;
			outBuffer[2] = 0;
			MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
		}
		if(rank==0){
			outBuffer[0] = 0;
			outBuffer[2] = 1;
			MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
		}
	}

	int black = 0;
	int tokenFlag = 0;
	int firstTime = 1;
	while(1){
		
		MPI_Recv(&recvBuffer[0],3,MPI_INT,MPI_ANY_SOURCE,rank,MPI_COMM_WORLD,&status);
		if(rank!=0){

			if(recvBuffer[2]==0){ // DIST
				if(recvBuffer[0]<m_dist){
					
					m_dist = recvBuffer[0];
					father = recvBuffer[1];

					for(i=0;i<outCount;i++){
						outBuffer[0] = m_dist+outWeight[i];
						outBuffer[1] = rank;
						outBuffer[2] = 0;
						MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
					}
					if(reactive==1)
						black = 1;
					if(tokenFlag==1){
						tokenFlag = 0;
						outBuffer[0] = black;
						outBuffer[2] = 1;
						MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
						black = 0;
					}

				}
			}
			else if(recvBuffer[2]==1){ // TOKEN
				if(m_dist<210000000){
					
					if(black==1){
						outBuffer[0] = 1;
					}
					else{
						outBuffer[0] = recvBuffer[0];
					}
					outBuffer[2] = 1;
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
					black = 0;
				
				}
				else{
					tokenFlag = 1;
				}
			}
			else if(recvBuffer[2]==2){ // TERMINATE
				break;
			}

		}
		else{

			if(recvBuffer[2]==0){

				if(recvBuffer[0]<m_dist){
					
					m_dist = recvBuffer[0];
					father = recvBuffer[1];
					for(i=0;i<outCount;i++){
						outBuffer[0] = m_dist+outWeight[i];
						outBuffer[1] = rank;
						outBuffer[2] = 0;
						MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
					}
				}


				if(firstTime&&rank!=source_id){
					firstTime = 0;
					outBuffer[0] = 0;
					outBuffer[2] = 1;
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
				}
			}
			else if(recvBuffer[2]==1){
				
				if(recvBuffer[0]==0){
					
					outBuffer[0] = 0;
					outBuffer[1] = 0;
					outBuffer[2] = 2;
					for(i=1;i<process_num;i++){
						MPI_Isend(&outBuffer[0],3,MPI_INT,i,i,MPI_COMM_WORLD,&reqs[outCount+i]);
					}
					break;
				}
				else{
					outBuffer[0] = 0;
					outBuffer[2] = 1;
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
				}
			}


		}

	}
	
	//printf("bbb %d\n",rank);
	//printf("rank%d father:%d\n",rank,father);
	int* father_array = (int*)malloc(sizeof(int)*(total_vertex+1)); 

	//MPI_Gather(&father,1,MPI_INT,&father_array[0],process_num,MPI_INT,0,MPI_COMM_WORLD);

	if(rank==0){

			//int* father_array = (int*) malloc(sizeof(int)*process_num); 
			father_array[0] = father;
			for(i=1;i<process_num;i++){
				//printf("%d\n",i);
				MPI_Recv(&father_array[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
			}

			FILE* fo = fopen(argv[3],"w");

			//ofstream output;
			//output.open(argv[3]);
			int* tempArray = (int*) malloc(sizeof(int)*total_vertex);
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
				for(j=count-1;j>=0;j--){
					fprintf(fo,"%d ",tempArray[j]);
					//output<<tempArray[j]<<" ";
				}
				fprintf(fo,"\n");
			}

			fclose(fo);

	}
	else{
		MPI_Send(&father,1,MPI_INT,0,rank,MPI_COMM_WORLD);
	}
	
	
	//printf("ccc %d\n",rank);

	MPI_Finalize();

	return 0;
}