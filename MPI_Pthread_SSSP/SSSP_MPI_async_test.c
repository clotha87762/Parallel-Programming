#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


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
	MPI_Request* reqs;
	MPI_Status status;
	double t1,t2;
	double startTime;
	double endTime;
	double ioTimer = 0;
	double commuTimer = 0;
	double cpu1,cpu2;
	double cpuTime;
	



	rc = MPI_Init(&argc,&argv);
	int i,j;
	 if(rc!= MPI_SUCCESS){
	printf("Error when initializing mpi \n");
  }

  	startTime = MPI_Wtime();

  	MPI_Comm_size(MPI_COMM_WORLD,&process_num);
  	MPI_Comm_rank(MPI_COMM_WORLD,&rank);



	source_id = atoi(argv[4]) - 1;
	


//	printf("rank%d process_num%d\n",rank,process_num);

	//if(rank==0){


		//ifstream input;
		FILE* fi;
		//printf("rank%d 01\n",rank);
		//input.open(argv[2]);
		//input >> total_vertex;
		//input >> edge_num;
		t1 = MPI_Wtime();
		fi = fopen(argv[2],"r");
		fscanf(fi,"%d %d",&total_vertex,&edge_num);
		t2 = MPI_Wtime();
		ioTimer += (t2-t1);

		reqs = (MPI_Request*)malloc(sizeof(MPI_Request)*(total_vertex+3));

		weight = (int**) malloc (sizeof(int*) * total_vertex);

		for(int i=0;i<total_vertex;i++)
			weight[i] = (int*) malloc(sizeof(int)*total_vertex);


		//printf("total_vertex%d  edge%d\n",total_vertex,edge_num);
		for(i=0;i<total_vertex;i++)
			for(j=0;j<total_vertex;j++)
				weight[i][j] = -1;

		for(i=0;i<total_vertex;i++){
			weight[i][i] = 0;
		}
		int from,to;
		//printf("www \n");

		t1 = MPI_Wtime();
		for(i=0;i<edge_num;i++){
			//input>>from;
			//input>>to;
			//input>>weight[from-1][to-1];
			fscanf(fi,"%d %d ",&from,&to);
			fscanf(fi,"%d ",&weight[from-1][to-1]);
			weight[to-1][from-1] = weight[from-1][to-1];
			//printf("weight[%d][%d]: %d\n",from-1,to-1,weight[from-1][to-1]);
		}
		t2 = MPI_Wtime();
		ioTimer += (t2-t1);

		//printf("rank%d 02\n",rank);
//	}
	
//	printf("000 %d\n",rank);
	
	//MPI_Bcast(&total_vertex,1,MPI_INT,0,MPI_COMM_WORLD);
	//t1 = MPI_Wtime();
	//MPI_Bcast(&edge_num,1,MPI_INT,0,MPI_COMM_WORLD);
	//for(i=0;i<total_vertex;i++)
	//MPI_Bcast(&weight[i],total_vertex,MPI_INT,0,MPI_COMM_WORLD);
	//t2 = MPI_Wtime();
	//commuTimer += (t2-t1);

	double bcastTime = commuTimer;

	cpu1 = MPI_Wtime();

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
	//printf("out %d  = %d\n",rank,outCount);
	//printf("in %d  = %d\n",rank,inCount);

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
	

	//printf("@@ %d\n",rank);
	// BUFFER : [0] = value  [1] = source [2] = what(0 = dist,1=token,2=terminate)   white = 0 , black = 1


	// 注意！！！！！！！！！   之後要回來完成
	//  process 數 = 1 得版本！！
	int sendCount = 0;
	int recvCount = 0;

	outBuffer[1] = rank;
	if(rank==source_id){
		//	printf("kkk %d\n",rank);
		for(i=0;i<outCount;i++){
			outBuffer[0] = m_dist+outWeight[i];
			outBuffer[1] = rank;
			outBuffer[2] = 0;
			MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
			sendCount+=outCount;
		}
		if(rank==0){
			outBuffer[0] = 0;
			outBuffer[2] = 1;
			MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
			sendCount++;
		}
	}

	int black = 0;
	int tokenFlag = 0;
	int firstTime = 1;

	while(1){
		
		t1 = MPI_Wtime();
		MPI_Recv(&recvBuffer[0],3,MPI_INT,MPI_ANY_SOURCE,rank,MPI_COMM_WORLD,&status);
		t2 = MPI_Wtime();
		commuTimer += (t2-t1);

		recvCount++;

		if(rank!=0){

			if(recvBuffer[2]==0){ // DIST
				if(recvBuffer[0]<m_dist){
					
					m_dist = recvBuffer[0];
					father = recvBuffer[1];

					for(i=0;i<outCount;i++){
						outBuffer[0] = m_dist+outWeight[i];
						outBuffer[1] = rank;
						outBuffer[2] = 0;
						
						t1 = MPI_Wtime();
						MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
						t2 = MPI_Wtime();
						commuTimer += (t2-t1);

						
					}
					sendCount+= outCount;
					if(reactive==1)
						black = 1;
					if(tokenFlag==1){
						tokenFlag = 0;
						outBuffer[0] = black;
						outBuffer[2] = 1;

						t1 = MPI_Wtime();
						MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
						t2 = MPI_Wtime();
						commuTimer += (t2-t1);
						black = 0;

						sendCount++;
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

					t1 = MPI_Wtime();
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
					t2 = MPI_Wtime();
					commuTimer += (t2-t1);
					black = 0;

					sendCount++;
				
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
						
						t1 = MPI_Wtime();
						MPI_Isend(&outBuffer[0],3,MPI_INT,outList[i],outList[i],MPI_COMM_WORLD,&reqs[i]);
						t2 = MPI_Wtime();
						commuTimer += (t2-t1);

					}
					sendCount+=outCount;
				}


				if(firstTime&&rank!=source_id){
					firstTime = 0;
					outBuffer[0] = 0;
					outBuffer[2] = 1;
					
					t1 = MPI_Wtime();
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
					t2 = MPI_Wtime();
					commuTimer += (t2-t1);
					sendCount++;
				}

			}
			else if(recvBuffer[2]==1){
				
				if(recvBuffer[0]==0){
					
					outBuffer[0] = 0;
					outBuffer[1] = 0;
					outBuffer[2] = 2;
					for(i=1;i<process_num;i++){
						
						t1 = MPI_Wtime();
						MPI_Isend(&outBuffer[0],3,MPI_INT,i,i,MPI_COMM_WORLD,&reqs[outCount+i]);
						t2 = MPI_Wtime();
						commuTimer += (t2-t1);
					}
					sendCount+=process_num;
					break;
				}
				else{
					outBuffer[0] = 0;
					outBuffer[2] = 1;
					
					t1 = MPI_Wtime();
					MPI_Isend(&outBuffer[0],3,MPI_INT,(rank+1)%process_num,(rank+1)%process_num,MPI_COMM_WORLD,&reqs[outCount]);
					t2 = MPI_Wtime();
					commuTimer += (t2-t1);
					sendCount++;
				}
			}


		}

	}
	cpu2 = MPI_Wtime();
	cpuTime = (cpu2-cpu1) - (commuTimer) + bcastTime ;
	//printf("bbb %d\n",rank);
	//printf("rank%d father:%d\n",rank,father);
	int* father_array = (int*)malloc(sizeof(int)*(total_vertex+1)); 
	//double* syncArray = (double*)malloc(sizeof(double)*total_vertex); 
	double* commuArray = (double*) malloc(sizeof(double)*total_vertex);
	double* cpuArray = (double*) malloc(sizeof(double)*total_vertex);
	double qq[3];
	if(rank==0){

			//int* father_array = (int*) malloc(sizeof(int)*process_num); 
			father_array[0] = father;

			t1 = MPI_Wtime();
			for(i=1;i<process_num;i++){
				//printf("%d\n",i);
				MPI_Recv(&qq[0],3,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
				father_array[i] = (double)qq[0];
				//syncArray[i] = qq[1];
				commuArray[i] = qq[1];
				cpuArray[i] = qq[2];
			}
			t2 = MPI_Wtime();
			commuTimer += (t2-t1);

			recvCount += (process_num - 1);


			t1 = MPI_Wtime();
			FILE* fo = fopen(argv[3],"w");
			t2 = MPI_Wtime();
			ioTimer += (t2-t1);

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

				t1 = MPI_Wtime();
				for(j=count-1;j>=0;j--){
					fprintf(fo,"%d ",tempArray[j]);
					//output<<tempArray[j]<<" ";
				}
				fprintf(fo,"\n");
				t2 = MPI_Wtime();
				ioTimer += (t2-t1);

			}
			t1 = MPI_Wtime();
			fclose(fo);
			t2 = MPI_Wtime();
			ioTimer += (t2-t1);

	}
	else{
		qq[0] = (double) father;
		//qq[1] = syncTimer;
		qq[1] = commuTimer;
		qq[2] = cpuTime;
		t1 = MPI_Wtime();
		MPI_Send(&qq[0],3,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		commuTimer += (t2-t1);

		sendCount++;
	}
	
	//printf("ccc %d\n",rank);

	endTime = MPI_Wtime();
	double totalTime = endTime - startTime;
	
	if(rank==0){
	double CPU , SYNC, COMMU;
	CPU = 0;
	SYNC = 0;
	COMMU = 0;
	for(i=0;i<process_num;i++){
		CPU += cpuArray[i];
		//SYNC += syncArray[i];
		COMMU += commuArray[i];
	}
	CPU = CPU / (double) process_num;
	//SYNC = SYNC / (double) process_num;
	COMMU = COMMU / (double) process_num;

	
		fprintf(stderr,"rank:%d |send count:%d | recv count: %d |total Time:%lf |compute Time:%lf |commuTime:%lf |Sync Time:%lf |IO Time:%lf\n",rank,sendCount,recvCount,totalTime,CPU ,COMMU,SYNC,ioTimer);
	}
	else{
			fprintf(stderr,"rank:%d |send count:%d | recv count: %d |total Time:%lf |compute Time:%lf |commuTime:%lf|IO Time:%lf\n",rank,sendCount,recvCount,totalTime,cpuTime ,commuTimer,ioTimer);

	}

	MPI_Finalize();

	return 0;
}