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
	double t1,t2;
	double startTime;
	double endTime;
	double ioTimer = 0;
	double commuTimer = 0;
	double syncTimer = 0;
	double cpu1,cpu2;
	double cpuTime;
	MPI_Request* reqs;
	MPI_Status status;

	rc = MPI_Init(&argc,&argv);
	int i,j;
	 if(rc!= MPI_SUCCESS){
	printf("Error when initializing mpi \n");
  }

  	startTime = MPI_Wtime();

  	MPI_Comm_size(MPI_COMM_WORLD,&process_num);
  	MPI_Comm_rank(MPI_COMM_WORLD,&rank);



	source_id = atoi(argv[4]) - 1;
	


	//if(rank==0){


		//ifstream input;
		FILE* fi;
		
		//input.open(argv[2]);
		//input >> total_vertex;
		//input >> edge_num;
		t1 = MPI_Wtime();
		fi = fopen(argv[2],"r");
		fscanf(fi,"%d %d",&total_vertex,&edge_num);
		t2 = MPI_Wtime();
		ioTimer += (t2-t1);


		reqs = (MPI_Request*)malloc(sizeof(MPI_Request)*total_vertex);

		weight = (int**) malloc (sizeof(int*) * total_vertex);

		for(int i=0;i<total_vertex;i++)
			weight[i] = (int*) malloc(sizeof(int)*total_vertex);

		
		for(i=0;i<total_vertex;i++)
			for(j=0;j<total_vertex;j++)
				weight[i][j] = -1;

		for(i=0;i<total_vertex;i++){
			weight[i][i] = 0;
		}
		int from,to;

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

		
	//}
	
	
	//MPI_Bcast(&total_vertex,1,MPI_INT,0,MPI_COMM_WORLD);
	//for(i=0;i<total_vertex;i++)
	//MPI_Bcast(&weight[i],total_vertex,MPI_INT,0,MPI_COMM_WORLD);
	
	//t1 = MPI_Wtime();
	//MPI_Bcast(&edge_num,1,MPI_INT,0,MPI_COMM_WORLD);
	//MPI_Bcast(&total_vertex,1,MPI_INT,0,MPI_COMM_WORLD);
	//t2 = MPI_Wtime();
	//commuTimer += (t2-t1);

	double bcastTime = commuTimer;

	cpu1 = MPI_Wtime();

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

	outList =(int*) malloc(sizeof(int)*outCount);
	outWeight = (int*) malloc(sizeof(int)*outCount);
	recvBuffer = (int*) malloc(sizeof(int)*inCount*2);
	outBuffer = (int*) malloc(sizeof(int)*outCount*2);
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

	int termination = 0;

	if(rank==source_id){
		for(i=0;i<outCount;i++){
			outBuffer[i*2 +0] = m_dist+outWeight[i];
		}
	}


	int sendCount = 0;
	int recvCount = 0;



	while(1){
		
		termination = 1;
		updateFlag = 0;

		t1 = MPI_Wtime();
		for(i=0;i<outCount;i++){
			MPI_Isend(&outBuffer[i*2],2,MPI_INT,outList[i],0,MPI_COMM_WORLD,&reqs[i]);
		}
		t2 = MPI_Wtime();
		commuTimer += (t2-t1);

		for(i=0;i<inCount;i++){
			
			t1 = MPI_Wtime();
			MPI_Recv(&recvBuffer[i*2],2,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
			t2 = MPI_Wtime();
			commuTimer += (t2-t1);

			if(recvBuffer[i*2]<m_dist){
				m_dist = recvBuffer[i*2];
				father = recvBuffer[i*2+1];
				updateFlag = 1;
			}
		}

		sendCount+=outCount;
		recvCount+=inCount;


		// 不知道這邊要不要@@
		//for(i=0;i<outCount;i++){
		//	MPI_Wait(&reqs[i],&status);
		//}

	
		for(i=0;i<outCount;i++){
			outBuffer[i*2 +0] = m_dist+outWeight[i];
		}


		t1 = MPI_Wtime();
		MPI_Allgather(&updateFlag,1,MPI_INT,allFlags,1,MPI_INT,MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		syncTimer += (t2-t1);

		for(i=0;i<process_num;i++){
			if(allFlags[i]==1){
				termination = 0;
				break;
			}
		}
		if(termination)
			break;

	}
	
	cpu2 = MPI_Wtime();
	cpuTime = (cpu2-cpu1) - (syncTimer+commuTimer) + bcastTime ;

	//MPI_Gather(&father,1,MPI_INT,&father_array[0],process_num,MPI_INT,0,MPI_COMM_WORLD);
	int* father_array = (int*)malloc(sizeof(int)*(total_vertex+1)); 
	double* syncArray = (double*)malloc(sizeof(double)*total_vertex); 
	double* commuArray = (double*) malloc(sizeof(double)*total_vertex);
	double* cpuArray = (double*) malloc(sizeof(double)*total_vertex);

	double qq[4];
	if(rank==0){

			//int* father_array = (int*) malloc(sizeof(int)*process_num); 
			father_array[0] = father;

			t1 = MPI_Wtime();
			for(i=1;i<process_num;i++){
				//printf("%d\n",i);
				MPI_Recv(&qq[0],4,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
				father_array[i] = (double)qq[0];
				syncArray[i] = qq[1];
				commuArray[i] = qq[2];
				cpuArray[i] = qq[3];
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
			fprintf(fo,"commu:%lf \n",commuTimer);

			t1 = MPI_Wtime();
			fclose(fo);
			t2 = MPI_Wtime();
			ioTimer += (t2-t1);

	}
	else{
		qq[0] = (double) father;
		qq[1] = syncTimer;
		qq[2] = commuTimer;
		qq[3] = cpuTime;
		t1 = MPI_Wtime();
		MPI_Send(&qq[0],4,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		commuTimer += (t2-t1);

		sendCount++;
	}
	
	
	endTime = MPI_Wtime();
	double totalTime = endTime - startTime;
	
	if(rank==0){
	double CPU , SYNC, COMMU;
	CPU = 0;
	SYNC = 0;
	COMMU = 0;
	for(i=0;i<process_num;i++){
		CPU += cpuArray[i];
		SYNC += syncArray[i];
		COMMU += commuArray[i];
	}
	CPU = CPU / (double) process_num;
	SYNC = SYNC / (double) process_num;
	COMMU = COMMU / (double) process_num;

	
		fprintf(stderr,"rank:%d |send count:%d | recv count: %d |total Time:%lf |compute Time:%lf |commuTime:%lf |Sync Time:%lf |IO Time:%lf\n",rank,sendCount,recvCount,totalTime,CPU ,COMMU,SYNC,ioTimer);
	}
	else{
			fprintf(stderr,"rank:%d |send count:%d | recv count: %d |total Time:%lf |compute Time:%lf |commuTime:%lf |Sync Time:%lf |IO Time:%lf\n",rank,sendCount,recvCount,totalTime,cpuTime ,commuTimer,syncTimer,ioTimer);

	}
	MPI_Finalize();

	return 0;
}