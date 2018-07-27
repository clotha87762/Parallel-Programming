#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <queue>
#include <pthread.h>
#include <vector>
using namespace std;

int thread_num;
int source_id;
int total_vertex;
int edge_num;
int* dist_array;
int* father_array;
queue<int> work_queue;
int** weight;
pthread_mutex_t queue_mutex;
int* flags;
bool* nowInQueue;
vector<int>* outList;

void* FindPath(void* id){

	int thread_id = *((int*) id);
	int nowVertex; 
	int i,j;
	bool breakChain = false;
	vector<int> updateList;
	//bool flag = false;
	//printf("thread id:%d\n",thread_id);

	while(1){
		
		nowVertex = -1;
		updateList.clear();
		//printf("aaa %d\n",thread_id);
		do{

				//printf("bbb %d\n",thread_id);

				if(!work_queue.empty()){

					pthread_mutex_lock(&queue_mutex);
					//printf("get lock %d\n",thread_id);
					if(!work_queue.empty()){
						
						nowVertex = work_queue.front();
						//printf("get vertex %d\n",nowVertex);
						work_queue.pop();
						nowInQueue[nowVertex] = false;
					
					}
					else{
						flags[thread_id] = 1;
					}
					pthread_mutex_unlock(&queue_mutex);
				}
				else{
					flags[thread_id] = 1;
				}
				
				 
				//printf("unlock %d \n",thread_id);
		

			for(i=0;i<thread_num;i++){

				if(flags[i]==0)
					break;
				if(i==thread_num-1){
					//printf("end1\n");
					return(NULL);
					printf("end2\n");
				}
				//printf("ffflag %d \n",thread_id);
			}

		

		}while(nowVertex<0);


		flags[thread_id] = 0;


		 for (vector<int>::iterator it = outList[nowVertex].begin() ; it != outList[nowVertex].end(); ++it){

		 	i = *it;
			//printf("qq w:%d\n",weight[nowVertex][i]);
			if(i==source_id)
				continue;
			//printf("a=%d  b=%d  \n",dist_array[nowVertex]+weight[nowVertex][i],dist_array[i] );
			if(dist_array[nowVertex]+weight[nowVertex][i]<dist_array[i]){  // 試要不要lock dis array
				//printf("OAO\n");
				updateList.push_back(i);
			}

		
		}
		pthread_mutex_lock(&queue_mutex);
		 for (vector<int>::iterator it = updateList.begin() ; it != updateList.end(); ++it){
						i = *it;
						if(dist_array[nowVertex]+weight[nowVertex][i]<dist_array[i]){
							
							dist_array[i] = dist_array[nowVertex] + weight[nowVertex][i];
							father_array[i] = nowVertex;
							
							if(!nowInQueue[i]){
								//printf("push%d\n",i);
								work_queue.push(i);
								nowInQueue[i] = true;
							}
						
						}
		}
		pthread_mutex_unlock(&queue_mutex);
		//printf("bbbid:%d\n",thread_id);
	}

	pthread_exit(NULL);

}

int main(int argc,char* argv[]){

	thread_num = atoi(argv[1]);
	source_id = atoi(argv[4]) - 1;
	//FILE* input = fopen(argv[2],"r");
	ifstream input;
	input.open(argv[2]);
	input >> total_vertex;
	input >> edge_num;

	//printf("total vertex %d\n",total_vertex);

	dist_array = new int[total_vertex];
	father_array = new int[total_vertex];
	weight = new int*[total_vertex];

	flags = new int[thread_num];
	nowInQueue = new bool[total_vertex];

	outList = new vector<int>[total_vertex];



	for(int i=0;i<total_vertex;i++)
		weight[i] = new int[total_vertex];

	for(int i=0;i<total_vertex;i++)
		for(int j=0;j<total_vertex;j++)
			weight[i][j] = -1;

	for(int i=0;i<total_vertex;i++){
		dist_array[i] = 210000000;
		father_array[i] = -1;
		nowInQueue[i] = false;
		weight[i][i] = 0;
	}
	for(int i=0;i<thread_num;i++){
		flags[i]   = 0;
	}
	dist_array[source_id] = 0;
	nowInQueue[source_id] = true;
	father_array[source_id] = source_id;

	int xx,yy,X,Y;
	for(X=0;X<total_vertex*10;X++)
		for(Y=0;Y<total_vertex*10;Y++){
			xx = X;
			yy = Y;
		}

	int from,to;
	for(int i=0;i<edge_num;i++){
		input>>from;
		input>>to;
		input>>weight[from-1][to-1];
		weight[to-1][from-1] = weight[from-1][to-1];
		outList[from-1].push_back(to-1);
		outList[to-1].push_back(from-1);
	}

	work_queue.push(source_id);

	pthread_mutex_init(&queue_mutex,NULL);





	pthread_t* threads;
	threads = new pthread_t[thread_num];

	int *k = new int[thread_num];
	for(int i=0;i<thread_num;i++){
		k[i] = i; 
		pthread_create(&threads[i],NULL,FindPath,(void*)&k[i]);  // 這邊到底要不要加&啊啊啊
	}

	//while(1);

	for(int i=0;i<thread_num;i++){
		pthread_join(threads[i],NULL);
		//printf("join%d\n",i);
	}

	//pthread_exit(NULL);


	//printf("111\n");



	ofstream output;
	output.open(argv[3]);
	int* tempArray = new int[total_vertex];
	int x,count;
	for(int i=1;i<total_vertex+1;i++){
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

	//printf("222\n" );

	delete[] dist_array;
	delete[] father_array;
	for(int i=0;i<total_vertex;i++){
		delete[] weight[i];
	}
	delete[] weight;
	delete[] tempArray;
	delete[] flags;
	delete[] nowInQueue;
	delete[] threads;
	pthread_mutex_destroy(&queue_mutex);

	return 0;
}
