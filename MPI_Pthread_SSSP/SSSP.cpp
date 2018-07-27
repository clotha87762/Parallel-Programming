#include <pthread.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>

int thread_num;
int source_id;
int total_vertex;
int edge_num;


int main(int argc,char*[] argv){

	//thread_num = atoi(argv[1]);
	//source_id = atoi(argv[4]);
	//FILE* input = fopen(argv[2],"r");
	ifstream input;
	input.open("input.txt");
	input >> total_vertex;
	input >> edge_num;
	cout<<total_vertex<<" "<<edge_num<<endl;


	pthread_t threads[]

	return 0;
}
