

global distance array[];
global working_queue;

find_path(){

	while(1){
	
		while(vertex!=null){
			
			{ 	// critical section
				get a vertex from queue
			}
			flag[this_thread] = 1;
			if(flags of all threads are 1, then return)
		}

		flag[this_thread]=0;

		{
			for(){
			calculate distance
			
				{ // critical 
				update distance from this vertex
				put the vertex into queue
				}
			}
		}
	}

}


main(){
	
	for(i=0;i<thread_num;i++)
	create_threads(find_path)

	for(i=0;i<thread_num;i++)
	thread_join(i)
}