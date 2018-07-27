

if(this process has fanout id < its ouwn id)
	reactive = true;

 

 rootflag= 1;

 if(rank==source id){
 	// MPI_Isend new dist to all fanout
 	if(rank==0)
 		send white token to process 1

 }

while(1){
	
	MPI_Recv()

	if(rank!=0){

		if(recv is dist){
			update distance , isend to all outgoing
			if(reactive)
				black = 1;
			if (flag){
				pass token to process x+1
				black = 0;
				set flag = false
			}
		}
		else if(recv is token){
			if(relaxed)
				send token to process x+1
				black = 0;
			else
				set flag = true
		}
		else{
			terminate();
		}
	}
	else{

		if(recv is dist){
			update distance , isend to all outgoing
			
			if (first time&&rank!=source){
				send white token to process 1
			}
		}
		else if(recv is token){
			if(recv white token)
				send terminate to all process;
				terminate();
			
			else{
				send white token to process 1;
			}	
		}

	

	}

	

}