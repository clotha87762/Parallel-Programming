
		
		startX = rank;
		if(rank<numX%processNum){
			widthCount = (numX/processNum)+1;
		}
		else{
			widthCount = (numX/processNum);
		}

	
		store = (int*) malloc(sizeof(int)*numY*widthCount);

		#pragma omp parallel num_threads(threadNum) private(i, j) 
		{
			#pragma omp for schedule(static) nowait
		    for(i=0;i<numY; i++) {
			    for(j=startX; j<numX; j+=processNum) {
			    	// Mandelbrot set calculation and use "store" array to store the result
			    }
		    }

		}

		if(rank==0){
			pixel = (int*)malloc(sizeof(int)*numX*numY);
		}
			int* displace = (int*) malloc(sizeof(int)*processNum);
			int* recvCount = (int*) malloc(sizeof(int)*processNum);
			displace[0] = 0;

		for(i=0;i<processNum;i++){
			if(i<numX%processNum){
				recvCount[i] = ((numX/processNum)+1) * numY;
			}
			else{
				recvCount[i] = (numX/processNum) * numY;
			}
			if(i!=processNum-1){
				displace[i+1] = displace[i] + recvCount[i];
			}
		}

	MPI_Gatherv(store,widthCount*numY,MPI_INT,pixel,recvCount,displace,MPI_INT,0,MPI_COMM_WORLD);