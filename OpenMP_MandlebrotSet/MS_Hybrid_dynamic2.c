 /* 
   Sequential Mandelbrot set
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef struct complextype
{
	double real, imag;
} Compl;


int main(int argc,char** argv)
{

    int rank,processNum;
    int threadNum;
    double left,right,bottom,top,resX,resY;
    double time1,time2;
    int startX ,  widthCount,i,j;
    int numX,numY;
    int xEnable;
    int* store;
    int* pixel;
    int chunkSize = 10;
    int processPointNum = 0;
    MPI_Init(&argc,&argv);

    time1 = MPI_Wtime();
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD,&processNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

 //   slaveNum = processNum-1;
    threadNum = atoi(argv[1]);


    left = atof(argv[2]);
    right = atof(argv[3]);
    bottom = atof(argv[4]);
    top = atof(argv[5]);
    numX = atoi(argv[6]);
    numY = atoi(argv[7]);
    resX = (right-left)/(double)numX;
    resY = (top-bottom)/(double)numY;
    if(!strcmp(argv[8],"enable")){
    	xEnable = 1;
    }
    else if(!strcmp(argv[8],"disable")){
    	xEnable = 0;
    }
    else{
    	return 0 ;
    }


     //initialization for a window
	Display *display = NULL;
	Window window = NULL;  
	int screen; 

	//which screen 

	GC gc = NULL;
	XGCValues values;
	long valuemask = 0;
	int width = numX;
	int height = numY;

		/* set window position */
	int x = 0;
	int y = 0;
		/* border width in pixels */
	int border_width = 0;


    if(xEnable&&rank==0){

		    
		/* open connection with the server */ 
		    display = XOpenDisplay(NULL);
		    if(display == NULL) {
			    fprintf(stderr, "cannot open display\n");
			    return 0;
		    }

		    screen = DefaultScreen(display);
		/* set window size */
		  

		/* create window */
		    window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width,
						BlackPixel(display, screen), WhitePixel(display, screen));
		/* create graph */

		
		    gc = XCreateGC(display, window, valuemask, &values);
		//XSetBackground (display, gc, WhitePixel (display, screen));
		    XSetForeground (display, gc, BlackPixel (display, screen));
		    XSetBackground(display, gc, 0X0000FF00);
		    XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
		
		/* map(show) the window */
		    XMapWindow(display, window);
		    XSync(display, 0);
	}


	printf("processNum:%d\n",processNum);
  	
  
  	startX = rank;
  	widthCount = 0;
  
	if(rank<numX%processNum){
		widthCount = (numX/processNum)+1;
	}
	else{
		widthCount = (numX/processNum);
	}

	
	store = (int*) malloc(sizeof(int)*numY*widthCount);

		#pragma omp parallel num_threads(threadNum) private(i, j) 
		{
			double start = omp_get_wtime();
			int threadPointNum = 0;
			#pragma omp for schedule(dynamic,chunkSize) nowait
		    for(i=0;i<numY; i++) {
			    for(j=startX; j<numX; j+=processNum) {
			    	Compl z, c;
		  			int repeats;
		  			double temp, lengthsq;
				    z.real = 0.0;
				    z.imag = 0.0;
				    c.real =  left + (double)j * resX;//((double)i - 400.0)/200.0; /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
				    c.imag =  bottom + (double)i * resY;//((double)j - 400.0)/200.0; /* So needs to scale the window */
				    repeats = 0;
				    lengthsq = 0.0;

				    threadPointNum++;

				    while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
					    temp = z.real*z.real - z.imag*z.imag + c.real;
					    z.imag = 2*z.real*z.imag + c.imag;
				    	z.real = temp;
				    	lengthsq = z.real*z.real + z.imag*z.imag; 
			    		repeats++;
		    		}
		   			store[((j-startX)/processNum*numY)+i] = repeats;

			    }
		    }

		    #pragma omp critical
		    {
		    	processPointNum += threadPointNum;
		    	double endTime  = omp_get_wtime();
		    	double elapseTime = (double)(endTime-start);
		    	printf("Process#:%d  thread#:%d  threadTime:%lf threadNumPoint:%d \n",rank,omp_get_thread_num(),elapseTime,threadPointNum);

		    }
		}


	startX = 0;
	//tempCount = 0;
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

	if(xEnable){
		int readCount = 0;
		int iStart = 0;
		if(rank==0){
			for(i=0;readCount<numX*numY;i+=processNum){

				if(i>=numX){
					i = ++iStart;
				}

				for(j=0;j<numY;j++){
					XSetForeground (display, gc,  1024 * 1024 * (pixel[readCount++] % 256));		
					XDrawPoint (display, window, gc, i, j);
				}


			}
			XFlush(display);
		}
	}
	




	time2 = MPI_Wtime();
	printf("rank:%d  processPoint#:%d time:%lf\n",rank,processPointNum,time2-time1);
	MPI_Finalize();
	//XFlush(display);
	if(xEnable){
		sleep(5);
	}
	return 0;
}
