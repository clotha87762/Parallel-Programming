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


int DATA_TAG = 0;
int END_TAG = 1;

int main(int argc,char** argv)
{

    int rank,processNum,slaveNum,threadNum;
    double left,right,bottom,top,resX,resY;
    double time1,time2;
    int numX,numY;
    int xEnable;

    MPI_Init(&argc,&argv);

    time1 = MPI_Wtime();
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD,&processNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    slaveNum = processNum-1;

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
	int processPointNum = 0;
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

	//int* pixel = (int*) malloc(sizeof(int)*);
	printf("processNum:%d\n",processNum);
    if(processNum==1&&rank==0){
    		

		#pragma omp parallel num_threads(threadNum) 
		{
			Compl z, c;
		    int repeats;
		    double temp, lengthsq;
		    int i, j;

			#pragma omp for schedule(dynamic,1) 
		    for(i=0;i<numX; i++) {
			    for(j=0; j<numY; j++) {
				    z.real = 0.0;
				    z.imag = 0.0;
				    c.real =  left + (double)i * resX;//((double)i - 400.0)/200.0; /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
				    c.imag =  bottom + (double)j * resY;//((double)j - 400.0)/200.0; /* So needs to scale the window */
				    repeats = 0;
				    lengthsq = 0.0;

				    while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
					    temp = z.real*z.real - z.imag*z.imag + c.real;
					    z.imag = 2*z.real*z.imag + c.imag;
				    	z.real = temp;
				    	lengthsq = z.real*z.real + z.imag*z.imag; 
			    		repeats++;
		    		}
 					if(xEnable){
 						#pragma omp critical
 						{
						    XSetForeground (display, gc,  1024 * 1024 * (repeats % 256));		
						  	XDrawPoint (display, window, gc, i, j);
			    		}
			    	}
			    }
		    }
		}

		    if(xEnable&&rank==0){
		        XFlush(display);
		      //  sleep(5);
	    	}
    }
    else{
	    if(rank==0){

	    	
			int row;
			int i;
			row=0;
			i=0;
			int count = 0;
			//int* toSend = (int*) malloc(sizeof(int)*2);
			int toSend;
			int* received = (int*) malloc(sizeof(int)*(numX+2)); // 0-> rank  1->row other->iter
			for(i=0;i<slaveNum;i++){
				MPI_Send(&row,1,MPI_INT,i+1,i+1,MPI_COMM_WORLD);
				row++;
				count++;
			}

			do{
				MPI_Recv(received,(numX+2),MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
				count--;
				if(row<numY){
					toSend = row;
					MPI_Send(&toSend,1,MPI_INT,received[1],received[1],MPI_COMM_WORLD);
					row++;
					count++;
				}
				else{
					toSend = -1;
					MPI_Send(&toSend,1,MPI_INT,received[1],received[1],MPI_COMM_WORLD);
				}

				if(xEnable){
					for(i=0;i<numX;i++){
						XSetForeground (display, gc,  1024 * 1024 * (received[2+i] % 256));		
					    XDrawPoint (display, window, gc, i, received[0]);
					}
				}

			}while(count>0);

			if(xEnable){
		        XFlush(display);
	    	}

	    }
	    else{
			int i;
			int received;
			int* toSend ;
		/* draw points */
			

			    //int* received = (int*) malloc(sizeof(int));
			    //int* received = (int*) malloc(sizeof(int)*2);
			    
			toSend = (int*) malloc(sizeof(int)*(numX+2));
			double* timeCounter = (double*) malloc(sizeof(double)*threadNum);
			int* pointCounter = (int*) malloc(sizeof(int)*threadNum);
			for(i=0;i<threadNum;i++){
				timeCounter[i] = 0;
				pointCounter[i] = 0;
			}

			while(1){
			   
			    MPI_Recv(&received,1,MPI_INT,0,rank,MPI_COMM_WORLD,&status);
			    if(received==-1)
			    	break;

				toSend[0] = received;
				toSend[1] = rank;

			#pragma omp parallel num_threads(threadNum) private(i) shared(toSend,received)
			{
				int threadPointNum = 0;
				double start = omp_get_wtime();
				#pragma omp for schedule(dynamic,1) nowait
			    for(i=0;i<numX; i++) {
				  		
				  		Compl z, c;
					    int repeats;
					    double temp, lengthsq;
					    z.real = 0.0;
					    z.imag = 0.0;
					    c.real =  left + (double)i * resX;//((double)i - 400.0)/200.0; /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
					    c.imag =  bottom + (double)received * resY;//((double)j - 400.0)/200.0; /* So needs to scale the window */
					    repeats = 0;
					    lengthsq = 0.0;

					    pointCounter[omp_get_thread_num()]++;

					    while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
						    temp = z.real*z.real - z.imag*z.imag + c.real;
						    z.imag = 2*z.real*z.imag + c.imag;
					    	z.real = temp;
					    	lengthsq = z.real*z.real + z.imag*z.imag; 
				    		repeats++;
			    		}
			    		toSend[2+i] = repeats;
				    
			     }
			  	 	
			  	 	timeCounter[omp_get_thread_num()]+=omp_get_wtime()-start;
			    	
			}
				MPI_Send(toSend,numX+2,MPI_INT,0,0,MPI_COMM_WORLD);
			 
			}
			for(i=0;i<threadNum;i++){
				printf("Process#:%d  thread#:%d  threadTime:%lf threadNumPoint:%d \n",rank,i,timeCounter[i],pointCounter[i]);
				processPointNum+=pointCounter[i];
			}

	    }   
	}

	time2 = MPI_Wtime();
	printf("rank:%d  processPoint#:%d time:%lf\n",rank,processPointNum,time2-time1);
	MPI_Finalize();
	//XFlush(display);
	if(xEnable)
	sleep(5);
	return 0;
}
