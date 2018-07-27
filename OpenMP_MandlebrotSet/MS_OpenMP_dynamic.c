 /* 
   Sequential Mandelbrot set
 */

#include <X11/Xlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
typedef struct complextype
{
	double real, imag;
} Compl;


int main(int argc,char** argv)
{

   // int rank,processNum,slaveNum;
	int threadNum;
    double left,right,bottom,top,resX,resY;
    double time1,time2;
    int numX,numY;
    int xEnable;
    int chunkSize = 200;

  //  MPI_Init(&argc,&argv);

   // time1 = MPI_Wtime();
   // MPI_Status status;

   // MPI_Comm_size(MPI_COMM_WORLD,&processNum);
   // MPI_Comm_rank(MPI_COMM_WORLD,&rank);

   // slaveNum = processNum-1;

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



    
    double start = omp_get_wtime();
    double end;

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


    if(xEnable){

		    
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


	//printf("processNum:%d\n",processNum);


		/* draw points */
		
		    int i, j;
		    //int* received = (int*) malloc(sizeof(int));
		    //int* received = (int*) malloc(sizeof(int)*2);
		    //int* toSend = (int*) malloc(sizeof(int)*3);
		    //MPI_Recv(received,2,MPI_INT,0,rank,MPI_COMM_WORLD,&status);

		#pragma omp parallel num_threads(threadNum) private(i, j)
		{
			int pointNum = 0;

			#pragma omp for schedule(dynamic,chunkSize) nowait
		    for(i=0;i<numY; i++) {
			    for(j=0; j<numX; j++) {
			    	Compl z, c;
		  			int repeats;
		  			double temp, lengthsq;
				    z.real = 0.0;
				    z.imag = 0.0;
				    c.real =  left + (double)j * resX;//((double)i - 400.0)/200.0; /* Theorem : If c belongs to M(Mandelbrot set), then |c| <= 2 */
				    c.imag =  bottom + (double)i * resY;//((double)j - 400.0)/200.0; /* So needs to scale the window */
				    repeats = 0;
				    lengthsq = 0.0;

				    pointNum++;

				    while(repeats < 100000 && lengthsq < 4.0) { /* Theorem : If c belongs to M, then |Zn| <= 2. So Zn^2 <= 4 */
					    temp = z.real*z.real - z.imag*z.imag + c.real;
					    z.imag = 2*z.real*z.imag + c.imag;
				    	z.real = temp;
				    	lengthsq = z.real*z.real + z.imag*z.imag; 
			    		repeats++;
		    		}
		   
					//MPI_Send(toSend,3,MPI_INT,0,0,MPI_COMM_WORLD);    
				    if(xEnable){
				    	#pragma omp critical
				    	{
				    		XSetForeground (display, gc,  1024 * 1024 * (repeats % 256));		
				    		XDrawPoint (display, window, gc, j, i);
			    		}
			    	}
			    }
		    }

		    #pragma omp critical
		    {
		    	int endTime;
				endTime = omp_get_wtime();
				double elapseTime = (double)(endTime-start);

				printf("pointNum: %d thread#:%d  Time: %lf \n",pointNum,omp_get_thread_num(),elapseTime);
		    }
		}

	     
	


		//end = omp_get_wtime();
		//double elapseTime = (double)(end-start);

		//printf("Time: %lf \n",elapseTime);
	//time2 = MPI_Wtime();
	//printf("rank:%d  time:%lf\n",rank,time2-time1);
	//MPI_Finalize();
	//XFlush(display);
	if(xEnable){
		XFlush(display);	
		sleep(5);	
	}

	return 0;
}
