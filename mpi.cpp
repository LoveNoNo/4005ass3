#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <mpi.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"


int n_body;
int n_iteration;
int size_l; // data size of slave thread
std::chrono::duration<double> time_total;

int my_rank;
int world_size;
//Logger l = Logger("sequential", n_body, bound_x, bound_y);


void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}


void update_position(double *x, double *y, double *xr, double *yr, double *vx, double *vy, int n, int rank) {
    //TODO: update position 
    for(int i=0; i<n; i++){
        int pos = (rank-1)*size_l + i; // x&y use pos as total index

        xr[i] = x[pos] + vx[i] * dt;
        yr[i] = y[pos] + vy[i] * dt;

        //collsion detection: check new pos, if collison, modify speed and cal right new pos

        // ball & ball collision
        for (int j=0; j<n_body; j++){ // MPI Notice: this also change to n-body
            if(j != i){
                
                if(pow(x[pos]-x[j],2) <= 2 && pow(y[pos]-y[j],2) <=2 ){ //distance less than 2R
                    // velocity & pos modify by collision in two cases
                    if( pow(x[pos]-x[j],2) <= pow(y[pos]-y[j],2)){ // more like col up & down
                        // back to origin pos
                        yr[i] = y[pos] - vy[i] * dt;

                        // inverse velocity
                        vy[i] = -vy[i]/2; 
                                          // MPI Notice! Here not change velocity of body j to avoid mesg pass

                        // re-calculate pos through inversed velocity
                        yr[i] = y[pos] + vy[i] * dt;
                                          // MPI Notice! Also here not change body j data

                    }else{ // more like left & right col
                        xr[i] = x[pos] - vx[i] * dt;

                        vx[i] = -vx[i]/2; // only inverse velocity X

                        xr[i] = x[pos] + vx[i] * dt;
                    }
                }
            }
        }

        // ball & wall collision
        if(xr[i] <= 0 || x[i] >= bound_x){ // check pos with boundary
            xr[i] = x[pos] - vx[i] * dt; // back to origin pos
            vx[i] = -vx[i]/2;         // inverse velocity X
            xr[i] = x[pos] + vx[i] * dt; // re-calculate pos X
        }
        if(yr[i] <=0 || y[i] >= bound_y){
            yr[i] = y[pos] - vy[i] * dt;
            vy[i] = -vy[i]/2;         // half velocity to avoid very high speed after col bound
            yr[i] = y[pos] + vy[i] * dt;
        }
        
    }

}


void update_velocity_MPI(double *m, double *x, double *y, double *vx, double *vy, int n, int rank) {
    //TODO: calculate force and acceleration, update velocity
    for(int i=0; i<n; i++){ //遍历所有body: Notice: every thread cal part of body
        double dvx = 0.0f; 
        double dvy = 0.0f;
        double dv_temp = 0.0f;
        int pos = (rank-1)*size_l + i; // confirm the data pos this thread cal in total data array 

        //求和计算所有其它body对此body的影响
        for(int j=0; j<n_body; j++){  // MPI Notice: here still sum all body F, so j < total n 
            if(j != pos){ 
                dv_temp = gravity_const * m[j] * dt / (pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2)); //向量 dv，还需要分解到x和y两个方向的dv
                dvx = dvx + dv_temp * ((x[j]-x[pos]) / sqrt(pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2))); 
                dvy = dvy + dv_temp * ((y[j]-y[pos]) / sqrt(pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2)));
            }           
        }
        
        vx[i] = vx[i] + dvx;
        vy[i] = vy[i] + dvy;
    }
}


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    size_l = n_body/(world_size-1);
    if(n_body % (world_size-1) != 0){
        printf("ERROR! Thread number minus one must be exact devide bu body mumber!");
    }

    int sendcounts[world_size];
    int recvcounts[world_size];
    int displs[world_size];

    for(int i=0; i<world_size; i++){
        if(i == 0){
            sendcounts[i] = 0;
            recvcounts[i] = 0;
            displs[i] = 0;
        }else{
            sendcounts[i] = size_l;
            recvcounts[i] = size_l;
            displs[i] = (i-1) * size_l;
        }
    }

    double* total_m = new double[n_body];
    double* total_x = new double[n_body];
    double* total_y = new double[n_body];
    double* total_vx = new double[n_body];
    double* total_vy = new double[n_body];

    double* local_m = new double[n_body]; // follow three are total data copy
    double* local_x = new double[n_body];
    double* local_y = new double[n_body];
    double* local_vx = new double[size_l]; // these two are local cal body data result
    double* local_vy = new double[size_l];

    double* local_all_vx = new double[n_body];
    double* local_all_vy = new double[n_body];

    double* local_result_x = new double[size_l];
    double* local_result_y = new double[size_l];


    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();


	if (my_rank == 0) {
		#ifdef GUI
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(500, 500); 
		glutInitWindowPosition(0, 0);
		glutCreateWindow("N Body Simulation MPI Implementation");
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, bound_x, 0, bound_y);
		#endif

        generate_data(total_m, total_x, total_y, total_vx, total_vy, n_body);

        //Logger l = Logger("sequential", n_body, bound_x, bound_y);

	} 
    MPI_Barrier(MPI_COMM_WORLD);

    for (int iter = 0; iter < n_iteration; iter++){
        if(my_rank == 0){

            // ====== MPI Master ======

            //std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            for(int threadID=1; threadID<world_size; threadID++){
                if(iter == 0){
                    MPI_Send(total_m, n_body, MPI_DOUBLE, threadID, 1, MPI_COMM_WORLD); 
                }
                MPI_Send(total_x, n_body, MPI_DOUBLE, threadID, 2, MPI_COMM_WORLD);
                MPI_Send(total_y, n_body, MPI_DOUBLE, threadID, 3, MPI_COMM_WORLD);
            }   

            // send velocity to thread One
            MPI_Send(total_vx, n_body, MPI_DOUBLE, 1, 4, MPI_COMM_WORLD);
            MPI_Send(total_vy, n_body, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);

            MPI_Recv(total_x, n_body, MPI_DOUBLE, 1, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(total_y, n_body, MPI_DOUBLE, 1, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // TODO End

            //std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
            //std::chrono::duration<double> time_span = t2 - t1;
            //time_total = time_total + time_span;

            //printf("Iteration %d, elapsed time: %.3f\n", iter, time_span);

            //l.save_frame(total_x, total_y);

            #ifdef GUI
            glClear(GL_COLOR_BUFFER_BIT);
            glColor3f(1.0f, 0.0f, 0.0f);
            glPointSize(2.0f);
            glBegin(GL_POINTS);
            double xi;
            double yi;
            for (int i = 0; i < n_body; i++){
                xi = total_x[i];
                yi = total_y[i];
                glVertex2f(xi, yi);
            }
            glEnd();
            glFlush();
            glutSwapBuffers();
            #else

            #endif

        }else{

            // ======= MPI Slave ========

            if(iter == 0){
                MPI_Recv(local_m, n_body, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // only once
            }

            MPI_Recv(local_x, n_body, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(local_y, n_body, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if(my_rank == 1){   
                MPI_Recv(local_all_vx, n_body, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(local_all_vy, n_body, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   
            }

            // scatter velocity parts to every slave thread
            MPI_Scatterv(local_all_vx, sendcounts, displs, MPI_DOUBLE, local_vx, sendcounts[my_rank], MPI_DOUBLE, 1, MPI_COMM_WORLD);
            MPI_Scatterv(local_all_vy, sendcounts, displs, MPI_DOUBLE, local_vy, sendcounts[my_rank], MPI_DOUBLE, 1, MPI_COMM_WORLD);

            update_velocity_MPI(local_m, local_x,local_y,local_vx, local_vy, size_l, my_rank);

            update_position(local_x, local_y, local_result_x, local_result_y, local_vx, local_vy, size_l, my_rank);
            
            MPI_Barrier(MPI_COMM_WORLD);

            //gather data part from each thread after cal
            MPI_Gatherv(local_vx, size_l, MPI_DOUBLE, local_all_vx, recvcounts, displs, MPI_DOUBLE, 1, MPI_COMM_WORLD);        
            MPI_Gatherv(local_vy, size_l, MPI_DOUBLE, local_all_vy, recvcounts, displs, MPI_DOUBLE, 1, MPI_COMM_WORLD);                   
            
            MPI_Gatherv(local_result_x, size_l, MPI_DOUBLE, local_x, recvcounts, displs, MPI_DOUBLE, 1, MPI_COMM_WORLD);        
            MPI_Gatherv(local_result_y, size_l, MPI_DOUBLE, local_y, recvcounts, displs, MPI_DOUBLE, 1, MPI_COMM_WORLD);        
            
            // send data back to master
            if(my_rank == 1){   
                MPI_Send(local_x, n_body, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
                MPI_Send(local_y, n_body, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);   
            }   
        }
    }

    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = t2 - t1;

    if(my_rank == 0){
        delete[] total_m;
        delete[] total_x;
        delete[] total_y;
        delete[] total_vx;
        delete[] total_vy;

        delete[] local_m;
        delete[] local_x;
        delete[] local_y;
        delete[] local_vx; 
        delete[] local_vy;

        delete[] local_all_vx ;
        delete[] local_all_vy;

        delete[] local_result_x;
        delete[] local_result_y;
    }

	if (my_rank == 0){
		printf("Student ID: 1170103321\n"); // replace it with your student id
		printf("Name: XU Jiale\n"); // replace it with your name
        printf("MPI Total time cost %.3f seconds\n", time_span);
		printf("Assignment 2: N Body Simulation MPI Implementation\n");
	}

	MPI_Finalize();

	return 0;
}

