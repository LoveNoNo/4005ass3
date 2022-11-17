#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <pthread.h>

#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_thd; // number of threads

int n_body;
int n_iteration;
int size_l;

double* m;
double* x;
double* y;
double* vx;
double* vy;
pthread_mutex_t mutex_v;
std::chrono::duration<double> time_total;

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



void update_position(double *x, double *y, double *vx, double *vy, int n, int td_id) {

    for(int i=0; i<n; i++){
        int pos = (td_id-1)*size_l + i; // x&y use pos as total index

        x[pos] = x[pos] + vx[i] * dt;
        y[pos] = y[pos] + vy[i] * dt;

        //collsion detection: check new pos, if collison, modify speed and cal right new pos

        // ball & ball collision
        for (int j=0; j<n_body; j++){ 
            if(j != i){
                // because not modify j in otehr thread, no data race
                if(pow(x[pos]-x[j],2) <= 2 && pow(y[pos]-y[j],2) <=2 ){ //distance less than 2R

                    if( pow(x[pos]-x[j],2) <= pow(y[pos]-y[j],2)){ // more like col up & down

                        y[pos] = y[pos] - vy[pos] * dt;

                        vy[pos] = -vy[pos]/2; 
                                          // Pthread Notice! Here not change velocity of body j to avoid data race

                        y[pos] = y[pos] + vy[pos] * dt;

                    }else{ // more like left & right col
                        x[pos] = x[pos] - vx[pos] * dt;

                        vx[pos] = -vx[pos]/2;

                        x[pos] = x[pos] + vx[pos] * dt;
                    }
                }
            }
        }

        // ball & wall collision
        if(x[pos]<=0 || x[pos]>= bound_x){// check pos with boundary
            x[pos] = x[pos] - vx[pos] * dt; // back to origin pos
            vx[pos] = -vx[pos]/2;         // inverse velocity X
            x[pos] = x[pos] + vx[pos] * dt; // re-calculate pos X
        }
        if(y[pos] <=0 || y[pos] >= bound_y){
            y[pos] = y[pos] - vy[pos] * dt;
            vy[pos] = -vy[pos]/2;         // half velocity to avoid very high speed after col bound
            y[pos] = y[pos] + vy[pos] * dt;
        }
        
    }

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n, int td_id) {
    
    for(int i=0; i<n; i++){
        double dvx = 0.0f; 
        double dvy = 0.0f;
        double dv_temp = 0.0f;
        int pos = (td_id-1)*size_l + i; // confirm the data pos this thread cal in total data array 

        for(int j=0; j<n_body; j++){  // MPI Notice: here still sum all body F, so j < total n 
            if(j != pos){ 
                dv_temp = gravity_const * m[j] * dt / (pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2)); //向量 dv，还需要分解到x和y两个方向的dv
                dvx = dvx + dv_temp * ((x[j]-x[pos]) / sqrt(pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2))); 
                dvy = dvy + dv_temp * ((y[j]-y[pos]) / sqrt(pow(x[j]-x[pos],2) + pow(y[j]-y[pos],2)));
            }           
        }
        
        // Every thread write own body list, will not occur data race
        vx[pos] = vx[pos] + dvx;
        vy[pos] = vy[pos] + dvy;
    }

}


typedef struct {
    //TODO: specify your arguments for threads
    int thd_id;
    double* m;
    double* x;
    double* y;
    double* vx;
    double* vy;
    //TODO END
} Args;


void* worker(void* args) {
    //Pthread TODO: procedure in each threads

    Args* my_arg = (Args*) args;
    int tid = my_arg->thd_id;
    // double* m = my_arg->m;
    // double* x = my_arg->x;
    // double* y = my_arg->y;
    // double* vx = my_arg->vx;
    // double* vy = my_arg->vy;

    update_velocity(m, x, y, vx, vy, size_l, tid);
    update_position(x, y, vx, vy, size_l, tid);


    // TODO END
}


void master(){
    m = new double[n_body];
    x = new double[n_body];
    y = new double[n_body];
    vx = new double[n_body];
    vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    pthread_t thds[n_thd];
    pthread_mutex_init(&mutex_v, NULL);

    Args args[n_thd];
    for (int thd = 0; thd < n_thd; thd++){
        args[thd].thd_id = thd;
        // args[thd].m = m;
        // args[thd].x = x;
        // args[thd].y = y;
        // args[thd].vx = vx;
        // args[thd].vy = vy;
    }

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        
        // PThread TODO: assign jobs
        
        for (int thd = 0; thd < n_thd; thd++) pthread_create(&thds[thd], NULL, worker, &args[thd]);
        for (int thd = 0; thd < n_thd; thd++) pthread_join(thds[thd], NULL);

        //TODO End

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = t2 - t1;
        time_total = time_total + time_span;

        printf("Iteration %d, elapsed time: %.3f\n", i, time_span);

        l.save_frame(x, y);

        #ifdef GUI
        glClear(GL_COLOR_BUFFER_BIT);
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(2.0f);
        glBegin(GL_POINTS);
        double xi;
        double yi;
        for (int i = 0; i < n_body; i++){
            xi = x[i];
            yi = y[i];
            glVertex2f(xi, yi);
        }
        glEnd();
        glFlush();
        glutSwapBuffers();
        #else

        #endif
    }

    printf("\nStudent ID: 117010332\n"); // replace it with your student id
    printf("Name: XU Jiale\n"); // replace it with your name
    printf("Total time cost %.3f seconds\n", time_total);
    printf("Assignment 2: N Body Simulation Pthread Implementation\n");

    // delete[] m;
    // delete[] x;
    // delete[] y;
    // delete[] vx;
    // delete[] vy;

    pthread_mutex_destroy(&mutex_v);
}


int main(int argc, char *argv[]) {
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);
    n_thd = atoi(argv[3]);
    size_l = n_body/(n_thd-1);

    if(n_body % (n_thd-1) != 0){
        printf("ERROR! Thread number minus one must be exact devide bu body mumber!");
    }

    #ifdef GUI
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Pthread");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

	return 0;
}

