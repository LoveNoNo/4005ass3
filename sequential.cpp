#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#ifdef GUI
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "./headers/physics.h"
#include "./headers/logger.h"

int n_body;
int n_iteration;
std::chrono::duration<double> time_total;

void generate_data(double *m, double *x,double *y,double *vx,double *vy, int n) {
    // TODO: Generate proper initial position and mass for better visualization
    srand((unsigned)time(NULL));
    //srandom(n_body);
    for (int i = 0; i < n; i++) {
        m[i] = rand() % max_mass + 1.0f;
        x[i] = 2000.0f + rand() % (bound_x / 4);
        y[i] = 2000.0f + rand() % (bound_y / 4);
        vx[i] = 0.0f;
        vy[i] = 0.0f;
    }
}



void update_position(double *x, double *y, double *vx, double *vy, int n) {
    //TODO: update position: x = x0 + vx * dt

    for(int i=0; i<n; i++){
        x[i] = x[i] + vx[i] * dt;
        y[i] = y[i] + vy[i] * dt;

        //collsion detection: check new pos, if collison, modify speed and cal right new pos

        // ball & ball collision
        for (int j=0; j<n; j++){
            if(j != i){
                
                if(pow(x[i]-x[j],2) <= 2 && pow(y[i]-y[j],2) <=2 ){ //distance less than 2R
                    // velocity & pos modify by collision in two cases
                    if( pow(x[i]-x[j],2) <= pow(y[i]-y[j],2)){ // more like col up & down
                        // back to origin pos
                        y[i] = y[i] - vy[i] * dt;
                        y[j] = y[j] - vy[j] * dt;

                        // inverse velocity
                        vy[i] = -vy[i]/2; // only inverse velocity Y
                        vy[j] = -vy[j]/2; // inverse both col body

                        // re-calculate pos through inversed velocity
                        y[i] = y[i] + vy[i] * dt;
                        y[j] = y[j] + vy[j] * dt;

                    }else{ // more like left & right col
                        x[i] = x[i] - vx[i] * dt;
                        x[j] = x[j] - vx[j] * dt;

                        vx[i] = -vx[i]/2; // only inverse velocity X
                        vx[j] = -vx[j]/2;

                        x[i] = x[i] + vx[i] * dt;
                        x[j] = x[j] + vx[j] * dt;
                    }
                }
            }
        }

        // ball & wall collision
        if(x[i] <= 0 || x[i] >= bound_x){ // check pos with boundary
            x[i] = x[i] - vx[i] * dt; // back to origin pos
            vx[i] = -vx[i]/2;         // inverse velocity X
            x[i] = x[i] + vx[i] * dt; // re-calculate pos X
        }
        if(y[i] <=0 || y[i] >= bound_y){
            y[i] = y[i] - vy[i] * dt;
            vy[i] = -vy[i]/2;         // half velocity to avoid very high speed after col bound
            y[i] = y[i] + vy[i] * dt;
        }
        
    }

}

void update_velocity(double *m, double *x, double *y, double *vx, double *vy, int n) {
    //TODO: calculate force and acceleration, update velocity: v = v0 + (G m' dt) / R^2


    for(int i=0; i<n; i++){ //遍历所有body
        double dvx = 0.0f; //记录所有其它body对该body的速度影响之和
        double dvy = 0.0f;
        double dv_temp = 0.0f;

        for(int j=0; j<n; j++){ //求和计算所有其它body对此body的影响
            if(j != i){ //自己和自己无需计算万有引力，否则分母的R方为0
                dv_temp = gravity_const * m[j] * dt / (pow(x[j]-x[i],2) + pow(y[j]-y[i],2)); //向量 dv，还需要分解到x和y两个方向的dv
                dvx = dvx + dv_temp * ((x[j]-x[i]) / sqrt(pow(x[j]-x[i],2) + pow(y[j]-y[i],2))); 
                dvy = dvy + dv_temp * ((y[j]-y[i]) / sqrt(pow(x[j]-x[i],2) + pow(y[j]-y[i],2)));
            }           
        }
        // 注意计算的时候要使用迭代前的参数值，所以速度迭代和位置迭代才要分成两步进行。先计算速度那么使用老的位置值
        vx[i] = vx[i] + dvx;
        vy[i] = vy[i] + dvy;
    }

    //consider collision in "update_position" function

}


void master() {
    double* m = new double[n_body];
    double* x = new double[n_body];
    double* y = new double[n_body];
    double* vx = new double[n_body];
    double* vy = new double[n_body];

    generate_data(m, x, y, vx, vy, n_body);

    Logger l = Logger("sequential", n_body, bound_x, bound_y);

    for (int i = 0; i < n_iteration; i++){
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        update_velocity(m, x, y, vx, vy, n_body);
        update_position(x, y, vx, vy, n_body);

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

    delete[] m;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    
}


int main(int argc, char *argv[]){
    
    n_body = atoi(argv[1]);
    n_iteration = atoi(argv[2]);

    #ifdef GUI
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(500, 500);
    glutCreateWindow("N Body Simulation Sequential Implementation");
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    gluOrtho2D(0, bound_x, 0, bound_y);
    #endif
    master();

    printf("\nStudent ID: 117010332\n"); // replace it with your student id
    printf("Name: XU Jiale\n"); // replace it with your name
    printf("Total time cost %.3f seconds\n", time_total);
    printf("Assignment 2: N Body Simulation Sequential Implementation\n");
    
    return 0;

}


