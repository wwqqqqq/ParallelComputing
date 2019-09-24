#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define G (6.67e-11) //N*m^2/kg^2
#define RADIUS 1 //cm
#define WEIGHT 10000 //kg
#define TIMESTEP 10000
#define TIME 60 //s
#define DELTA 1e-4

const int num = {64, 256};

//delta v = a delta t
//m a  = F
//F = GMm/r^2

double F_x[256];
double F_y[256];
double loc_x[256];
double loc_y[256];
double v_x[256];
double v_y[256];
void Init(int N);
void compute_force(int N);
void compute_velocities(int N);
void compute_positions(int N);

int main(int argc, char**argv) {
    clock_t st, ed;
    for(int i=0;i<2;i++) {
        int N = num[i];
        Init(N);
        st = clock();
        for(int j=0;j<TIME*TIMESTEP;j++) {
            compute_force(N);
            compute_velocities(N);
            compute_positions(N);
        }
        ed = clock();
        PrintPosition(N);
        printf("time = %lfms\n",(double)(ed-st)*1000/CLOCKS_PER_SEC);
    }
}

void Init(int N) {
    int width = (N==64)?8:16;
    for(int i=0;i<N;i++) {
        v_x[i] = 0;
        v_y[i] = 0;
        F_x[i] = 0;
        F_y[i] = 0;
        //小球初始间隔为1cm
        loc_x[i] = i/width;
        loc_y[i] = i%width;
    }
}

double calculate_force(double r) {
    if(r<RADIUS*2) r = RADIUS*2;
    r = r / 100; //m
    return G*WEIGHT*WEIGHT/(r*r); //N
}

void compute_force(int N) {
    for(int i=0;i<N;i++) {
        F_x[i] = 0;
        F_y[i] = 0;
    }
    for(int i=0;i<N;i++) {
        for(int j=i+1;j<n;j++) {
            double distance_x = loc_x[i]-loc_x[j];
            double distance_y = loc_y[i]-loc_y[j];
            double r = sqrt(distance_x*distance_x + distance_y*distance_y); //cm
            double F = calculate_force(r);
            F_x[i] -= distance_x / r * F;
            F_y[i] -= distance_y / r * F;
            F_x[j] += distance_x / r * F;
            F_y[j] += distance_y / r * F;
        }
    }
}
// delta_v = delta_t * a
void compute_velocities(int N) {
    for(int i=0;i<N;i++) {
        double a_x = F_x[i]/WEIGHT; //m*s^-2
        double a_y = F_y[i]/WEIGHT; //m*s^-2
        v_x[i] += DELTA * a_x;
        v_y[i] += DELTA * a_y;
    }
}
//d = v * t
void compute_positions(int N) {
    for(int i=0;i<N;i++) {
        loc_x[i] += DELTA * v_x[i];
        loc_y[i] += DELTA * v_y[i];
    }
}
