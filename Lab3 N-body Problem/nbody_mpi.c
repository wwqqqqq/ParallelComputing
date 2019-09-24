#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define G (6.67e-11) //N*m^2/kg^2
#define RADIUS 1 //cm
#define WEIGHT 10000 //kg
#define TIMESTEP 1000
#define TIME 60 //s
#define DELTA 1e-3

const int num[] = {64, 256};

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
void compute_force(int N, int myid, int size, double* temp_x, double* temp_y);
void compute_velocities(int beg, int end);
void compute_positions(int beg, int end);
void print_positions(int N);

int main(int argc, char** argv) {
    int size, myid;
    double st, ed;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for(int i=0;i<2;i++) {
        int N = num[i];
        int beg = myid*N/size;
        int end = (myid+1)*N/size-1;
        Init(N);
        if(myid==0) {
            st = MPI_Wtime();
            //print_positions(N);
        }
        for(int j=0;j<TIME*TIMESTEP;j++) {
            double temp_x[256];
            double temp_y[256];
            compute_force(N, myid, size, temp_x, temp_y);
            MPI_Reduce(temp_x, F_x, 256, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(temp_y, F_y, 256, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(F_x, 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(F_y, 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            compute_velocities(beg, end);
            compute_positions(beg, end);
            if(myid == 0) {
                if(size > 1)
                for(int id=1;id<size;id++) {
                    int tag1 = id+2*size*(j+i*TIMESTEP*TIME);
                    int tag2 = tag1+size;
                    MPI_Recv(loc_x+id*N/size, N/size, MPI_DOUBLE, id, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(loc_y+id*N/size, N/size, MPI_DOUBLE, id, tag2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            else {
                int tag1 = myid+2*size*(j+i*TIMESTEP*TIME);
                int tag2 = tag1+size;
                MPI_Send(loc_x+myid*N/size, N/size, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD);
                MPI_Send(loc_y+myid*N/size, N/size, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);
            }
            MPI_Bcast(loc_x, 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(loc_y, 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
        if(myid == 0){
            ed = MPI_Wtime();
            //print_positions(N);
            printf("N = %d\n",N);
            printf("(%lf, %lf)\n",loc_x[N-1],loc_y[N-1]);
            printf("time = %lfms\n",(ed-st)*1000);
        }
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

void compute_force(int N, int myid, int size, double* temp_x, double* temp_y) {
    for(int i=0;i<N;i++) {
        temp_x[i] = 0;
        temp_y[i] = 0;
    }
    for(int i=myid;i<N;i+=size) {
        for(int j=i+1;j<N;j++) {
            double distance_x = loc_x[i]-loc_x[j];
            double distance_y = loc_y[i]-loc_y[j];
            double r = sqrt(distance_x*distance_x + distance_y*distance_y); //cm
            double F = calculate_force(r);
            temp_x[i] -= distance_x / r * F;
            temp_y[i] -= distance_y / r * F;
            temp_x[j] += distance_x / r * F;
            temp_y[j] += distance_y / r * F;
        }
    }
}
// delta_v = delta_t * a
void compute_velocities(int beg, int end) {
    double a_x, a_y;
    for(int i=beg;i<=end;i++) {
        a_x = F_x[i]/WEIGHT; //m*s^-2
        a_y = F_y[i]/WEIGHT; //m*s^-2
        v_x[i] += DELTA * a_x;
        v_y[i] += DELTA * a_y;
    }
}
//d = v * t
void compute_positions(int beg, int end) {
    for(int i=beg;i<=end;i++) {
        loc_x[i] += DELTA * v_x[i];
        loc_y[i] += DELTA * v_y[i];
    }
}

void print_positions(int N) {
    printf("************* N = %d *************\n",N);
    for(int i=0;i<N;i++) {
        printf("(%lf, %lf)\n",loc_x[i],loc_y[i]);
    }
    printf("**************************************\n\n");
}