#include <stdio.h>
#include "mpi.h"
#include <time.h>

#define TOTAL_TIMES 1000

const int num[]={1000,10000,50000,100000};

//pi=4/(1+x^2)在0到1上的积分

int main(int argc, char **argv) {
    double pi, h, sum, x, start, end, total;
    int size, myid;
    int i;
    int n;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("size = %d\n",size);
    for(i=0;i<4;i++)
    {
        n = num[i];
        total = 0.0;
        for(int times=0;times<TOTAL_TIMES;times++)
        {
            if(myid == 0)
                start = MPI_Wtime();
            h = 1.0/n;
            sum=0.0;
            int j;
            for(j=myid;j<n;j+=size) {
                x = h*(j+0.5);
                sum += 4.0/(1+x*x);
            }
            sum = sum*h;
            MPI_Reduce(&sum, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid==0) {
                end = MPI_Wtime();
                total += end-start;
            }
        }
        if(myid == 0) {
            //end = MPI_Wtime();
            printf("n = %d\ttime = %.10lf\npi = %lf\n\n",n,(double)total*1000/TOTAL_TIMES,pi);
        }
    }
    MPI_Finalize();
    return 0;
}