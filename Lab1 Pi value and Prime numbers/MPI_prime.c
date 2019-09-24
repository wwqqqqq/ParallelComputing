#include <stdio.h>
#include "mpi.h"
#include <time.h>
#include <math.h>
#define TOTAL_TIMES 10

const int num[]={1000,10000,100000,500000};

int isPrime(int n) {
    int i;
    int s = sqrt((double)n);
    for(i=2;i<=s;i++) {
        if(n%i==0) return 0;
    }
    return 1;
}

int main(int argc, char **argv) {
    double start, end,total;
    int size, myid;
    int i,n,count;
    int sum;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    printf("size = %d\n",size);
    for(i=0;i<4;i++)
    {
        count=0;
        sum=0;
        n = num[i];
        total = 0.0;
        for(int times=0;times<TOTAL_TIMES;times++)
        {
            if(myid == 0)
                start = MPI_Wtime();
            for(int j=myid;j<=n;j+=size) {
                count += isPrime(j);
            }
            MPI_Reduce(&count, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            if(myid == 0) {
                end = MPI_Wtime();
                total += end-start;
            }
        }
        if(myid == 0) {
            //end = MPI_Wtime();
            printf("n = %d\ttime = %.10lf\ncount = %d\n\n",n,(total)*1000/TOTAL_TIMES,sum);
        }
    }
    MPI_Finalize();
    return 0;
}