#include <stdio.h>
#include <time.h>
#include <omp.h>
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

int main(void) {
    int i,n,count;
    clock_t start,end,total;
    int threadnum;
    for(threadnum=1;threadnum<=8;threadnum*=2)
    {
        
        for(i=0;i<4;i++) {
            n=num[i];
            count=0;
            total = 0;
            for(int times=0;times<TOTAL_TIMES;times++)
            {
                start=clock();

                omp_set_num_threads(threadnum);
                #pragma omp parallel for reduction(+:count)
                for(int j=2;j<=n;j++) {
                    count+=isPrime(j);
                }
                
                end=clock();
                total += end-start;
            }
            printf("threadnum = %d\tn = %d\ttime = %lfms\ncount = %d\n\n",threadnum,n,(double)(total)/CLOCKS_PER_SEC*1000/TOTAL_TIMES,count);
        }
    }
    return 0;
}