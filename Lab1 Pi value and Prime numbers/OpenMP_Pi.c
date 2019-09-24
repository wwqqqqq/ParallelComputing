#include <stdio.h>
#include <time.h>
#include <omp.h>

#define TOTAL_TIMES 1000

const int num[]={1000,10000,50000,100000};

int main(void) {
    clock_t start, end, total_time;
    double h, pi, sum;
    int i,n;
    int threadnum;
    for(threadnum=1;threadnum<=8;threadnum*=2)
    {
        omp_set_num_threads(threadnum);
        for(i=0;i<4;i++) {
            n=num[i];
            h=1.0/n;
            int times;
            total_time=0;
            for(times=0;times<TOTAL_TIMES;times++)
            {
                
                start=clock();
                sum=0.0;
                #pragma omp parallel for reduction(+:sum)
                for(int j=0;j<n;j++) {
                    double x = h*(j+0.5);
                    sum += 4.0/(1.0+x*x);
                }
                pi = h*sum;
                end=clock();
                
                total_time+=end-start;
            }
            printf("threadnum = %d\tn = %d\ttime = %lfms\npi = %.10lf\n\n",threadnum,n,(double)total_time/CLOCKS_PER_SEC*1000/TOTAL_TIMES,pi);
        }
    }
    return 0;
}