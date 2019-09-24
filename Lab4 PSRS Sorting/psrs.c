#include <stdio.h>
#include "mpi.h"
#include <time.h>
#include <stdlib.h>
#define MAXKEY RAND_MAX

const int scale[] = {64, 1000, 10000, 30000};

/*
** phase 1: initialization
set up p processors, let the root processor, 0, get data of size n.
** phase 2: scatter data, local sort and regular samples collected
scatter the data values to the p processors. Each processors sorts its local data set, roughtly of size n/p, using quicksort.
each processor chooses p sample points, in a very regular manner, from its locally sorted data.
** phase 3: gather and merge samples, choose and broadcast p-1 pivots
the root processor, 0, gathers the p sets of p sampel points. 
it is important to realize that each set of these p sample points is sorted. These p sets are sorted using multimerge.
from these p^2 sorted points, p-1 pivot values are regularly chosen and are broadcast to the other p-1 processors.
** phase 4: local data is partitioned
each of the p processors partitions its local sorted data, roughly of size n/p, into p classes using the p-1 pivot values.
** phase 5: all *ith* classes are gathered and merged
processor i gathers the ith class of data from every other processor.
each of these classes is sorted using multimerge.
** phase 6: root processor collects all the data
the root processor gathers all the data and assembles the sorted list of n values.
*/

int data[5000000];
int s_data[5000000];
int recv_data[7][5000000];
int sample_point[8];
int points[64];
int t[8][625000];
int pivot[8];

void init(int N) {
    srand((unsigned int)time(NULL));
    for(int i=0;i<N;i++)
        data[i] = rand()%10000;
}

int my_cmp(const void* ele1, const void* ele2) {
    return *((int*)ele1) - *((int*)ele2);
}

void choose_sample_point(int N, int p, int myid) {
    int w = N/p/p;
    if(myid == 0) {
        for(int i=0;i<p;i++)
            t[0][i] = s_data[i*w];
        t[0][p] = MAXKEY;
    }
    for(int i=0;i<p;i++)
        sample_point[i] = s_data[i*w];
}

int merge(int* result, int len) {
    int i=0,j=0;
    int count=0;
    while(count<len*2 && i<len && j<len) {
        if(t[0][i]<t[1][j]) {
            result[count] = t[0][i];
            i++;
        }
        else if(t[0][i]>t[1][j] || t[0][i]!=MAXKEY) {
            result[count] = t[1][j];
            j++;
        }
        else
            break;
        count++;
    }
    result[count] = MAXKEY;
    return count;
}

int MultiMerge(int p, int* result, int len) {
    if(p == 2) return merge(result,len);
    int val[p*2];
    int pos1[p*2];
    int pos2[p*2];
    int count = 0;
    for(int i=0;i<p*2;i++) val[i] = -1;
    for(int i=0;i<p;i++) {
        int p1 = i;
        int p2 = 0;
        int ind = p1+p;
        int temp = t[p1][0];
        val[ind] = t[p1][0];
        pos1[ind] = p1;
        pos2[ind] = 0;
        ind /= 2;
        while(ind>=0) {
            if(temp>val[ind]) {
                int last_t = temp;
                int last_p1 = p1;
                int last_p2 = p2;
                temp = val[ind];
                p1 = pos1[ind];
                val[ind] = last_t;
                pos1[ind] = last_p1;
                pos2[ind] = last_p2;
            }
            if(ind==0) break;
            ind = ind/2;
        }
    }
    while(1) {
        if(val[0]==MAXKEY) break;
        result[count] = val[0];
        count++;
        if(count==p*p) break;
        int p1 = pos1[0];
        int p2 = pos2[0];
        int ind = p + p1;
        p2++;
        if(p1>p || p2>len) {
            val[ind] = MAXKEY;
        }
        else val[ind] = t[p1][p2];
        pos2[ind] = p2;
        int temp = val[ind];
        ind /= 2;
        while(ind>=0) {
            if(temp>val[ind]) {
                int last_t = temp;
                int last_p1 = p1;
                int last_p2 = p2;
                temp = val[ind];
                p1 = pos1[ind];
                p2 = pos2[ind];
                val[ind] = last_t;
                pos1[ind] = last_p1;
                pos2[ind] = last_p2;
            }
            if(ind==0) break;
            ind = ind/2;
        }
    }

    result[count] = MAXKEY;
    return count;
}

void choose_pivot(int p) {
    for(int i=0;i<p-1;i++)
        pivot[i] = points[(i+1)*p];
}

int class[8][625000];

void partition(int N, int p) {
    int cls = 0;
    int count = 0;
    for(int i=0;i<N;) {
        if(s_data[i]<pivot[cls]) {
            class[cls][count] = s_data[i];
            count++;
            i++;
        }
        else {
            class[cls][count] = MAXKEY;
            count = 0;
            cls++;
            if(cls==p) {
                for(;i<N;i++)
                    class[p-1][count++] = s_data[i];
                class[p-1][count] = MAXKEY;
                break;
            }
        }
    }
    if(cls<p && count<625000)
        class[cls][count] = MAXKEY;
    for(cls++;cls<p;cls++) class[cls][0] = MAXKEY;
}

int main(int argc, char** argv) {
    double start, end;
    int size, myid;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for(int i=0;i<4;i++) {
        int N = scale[i];
        if(size == 1) {
            start = MPI_Wtime();
            qsort(s_data, N/size, sizeof(int), my_cmp);
            end = MPI_Wtime();
            printf("N = %d\ttime = %.10lf\n\n",N,(end-start)*1000);
            continue;
        }
        //phase 1
        if(myid == 0) {
            printf("begin N = %d\n",N);
            init(N);
            start = MPI_Wtime();
        }
        //phase 2
        MPI_Scatter(data, N/size, MPI_INT, s_data, N/size, MPI_INT, 0, MPI_COMM_WORLD);
        qsort(s_data, N/size, sizeof(int), my_cmp);
        choose_sample_point(N, size, myid);
        //phase 3
        if(myid == 0) {
            for(int id=1;id<size;id++) {
                MPI_Recv(t[id], size, MPI_INT, id, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                t[id][size] = MAXKEY;
            }
            MultiMerge(size, points, size+1);
            choose_pivot(size);
        }
        else {
            MPI_Send(sample_point, size, MPI_INT, 0, myid, MPI_COMM_WORLD);
        }
        MPI_Bcast(pivot, size-1, MPI_INT, 0, MPI_COMM_WORLD);
        pivot[size-1] = MAXKEY;
        //phase 4
        partition(N, size);
        //phase 5
        for(int id=0;id<size;id++) {
            if(id!=myid) {
                MPI_Send(class[id],N/size,MPI_INT,id,(myid+1)*size+id,MPI_COMM_WORLD);
            }
        }
        for(int id=0;id<size;id++) {
            if(id==myid) continue;
            MPI_Recv(t[id],N/size,MPI_INT,id,(id+1)*size+myid,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        
        int j;
        for(j=0;j<N/size && class[myid][j]<MAXKEY;j++)
            t[myid][j] = class[myid][j];
        if(j<N/size)
            t[myid][j] = MAXKEY;
        int len = MultiMerge(size, s_data, N/size);
        //phase 6
        if(myid == 0) {
            int count = 0;
            for(j=0;count<N && s_data[j]<MAXKEY;count++,j++) {
                data[count] = s_data[j];
            }
            for(int id=1;id<size;id++) {
                MPI_Recv(recv_data[id-1],N,MPI_INT,id,id+size*size+size,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(j=0;count<N && j<N && recv_data[id-1][j]<MAXKEY;j++) {
                    data[count++] = recv_data[id-1][j];
                }
            }
            end = MPI_Wtime();
            printf("N = %d\ttime = %.10lf\n\n",N,(end-start)*1000);

        }
        else {
            MPI_Send(s_data,N,MPI_INT,0,myid+size*size+size,MPI_COMM_WORLD);
        }

    }
    return 0;
}

