#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#define FAREST_DISTANCE 80000

typedef struct car {
    int velocity;
    struct car* next;
} Car;

typedef struct {
    Car *head;
    int count;
} CarList;

#define MAXV 40
#define P 0.25

const int car_num[] = {100000, 500000, 1000000};
const int cycles[] = {2000, 500, 300}; 

CarList road[FAREST_DISTANCE];
int have_car[FAREST_DISTANCE];
int temp[FAREST_DISTANCE];

void calculateNextV(int myid) {
    srand((unsigned int)time(NULL));
    for(int i=0;i<myid*500;i++)
        rand();
    int j;
    for(j=FAREST_DISTANCE-1;j>=0;j--)
        if(temp[j]>0) break;
    for(int i=j;i>=0;i--) {
        if(road[i].count==0) {
            if(temp[i]>0) j=i;
            continue;
        }
        Car *head = road[i].head;
        Car *p;
        if(i==j) {
            //no cars in the front.
            for(p=head;p!=NULL;p=p->next) {
                if(p->velocity < MAXV) p->velocity++;
                if(rand() <= P * RAND_MAX)  {
                    p->velocity--;
                }
            }
        }
        else {
            int distance = j-i;
            for(p=head;p!=NULL;p=p->next) {
                if(distance<=p->velocity) {
                    p->velocity = (distance-1<0)?0:(distance-1);
                }
                else if(p->velocity < MAXV) p->velocity++;
                if(p->velocity>0 && rand() <= P * RAND_MAX)  {
                    p->velocity--;
                }
            }
        }
        j = i;
    }
}

void moveCars() {
    for(int i=FAREST_DISTANCE-1;i>=0;i--) {
        have_car[i] = 0;
        if(road[i].count==0) {
            continue;
        }
        Car *p=road[i].head;
        Car *q=p;
        while(p!=NULL) {
            if(p->velocity==0) {
                //not moving
                q = p;
                p = p->next;
                have_car[i]++;
            }
            else {
                road[i].count--;
                if(p==road[i].head) {
                    road[i].head = p->next;
                }
                else {
                    q->next = p->next;
                }
                //insert p into target list.
                int new_loc = i+p->velocity;
                if(new_loc>=FAREST_DISTANCE) {
                    printf("out of range!\n");
                    exit(1);
                }
                else {
                    Car* next_p = p->next;
                    if(road[new_loc].head!=NULL)
                    {    
                        Car* temp=road[new_loc].head->next;
                        road[new_loc].head->next = p;
                        p->next = temp;
                    }
                    else {
                        road[new_loc].head = p;
                        p->next = NULL;
                    }
                    p = next_p;
                    road[new_loc].count++;
                    have_car[new_loc]++;
                }
            }
        }
    }

}

int main(int argc, char **argv) {
    double st, ed;
    int size, myid;
    FILE *out_file;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for(int i=0;i<3;i++) {
        //free loc
        for(int loc=0;loc<FAREST_DISTANCE;loc++) {
            if(road[loc].count>0) {
                road[loc].count = 0;
                Car *head = road[loc].head;
                while(head!=NULL) {
                    Car *temp=head->next;
                    free(head);
                    head = temp;
                }
            }
            have_car[loc] = 0;
            road[loc].head = NULL;
        }
        have_car[0] = car_num[i];
        int CarNum = car_num[i]/size;
        road[0].count = CarNum;
        Car *tail;
        for(int j=0;j<CarNum;j++) {
            Car *temp=(Car*)malloc(sizeof(Car));
            temp->next = NULL;
            temp->velocity = 0;
            if(j==0) {
                tail = temp;
                road[0].head = temp;
            }
            else {
                tail->next = temp;
                tail = temp;
            }
        }

        if(myid==0) {
            st = MPI_Wtime();
            printf("car num = %d\tcycles = %d\n",car_num[i],cycles[i]);
        }
        for(int j=0;j<cycles[i];j++) {
            calculateNextV(myid);
            moveCars();
            MPI_Reduce(have_car, temp, FAREST_DISTANCE, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Bcast(temp, FAREST_DISTANCE, MPI_INT, 0, MPI_COMM_WORLD);
            int count=0,count2=0;
            for(int loc=0;loc<FAREST_DISTANCE;loc++) {
                count+=temp[loc];
                count2+=have_car[loc];
            }
        }
        if(myid==0) {    
            ed = MPI_Wtime();
            printf("time = %lfms\n",(double)(ed-st)*1000);
            char filename[]="xxx.txt";
            sprintf(filename, "%d_%d.txt",i,size);
            out_file = fopen(filename,"w");
            if(out_file==NULL) {
                printf("open file fail\n");
                exit(2);
            }
            int count=0;
            for(int loc=0;loc<FAREST_DISTANCE;loc++) {
                if(temp[loc]>0) {
                    fprintf(out_file,"%d %d\n",loc,temp[loc]);
                    count+=temp[loc];
                }
            }
            printf("\n************\ncount=%d\n************\n\n",count);
            fclose(out_file);
        } 

    }
    return 0;
}