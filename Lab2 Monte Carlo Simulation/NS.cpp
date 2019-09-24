#include <stdio.h>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

using namespace std;

#define MAXV 40
#define P 0.25

const int car_num = {100000, 500000, 1000000};
const int cycles = {2000, 500, 300}; 

typedef struct{
    int loc;
    int v;
}car;

vector<car> cars;

int my_cmp(car A, car B) {
    return A.loc-B.loc;
}

void calculateNextV(vector<car>& cars) {
    int n = cars.size();
    srand((unsigned int)time(NULL));
    for(int i=0;i<n;i++) {
        if(i==n-1) {
            cars[i].v++;
        }
        else {
            int j;
            for(j=i+1;j<n;j++)
                if(cars[j].loc != cars[i].loc) break;
            if(cars[j].loc-cars[i].loc<=cars[i].v) {
                cars[i].v = cars[j].loc-cars[i].loc-1;
            }
            else cars[i].v++;
        }
        if(cars[i].v > MAXV) cars[i].v = MAXV;
        if(rand() <= p * RAND_MAX && cars[i].v>0) cars[i].v--; 
    }
    sort(cars.begin(), cars.end(), my_cmp);
}

void moveCars(vector<car>& cars) {
    int n = cars.size();
    for(int i=0;i<n;i++) 
        cars[i].loc += cars[i].v;
}

int main(void) {
    clock_t st, ed;
    for(int i=0;i<3;i++) {
        cars.erase(cars.begin(), cards.end());
        for(int j=0;j<car_num[i];j++) {
            car temp;
            temp.loc = 0; temp.v = 0;
            cars.push_back(temp);
        }
        st = clock();
        for(int j=0;j<cycles;j++) {
            calculateNextV(cars);
            moveCars(cars);
        }
        ed = clock();
        printf("car num = %d\tcycles = %d\ttime = %lfms\n",car_num[i],cycles[i],(double)(ed-st)/CLOCKS_PER_SEC*1000);

    }
}