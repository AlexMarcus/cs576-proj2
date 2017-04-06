#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <math.h>

#define NUM (1<<12)

using namespace std;

typedef struct{
        float x,y;
} point;

float ** generate_set(point a, point b, point c, point *points);
float get_distance(point a, point b);
__global__ void trilateration(point *a, point *b, point *c, float ** dv, point * pts);



int main(int argc, char *argv[]){
    srand(time(NULL));
    cout << NUM << endl;
    point *results =(point *) malloc((NUM/4) * (sizeof(point)));
    point *points =(point *) malloc((NUM/4) * (sizeof(point)));

    point a = {3.4,-2.4};
    point b = {5.6,1.23};
    point c = {-3.8,5.4};
    
    float ** distance_vector = generate_set(a,b,c,points);

    float ** dv;
    point * da;
    point * db;
    point * dc;
    point * pts;

    /*cudaMalloc(&da, sizeof(point *));
    cudaMalloc(&db, sizeof(point *));
    cudaMalloc(&dc, sizeof(point *));
    cudaMalloc((void **)&pts, (NUM/4) * sizeof(point));
    cudaMalloc((void **)&dv, NUM*sizeof(float *));

    for(int i = 0; i < NUM; i++){
    	    cudaMalloc(&dv[i], 3*sizeof(float));
    }

    cudaMemcpy(dv, distance_vector, NUM * sizeof(float*), cudaMemcpyHostToDevice);
    cudaMemcpy(da, &a, sizeof(point),cudaMemcpyHostToDevice);
    cudaMemcpy(db, &b, sizeof(point),cudaMemcpyHostToDevice);
    cudaMemcpy(dc, &c, sizeof(point),cudaMemcpyHostToDevice);*/

    cudaMallocManaged(&da, sizeof(point *));
    cudaMallocManaged(&db, sizeof(point *));
    cudaMallocManaged(&dc, sizeof(point *));
    cudaMallocManaged(&pts, (NUM/4) * sizeof(point));
    cudaMallocManaged(&dv, NUM * sizeof(float *));
    /*for(int i = 0; i < NUM; i++){
	cudaMallocManaged(&dv[i], 3*sizeof(float));
    }

    *da = a;
    *db = b;
    *dc = c;
    for(int i = 0; i < NUM; i++){
	for(int j = 0; j < 3; j++){
		dv[i][j] = distance_vector[i][j];
	}
    }*/

    trilateration<<<1,1>>>(da,db,dc,dv,pts);
    cudaDeviceSynchronize();
    
    //cudaMemcpy(results, pts, (NUM/4) * sizeof(point),cudaMemcpyDeviceToHost);

    /*for(int i = 0; i < NUM/4; i++){
	if(results[i].x != 0)
		cout << results[i].x << ", " << results[i].y << "\n";
    }*/

    for(int i = 0; i < NUM/4; i++){
	if(pts[i].x == 32)
		cout << pts[i].x << ", " << pts[i].y << "\n";
    }
    return 0;    
}

float ** generate_set(point a, point b, point c, point *points){

    float ** dist = (float **) malloc(NUM * sizeof(float *)); 
    int i,j;
    for(j = 0; j < NUM; j++){
    	  dist[j] = (float *) malloc(3 * sizeof(float));
	  for(i = 0; i < 3; i++){
	  	dist[j][i] = 0;
	}	     
    }

    srand(time(NULL));
    float x_ave = 0, y_ave = 0;
    point next;
    next.x = 0;
    next.y = 0;
    for(i = 0; i < NUM; i++){
    	  dist[i][0] = get_distance(a,next);
	  dist[i][1] = get_distance(b,next);
	  dist[i][2] = get_distance(c,next);


	  cout << dist[i][0] << "," << dist[i][1] << "," << dist[i][2] << endl;
	  
	  if(i%4 == 0 && i != 0){
	  	point t;
		t.x = x_ave/4;
		t.y = y_ave/4;
		points[(i/4)-1] = t;
		x_ave = 0;
		y_ave = 0;
	  }	  

	  x_ave	+= next.x;
	  y_ave	+= next.y;

	  //get new point
	  float temp = (rand() % 20000);
	  float delta_x = (temp / 100000);
	  temp = (rand() % 20000);
	  float delta_y = (temp / 100000) - .2;
	  
	  next.x += delta_x;
	  next.y += delta_y;
	  
    }

    return dist; 
    	  
}

float get_distance(point a, point b){
      float distance = sqrt((pow((a.x - b.x),2) + pow((a.y - b.y),2)));
      return distance;
}

float norm(point p){
	return pow(pow(p.x,2) + pow(p.y,2), .5);
}

__global__ void trilateration(point *a, point *b, point *c, float ** dv, point * pts){

	   int i = threadIdx.x;
	   pts[i].x = 32;
	   pts[i].y = 123;
	   /*float xa = a->x;
	   float ya = a->y;
	   float xb = b->x;
	   float yb = b->y;
	   float xc = c->x;
	   float yc = c->y;
	   float ra = dv[i][0];
	   float rb = dv[i][1];
	   float rc = dv[i][2];

	   	float S = (pow(xc, 2) - pow(xb, 2) + pow(yc, 2) - pow(yb, 2) + pow(rb, 2) - pow(rc, 2)) / 2;
		float T = (pow(xa, 2) - pow(xb, 2) + pow(ya, 2) - pow(yb, 2) + pow(ra, 2) - pow(rc, 2)) / 2;
		float y = ((T * (xb - xc)) - (S * (xb - xa))) / (((ya, yb) * (xb - xc)) - ((yc - yb) * (xb - xa)));
		float x = ((y * (ya)) - T) / (xb - xa);
		point ret;
		ret.x = x;
		ret.y = y;
		pts[i] = ret;*/

}
