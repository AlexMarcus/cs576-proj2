#include <stdio.h>
#include <iostream>
#include <cuda_runtime.h>
#include <cuda.h>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <math.h>

#define NUM (1<<12)
#define U 2
#define V 16

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
    //point *results =(point *) malloc((NUM/4) * (sizeof(point)));
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
    for(int i = 0; i < NUM; i++){
	cudaMallocManaged(&dv[i], 3*sizeof(float));
    }

    *da = a;
    *db = b;
    *dc = c;
    for(int i = 0; i < NUM; i++){
	for(int j = 0; j < 3; j++){
		dv[i][j] = distance_vector[i][j];
	}
    }

    trilateration<<<U,V>>>(da,db,dc,dv,pts);
    cudaDeviceSynchronize();
    
    //cudaMemcpy(results, pts, (NUM/4) * sizeof(point),cudaMemcpyDeviceToHost);

    /*for(int i = 0; i < NUM/4; i++){
	if(results[i].x != 0)
		cout << results[i].x << ", " << results[i].y << "\n";
    }*/

    for(int i = 0; i < 20; i++){
		cout << pts[i].x << ", " << pts[i].y << "\n";
    }

	/*
	first points
	0.170442, -0.212715
	0.642852, -0.825177
	1.1408, -1.19354
	1.42159, -1.6076
	1.62658, -2.18829
	1.77826, -2.66518
	2.17155, -2.97643
	2.65642, -3.23502
	3.06065, -3.62799
	3.42531, -4.1356
	3.96662, -4.42709
	4.41263, -4.69845
	4.93671, -4.96385
	5.46347, -5.28732
	6.06384, -5.70538
	6.48004, -6.27095
	7.05871, -6.8165
	7.63242, -7.21735
	8.05061, -7.52002
	8.4322, -7.82146
	8.73583, -8.10203
	9.03263, -8.45821
	9.40358, -8.86348
	9.6679, -9.13161
 

	*/


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


	  //cout << dist[i][0] << "," << dist[i][1] << "," << dist[i][2] << endl;
	  
	  if(i%4 == 0 && i != 0){
	  	point t;
		t.x = x_ave/4;
		t.y = y_ave/4;
		points[(i/4)-1] = t;
		x_ave = 0;
		y_ave = 0;
		//if(i < 100) cout << t.x << ", " << t.y << endl;
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

	   int i = blockIdx.x * blockDim.x + threadIdx.x;
		   
	   int j,k;
	   for(j =0; j < ((NUM/4)/(U*V));j++){
	   	 float ave_y = 0, ave_x = 0;
	   	 for(k = 0; k < 4; k++){
	   	        float xa = a->x;
		 	float ya = a->y;
	   	 	float xb = b->x;
	   	 	float yb = b->y;
	   	 	float xc = c->x;
	   	 	float yc = c->y;
	   	 	float ra = dv[i+ j*U*V + i][0];
	   	 	float rb = dv[i+ j*U*V][1];
	   	 	float rc = dv[i+ j*U*V][2];

			float A = -2*xa + (2*xb);
	   		float B = -2*ya + (2*yb);
	   		float C = pow(ra, 2) - pow(rb,2) - pow(xa, 2) + pow(xb, 2) - pow(ya, 2) + pow(yb, 2);
	   		float D = -2*xb + (2*xc);
   	   		float E = -2*yb + (2*yc);
	   		float F = pow(rb, 2) - pow(rc, 2) - pow(xb, 2) + pow(xc, 2) - pow(yb, 2) + pow(yc, 2);

			ave_x += (C*D - F*A) / (B*D - E*A);
			ave_y += (A*E - D*B) / (C*E - F*B);
		}
	   	point ret;
	   	ret.x = ave_x/4;
	   	ret.y = ave_y/4;
	   	pts[i + j*U*V] = ret;
	   	syncthreads();
	}
}
