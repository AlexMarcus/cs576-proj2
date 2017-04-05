#include <iostream>
#include <cuda_runtime.h>
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

int main(int argc, char *argv[]){
    srand(time(NULL));
    cout << NUM << endl;

    point *points =(point *) malloc((NUM/4) * (sizeof(point)));

    point a = {3.4,-2.4};
    point b = {5.6,1.23};
    point c = {-3.8,5.4};
    
    float ** distance_vector = generate_set(a,b,c,points);

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
	  float delta_x = .1 - (temp / 100000);
	  temp = (rand() % 20000);
	  float delta_y = .1 - (temp / 100000);
	  
	  next.x += delta_x;
	  next.y += delta_y;
	  
    }

    return dist; 
    	  
}

float get_distance(point a, point b){
      float distance = sqrt((pow((a.x - b.x),2) + pow((a.y - b.y),2)));
      return distance;
}     