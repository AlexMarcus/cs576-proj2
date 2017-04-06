#include <iostream>
#include <math.h>

using namespace std;

struct point{
	float x, y;
};

point trilateration(point a, point b, point c, float ra, float rb, float rc);

int main(){
	point a = {-2, 4};
	point b = {3, 1};
	point c = {-1, 2};
	point answer = trilateration(a, b, c, 4.472, 3.162,2.236);
	cout << answer.x << endl;
	cout << answer.y;
return 0;
}

point trilateration(point a, point b, point c, float ra, float rb, float rc){
	   float xa = a.x;
	   float ya = a.y;
	   float xb = b.x;
	   float yb = b.y;
	   float xc = c.x;
	   float yc = c.y;

	   	float S = (pow(xc, 2) - pow(xb, 2) + pow(yc, 2) - pow(yb, 2) + pow(rb, 2) - pow(rc, 2)) / 2;
		float T = (pow(xa, 2) - pow(xb, 2) + pow(ya, 2) - pow(yb, 2) + pow(ra, 2) - pow(rc, 2)) / 2;
		float y = ((T * (xb - xc)) - (S * (xb - xa))) / (((ya, yb) * (xb - xc)) - ((yc - yb) * (xb - xa)));
		float x = ((y * (ya)) - T) / (xb - xa);
		point ret;
		ret.x = x;
		ret.y = y;
		return ret;
}
