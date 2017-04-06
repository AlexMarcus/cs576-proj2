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
	point answer = trilateration(a, b, c, 4.24, 2, 2.236);
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

		float A = -2*xa + (-2*xb);
		float B = -2*ya + (-2*yb);
		float C = pow(ra, 2) - pow(rb, 2) - pow(xa, 2) + pow(xb, 2) - pow(ya, 2) + pow(yb, 2);
		float D = -2*xb + (-2*xc);
		float E = -2*yb + (-2*yc);
		float F = pow(rb, 2) - pow(rc, 2) - pow(xb, 2) + pow(xc, 2) - pow(yb, 2) + pow(yc, 2);
	   	/*float S = (pow(xc, 2) - pow(xb, 2) + pow(yc, 2) - pow(yb, 2) + pow(rb, 2) - pow(rc, 2)) / 2;
		float T = (pow(xa, 2) - pow(xb, 2) + pow(ya, 2) - pow(yb, 2) + pow(ra, 2) - pow(rc, 2)) / 2;
		float y = ((T * (xb - xc)) - (S * (xb - xa))) / (((ya, yb) * (xb - xc)) - ((yc - yb) * (xb - xa)));
		float x = ((y * (ya)) - T) / (xb - xa);*/
		point ret;
		ret.x = (C*D - F*A) / (B*D - E*A);
		ret.y = (A*E - D*B)/(C*E - F*B);
		return ret;
}
