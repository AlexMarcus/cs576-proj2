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
	point c = {-1, -2};
	point answer = trilateration(a, b, c, 7.8102, 3, 4);
	cout << fixed <<answer.x << endl;
	cout <<fixed <<  answer.y << endl;
	return 0;
}

point trilateration(point a, point b, point c, float ra, float rb, float rc){
	   float xa = a.x;
	   float ya = a.y;
	   float xb = b.x;
	   float yb = b.y;
	   float xc = c.x;
	   float yc = c.y;

	   float numerator = ((xb - xa) * (xc * xc + yc * yc - rc*rc) +
				(xa - xc) * (xb * xb + yb * yb - rb * rb) +
				(xc - xb) * (xa * xa + ya * ya - ra * ra));
	   float denominator = (2 * (yc *(xb - xa) + yb * (xa - xc) + ya * (xc - xb)));
	   float y = numerator/denominator;
	   float x = (rb * rb + xa * xa + ya * ya - ra * ra - xb * xb - yb * yb - 2*(ya - yb) * y) / (2*(xa -xb));
		point ret;
		ret.x = x;
		ret.y = y;
		return ret;
}
