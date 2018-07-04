#include <math.h>
#include "image_features.h"

point_2d_64f_s point_2d_create(double x, double y)
{
	point_2d_64f_s p;
	p.x = x;
	p.y = y;
	return p;
}
double descr_dist_sq(struct feature* f1, struct feature* f2)
{
	double diff, dsq = 0;
	double *descr1, *descr2;
	int i, d;
	d = f1->d;
	if (f2->d != d)
	{
		return 1e9;
	}
	descr1 = f1->descr;
	descr2 = f2->descr;
	for (i = 0; i < d; i++)
	{
		diff = descr1[i] - descr2[i];
		dsq += diff * diff;
	}
	return dsq;
}
