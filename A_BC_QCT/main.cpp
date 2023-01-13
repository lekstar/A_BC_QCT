#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "myrand.h"

#define DIM_RAND	1000

int main() {

	

	int i;
	double x[DIM_RAND], y[DIM_RAND], x_mean = 0., y_mean = 0., x2 = 0., y2 = 0., xy = 0.;

	mysrand();

	for (i = 0; i < DIM_RAND; i++) {
		x[i] = myrand();
		y[i] = myrand();

		x_mean += x[i];
		y_mean += y[i];
	}

	x_mean /= DIM_RAND;
	y_mean /= DIM_RAND;

	for (i = 0; i < DIM_RAND; i++) {
		x2 += (x[i] - x_mean) * (x[i] - x_mean);
		y2 += (y[i] - y_mean) * (y[i] - y_mean);
		xy += (x[i] - x_mean) * (y[i] - y_mean);
	}

	printf("cov(x;y)=%.10f\n", xy/(sqrt(x2*y2)));

	return 0;
}