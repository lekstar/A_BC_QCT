#ifndef TISE1D_SOL_H
#define TISE1D_SOL_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <mkl.h>

int TISE1D_sol(double* V, double x1, double xN, int N, int order, double Emin, int iElim);

#endif
