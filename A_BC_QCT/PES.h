#ifndef PES_H
#define PES_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include "pnnet.h"

#ifndef sqr(x)
	#define sqr(x)	((x)*(x))
#endif

#define N_ATOM	3
#define N_COORDQ	3
#define N_COORDQ2	9

struct rovib_levels {
	double** Evj;
	int v_max;
	int* Jimax;
};

typedef struct rovib_levels rovib_levels;

extern const int N_neuron1, n_link1[], N_neuron2, n_link2[];

extern const double weights1[], bias1[], weights2[], bias2[];

extern double mass[3], Ze[3];

extern char name_particle[3][15];

extern double D_lim;

extern rovib_levels EvJ_AB;

extern const double Ev_AB[];
extern const int dim_Ev_AB;
extern const double Ev_BC[];
extern const int dim_Ev_BC;
extern const double Ev_AC[];
extern const int dim_Ev_AC;
extern double Te_AB;
extern double Te_AC;
extern double Te_BC;
extern double re_AB;
extern double re_AC;
extern double re_BC;

extern NNet net1, net2;
extern NNet *net;
extern int Num_net;

int init_U();
int init_U(char* pes_file);
int deinit_U();

double U_dist(double* r);
double U_grad_dist(double* r, double* grad);
double U_grad_hess_dist(double* r, double* grad, double* hess);

double U_cart(double* r);
double U_grad_cart(double* r, double* grad);
double U_grad_hess_cart(double* r, double* grad, double* hess);

int U_properties();

int bound_states(double x1, double xN, int N, int J, int order, int r_type, int* N_levels);
int calc_rovib_levels();

int IJ(int n, int k, int *i, int *j);
int IJK(int n, int k, int i);

#endif
