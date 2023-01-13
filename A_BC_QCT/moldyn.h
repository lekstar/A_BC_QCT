#ifndef MOLDYN_H
#define MOLDYN_H

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <time.h>

//#include <omp.h>

#include "myrand.h"


#define eVtoE 	9.64853326e-3	// from 1 eV to aem*A^2/(fs^2)

int H_eq_A_BC(double* Q, double* P, double* k, double* K, double* V);

int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps);
int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps, int* react);
int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps, int** react);

int calc_Ev(double* data, int dim, int nstep, double* Ev, double* R);
int calc_vJ(int n, double* data, int dim, int nstep, int* vJ);

int create_init_coord(double* Q, double* P, double* x, double* Rmp0, double b, double Tv, double Vr, int J);

int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps);
int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps, int* react);
int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps, int** react);

int traject_A_BC_print(double* Q_start, double* P_start, double step, double t_fin, double eps);
int traject_A_BC_print(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps);

int find_Rmp_BC(int J, double E, double* Rmp);

double find_Tv_BC(double Rp, int J);
double find_Tv_BC(int N, double Rp, int J);
double find_Tv_BC_par(int N, double Rp, int J, int num_threads);

double find_b_max_A_BC(int N, double* Rmp, int J, double b_start, double Vr, double Tv);
double find_b_max_A_BC_par(int N, double* Rmp, int J, double b_start, double Vr, double Tv, int num_threads);

int P_r_b_fix(int N, double* Rmp, int J, double b, double Vr, double Tv);
int P_r_b_fix(int N, double* Rmp, int J, double b, double Vr, double Tv, int* react);
int P_r_b_fix_par(int N, double* Rmp, int J, double b, double Vr, double Tv, int* react, int num_threads);

int P_r_b_max(int N, double* Rmp, int J, double b_max, double Vr, double Tv);
int P_r_b_max(int N, double* Rmp, int J, double b_max, double Vr, double Tv, int** react);
int P_r_b_max_par(int N, double* Rmp, int J, double b_max, double Vr, double Tv, int num_threads);
int P_r_b_max_par(int N, double* Rmp, int J, double b_max, double Vr, double Tv, int** react, int num_threads);

int k_A_BC_QCT(int Ntr, int v1, int J1, double Ej, double* Ev, int dim_mesh, int num_threads);
int k_A_BC_QCT(int Ntr, int v1, int J1, double Ej, double* Vr, double* bmax1, int dim_mesh, int num_threads);
int k_A_BC_QCT(int Ntr, int v1, double* Ej, int dimJ, double* Etran, int dim_mesh, int num_threads);

#endif