#ifndef PES_ANALYSIS_H
#define PES_ANALYSIS_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdbool.h>
#include <conio.h>
#include "pnnet.h"
#include "eispack.h"
#include "PES.h"
#include "symmetry.h"

extern double	kb_c, amu_kg, hartree_J, eV_J, h_c, h1_c, c_c, Na, bohr_m, bohr_A;

int freq(double* hess, double* mass);
int locate_min(double* r1, int NSTEP, double OPTTOL);
int locate_CP(double* r1, int NSTEP, double OPTTOL, bool CP);

int freq_analysis(int num_atom, double* charge, double* coord, double* hess, double* eigene, double* eigenh, double* red_m, double B0[3], double* sigma);
int hess_mw_eigensystem(double* r, double* hess, double* eigenv, double* eigenh);
int hess_eigensystem(double* r, double* hess, double* eigenv, double* eigenh);

int MEP_GS2(double* r, double stride, int maxiter, bool TS, bool dir, double opttol);

int k_A_BC_CTST(double* r0, double* T_list, int N_T);

int group_order(char* group);

int e_ijk(int i, int j, int k);


#endif
