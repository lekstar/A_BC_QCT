#include "PES_analysis.h"

double	kb_c = 1.380649e-23,				// J/K, Boltsmann constant
		amu_kg = 1.6605390666e-27,			// 1 amu to 1 kg
		hartree_J = 4.3597447222071e-18,	// 1 hartree to 1 J
		eV_J = 1.602176634e-19,				// 1 eV to 1 J
		h_c = 6.62607015e-34,				// J*s, Plank constant
		h1_c = 1.054571817e-34,				// J*s, Dirak constant
		c_c = 299792458,					// m/s, speed of light in vacuum
		Na = 6.02214076e23,					// Avogadro number
		bohr_m = 0.52917721067e-10,			// 1 Bohr to 1 meter
		bohr_A = 0.52917721067;				// 1 Bohr to 1 Angstrem

double determ_matrix( double *p, int n ){
	double buf, max, *a, q, r;
	int i, j, k, l = 0, m, kmax;
	
	a = (double*)malloc(sizeof(double)*n*n);
	
	for(i=0; i<n*n; i++)
		a[i] = p[i];
	
	for(i=0; i<n-1; i++){
		max = a[i*n+i];
		kmax = i;
		for(j=i+1; j<n; j++){
			if(fabs(max)<fabs(a[j*n+i])){
				kmax = j;
				max = a[j*n+i];
			}
			
			if(kmax!=i) l++;
		}
		
		for(j=0; j<n; j++){
			buf = a[i*n+j];
			a[i*n+j] = a[kmax*n+j];
			a[kmax*n+j] = buf;
		}
		
		for(k=i+1; k<n; k++){
			q = a[k*n+i];
			r = a[i*n+i];
			if(r == 0.) return 0.;
			for(j=i; j<n; j++)
				a[k*n+j] -= a[i*n+j]*q/r;	
		}
	}
	
//	for(i=0; i<n; i++)
//		for(j=0; j<n; j++)
//			if(j==n-1) printf("%7.4f\n", a[n*i+j] );
//			else printf("%7.4f\t", a[n*i+j] );	
	
	if(l%2==1)
		buf = -1.;
	else 
		buf = 1.;
	
	for(i=0; i<n; i++)
		buf *= a[i*n+i];
			
	free(a);
	
	return buf;
}

int invert_matrix( double *p, double *inv, int n ){
	double *a, max, buf, q, r;
	int i, j, k, kmax, l = 0;
	
	if (fabs(determ_matrix(p, n))<1e-100){
		printf("\n\nMatrix is singular\n\n");
		return 1;
	}
	
	a = (double*)malloc(sizeof(double)*n*n);
	
	for(i=0; i<n*n; i++)
			a[i] = p[i];
	
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			if(i==j) inv[i*n+j] = 1.;
			else inv[i*n+j] = 0.;
	
	for(i=0; i<n; i++){
		max = a[i*n+i];
		kmax = i;
		for(j=i+1; j<n; j++){
			if(fabs(max)<fabs(a[j*n+i])){
				kmax = j;
				max = a[j*n+i];
			}
			
			if(kmax!=i) l++;
		}
		
		for(j=0; j<n; j++){
			buf = a[i*n+j];
			a[i*n+j] = a[kmax*n+j];
			a[kmax*n+j] = buf;
			
			buf = inv[i*n+j];
			inv[i*n+j] = inv[kmax*n+j];
			inv[kmax*n+j] = buf;
		}
		
		for(k=0; k<n; k++){
			q = a[k*n+i];
			r = a[i*n+i];
//			if(r == 0.) return 0.;
			if(k!=i)
				for(j=0; j<n; j++){
					a[k*n+j] -= a[i*n+j]*q/r;
					inv[k*n+j] -= inv[i*n+j]*q/r;
				}
		}
	}
	
	for(i=0; i<n; i++){
		r = a[i*n+i];
		for(j=0; j<n; j++){
			a[i*n+j] /= r;
			inv[i*n+j] /= r;
		}
			
	}
	
//	for(i=0; i<n; i++)
//		for(j=0; j<n; j++)
//			if(j==n-1) printf("%7.4e\n", a[n*i+j] );
//			else printf("%7.4e\t", a[n*i+j] );
	
	free(a);
	
	return 0;
}

int freq(double* hess, double* mass) {
	double* hess0, * eigenr, * eigenvec, * red_mass, buf;

	int i, j, a, b, N;

	N = N_ATOM;

	hess0 = (double*)calloc(3 * N * (6 * N + 2), sizeof(double));
	eigenr = hess0 + 9 * N * N;
	red_mass = eigenr + 3 * N;
	eigenvec = red_mass + 3 * N;

	for (i = 0; i < N; i++)
		for (a = 0; a < 3; a++)
			for (j = 0; j < N; j++)
				for (b = 0; b < 3; b++)
					hess0[3 * N * (3 * i + a) + 3 * j + b] = hess[3 * N * (3 * i + a) + 3 * j + b] / sqrt(mass[i] * mass[j]);

	if (rs(3 * N, hess0, eigenr, true, eigenvec) != 0)
		return 1;

	for (j = 0; j < 3 * N; j++)
		for (i = 0; i < N; i++)
			for (a = 0; a < 3; a++)
				eigenvec[3 * N * j + 3 * i + a] /= sqrt(mass[i]);

	for (j = 0; j < 3 * N; j++) {
		buf = 0.;
		for (i = 0; i < N; i++)
			for (a = 0; a < 3; a++)
				red_mass[j] += eigenvec[3 * N * j + 3 * i + a] * eigenvec[3 * N * j + 3 * i + a];

		red_mass[j] = 1. / red_mass[j];
	}

	printf("\n\tFREQUENCIES ANALYSIS\n\n");

	for (i = 0; i < 3 * N; i++)
		if (eigenr[i] < 0)
			printf("\t%10.3fi", 521.4495797 * sqrt(-eigenr[i]));		// from eV/(aem A^2) to cm^-1
		else
			printf("\t%10.3f", 521.4495797 * sqrt(eigenr[i]));			// from eV/(aem A^2) to cm^-1

	//	for(i=0; i<3*N; i++)
	//		if(eigenr[i]<0)
	//			printf("%f I\n",5140.27701398*sqrt(-eigenr[i]));		// from hartree/(aem Bohr^2) to cm^-1
	//		else
	//			printf("%f\n",5140.27701398*sqrt(eigenr[i]));			// from hartree/(aem Bohr^2) to cm^-1

	printf("\n\n");

	for (i = 0; i < 3 * N; i++)
		printf("\t%10.3f", red_mass[i]);

	printf("\n\n");

	for (i = 0; i < 3 * N; i++)
		for (j = 0; j < 3 * N; j++)
			if (j == 3 * N - 1) printf("\t%10.5e\n", eigenvec[3 * N * j + i]);
			else printf("\t%10.5e", eigenvec[3 * N * j + i]);

	printf("\n");

	free(hess0);

	return 0;
}

int locate_min(double* r1, int NSTEP, double OPTTOL) {
	double* grad, * hess0, * eigenr, * eigenvec, * dr;
	int i, j, k, N;

	N = N_ATOM;

	grad = (double*)calloc(9 * N * (2 * N + 1), sizeof(double));

	hess0 = grad + 3 * N;
	eigenr = hess0 + 9 * N * N;
	eigenvec = eigenr + 3 * N;
	dr = eigenvec + 9 * N * N;

	if ((OPTTOL <= 0.) || (OPTTOL > 1e-1))
		OPTTOL = 1.0e-3;

	for (k = 0; k < NSTEP; k++) {
		double U;
		//		U_grad_hess_LEPS_dist(r, grad, hess0);
		//U_grad_hess_net_cart(r1, grad, hess0);
		U = U_grad_hess_cart(r1, grad, hess0);
		double buf = 0., grms;

		grms = 0.;
		for (i = 0; i < 3 * N; i++) {
			if (buf < fabs(grad[i])) buf = fabs(grad[i]);
			grms += grad[i] * grad[i];
			dr[i] = 0.;
		}

		grms = sqrt(grms / (3 * N));

		printf("NSTEP=\t%d\n", k + 1);
		printf("Energy=\t%.10f\teV\n", U);
		printf("Max.gradient=\t%7.3e\tRMS_gradient=\t%7.3e\n\n", buf, grms);

		if ((buf < OPTTOL) && (grms < OPTTOL / 3.)) break;

		//		invert_matrix(hess0, a1, DIMD);

		rs(3 * N, hess0, eigenr, true, eigenvec);
		//hess_mw_eigensystem(r1, hess0, eigenr, eigenvec);

		for (i = 0; i < 3 * N; i++)
			for (j = 0; j < 3 * N; j++)
				if (fabs(eigenr[i]) > 1e-2) dr[j] += eigenvec[i * 3 * N + j] * grad[j] * eigenvec[i * 3 * N + j] / eigenr[i];

		for (i = 0; i < 3 * N; i++)
			r1[i] -= 0.005 * dr[i];

		printf("\tCOORDIANTS OF ATOMS (IN ANGSTROM)\n\n");
		for (i = 0; i < N; i++)
			printf("\tATOM%d\t%17.10f%17.10f%17.10f\n", i + 1, r1[3 * i], r1[3 * i + 1], r1[3 * i + 2]);

		printf("\n");
	}

	if (k >= NSTEP) {
		printf("\n\tGEOMETRY WAS NOT CONVERGED\n\n\tPLEASE TRY WITH ANOTHER POINT OR INCREASE ITERATION STEPS\n\n");
		free(grad);
		return 1;
	}

	free(grad);

	return 0;
}

int e_ijk(int i, int j, int k) {
	if (i == j) return 0;
	if (i == k) return 0;
	if (j == k) return 0;

	if ((i == 0) && (j == 1) && (k == 2))
		return 1;
	if ((i == 1) && (j == 2) && (k == 0))
		return 1;
	if ((i == 2) && (j == 0) && (k == 1))
		return 1;

	if ((i == 1) && (j == 0) && (k == 2))
		return -1;
	if ((i == 0) && (j == 2) && (k == 1))
		return -1;
	if ((i == 2) && (j == 1) && (k == 0))
		return -1;

	printf("\n\nIncorrect index of Levi-Civita tensor\n\n");
	return 0;
}

int freq_analysis(int num_atom, double* charge, double* coord, double* hess, double* eigene, double* eigenh, double* red_m, double B0[3], double* sigma) {
	if (num_atom == 1) {
		B0[0] = 0.;
		B0[1] = 0.;
		B0[2] = 0.;
		*sigma = 1.;
		return 0;
	}

	if (num_atom <= 0) {
		printf("\nAtom number is not positive. Critical error.\n");
		return 1;
	}

	double Rc[3], M, Ixx[3], II[9], Ivec[9], * D, * Pa, a1[6], dbuf, * hess1;
	char PG[10];
	int i, j, k, l, d, d1, * tr_rot, Nvib = 3 * num_atom - 6;

	printf("\n\n");

	point_group_calc(num_atom, charge, coord, PG, 0.02);
	*sigma = (double)group_order(PG);

	if (num_atom == 2)
		Nvib = 1;
	else
		if ((strcmp("Cinfv", PG) == 0) || (strcmp("Dinfh", PG) == 0))
			Nvib++;

//	Nvib++;



	//mass = (double*)calloc(num_atom, sizeof(double));
	tr_rot = (int*)calloc(3 * num_atom, sizeof(int));
	hess1 = (double*)calloc(9 * num_atom * num_atom, sizeof(double));
	D = (double*)calloc(3 * num_atom * 6, sizeof(double));
	Pa = (double*)calloc(3 * num_atom, sizeof(double));

	M = 0.;
	for (i = 0; i < num_atom; i++) {
		//mass[i] = mass_atom(charge[i]);
		M += mass[i];
		//		printf("%8.3f", mass_ts[i]);
	}

	for (j = 0; j < 3; j++) {
		Rc[j] = 0.;
		for (i = 0; i < num_atom; i++)
			Rc[j] += mass[i] * coord[3 * i + j];

		Rc[j] = Rc[j] / (M * bohr_A);
	}

	for (i = 0; i < 9; i++)
		II[i] = 0.;

	for (i = 0; i < num_atom; i++) {
		II[0] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));	// Ixx
		II[4] += mass[i] * ((coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));		// Iyy
		II[8] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]));		// Izz

		II[1] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixy
		II[2] -= mass[i] * (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixz
		II[5] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 2] / bohr_A - Rc[2]);	// Iyz
	}

	II[3] = II[1];		// Iyx
	II[6] = II[2];		// Izx
	II[7] = II[5];		// Izy

	rs(3, II, Ixx, true, Ivec);

	printf("\n   The moments of inertia (in a.m.u.*Bohr^2)\n%10.3f%10.3f%10.3f\n", Ixx[0], Ixx[1], Ixx[2]);

	for (i = 0; i < 3; i++)
		if (Ixx[i] > 1e-1) B0[i] = 0.5 * h1_c * h1_c / (Ixx[i] * amu_kg * bohr_m * bohr_m) / kb_c;
		else B0[i] = 0.;
	printf("\n   The rotational constants B=h^2/(8 Pi^2 I k_b) (in K)\n%11.5f%11.5f%11.5f\n", B0[0], B0[1], B0[2]);

	printf("\n   The rotational symmetry number %6.1f\n\n", *sigma);

	for (i = 0; i < num_atom; i++)
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				if (k == j) D[3 * num_atom * j + 3 * i + k] = sqrt(mass[i]);
				else D[3 * num_atom * j + 3 * i + k] = 0.;

	for (i = 3; i < 6; i++)
		for (j = 0; j < num_atom; j++)
			for (k = 0; k < 3; k++)
				for (d = 0; d < 3; d++)
					for (d1 = 0; d1 < 3; d1++)
						D[3 * num_atom * i + 3 * j + k] += sqrt(mass[j]) * (e_ijk(k, d, d1) * Ivec[3 * (i - 3) + d] * (coord[3 * j + d1] / bohr_A - Rc[d1]));

	//	printf("\n");
	//	for(j=0; j<6; j++){
	//		for(i=0; i<3*num_atom; i++)
	//			printf("%7.3f", D[3*num_atom*j+i]);
	//		printf("\n");
	//	}

	for (j = 0; j < 6; j++) {
		dbuf = 0.;
		for (i = 0; i < 3 * num_atom; i++)
			dbuf += D[3 * num_atom * j + i] * D[3 * num_atom * j + i];

		for (i = 0; i < 3 * num_atom; i++)
			if (dbuf >= 1e-3) D[3 * num_atom * j + i] /= sqrt(dbuf);
	}
	//	printf("\n");
	//	for(j=0; j<6; j++){
	//		for(i=0; i<3*num_atom; i++)
	//			printf("%7.3f", D[3*num_atom*j+i]);
	//		printf("\n");
	//	}

	for (i = 0; i < 3 * num_atom; i++)
		for (j = 0; j < 3 * num_atom; j++)
			hess1[i * 3 * num_atom + j] = hess[i * 3 * num_atom + j] / (sqrt(mass[i / 3]) * sqrt(mass[j / 3]));

	//	for(i=0; i<9*num_atom*num_atom; i++)
	//		hess1[i] = hess[i];

	rs(3 * num_atom, hess1, eigene, true, eigenh);

	for (i = 0; i < 3 * num_atom; i++) {
		for (k = 0; k < 6; k++)
			a1[k] = 0.;
		for (j = 0; j < 3 * num_atom; j++)
			for (k = 0; k < 6; k++)
				a1[k] += eigenh[3 * num_atom * i + j] * D[3 * num_atom * k + j];


		for (j = 0; j < 3 * num_atom; j++) {
			dbuf = 0.;
			for (k = 0; k < 6; k++)
				dbuf += a1[k] * D[3 * num_atom * k + j];
			Pa[i] += sqr(eigenh[3 * num_atom * i + j] - dbuf);
		}

		Pa[i] = sqrt(Pa[i]);
	}


	for (i = 0; i < 3 * num_atom; i++)
		tr_rot[i] = i;

	for (i = 0; i < 3 * num_atom; i++) {
		for (j = 0; j < 3 * num_atom - 1; j++) {
			if (fabs(Pa[tr_rot[j]]) < fabs(Pa[tr_rot[j + 1]])) {
				int ibuf = tr_rot[j];				// создали дополнительную переменную
				tr_rot[j] = tr_rot[j + 1];			// меняем местами
				tr_rot[j + 1] = ibuf;				// значения элементов
			}
		}
	}

	free(D);

	for (i = 0; i < 3 * num_atom; i++)
		if (eigene[i] >= 0) eigene[i] = sqrt(eigene[i]) * 521.4495797;
		else eigene[i] = -sqrt(-eigene[i]) * 521.4495797;

	for (i = Nvib; i < 3 * num_atom; i++)
		eigene[tr_rot[i]] = 0.;

	for (i = Nvib; i < 3 * num_atom; i++) {
		for (j = Nvib; j < 3 * num_atom - 1; j++) {
			if (tr_rot[j] > tr_rot[j + 1]) {
				int ibuf = tr_rot[j];				// создали дополнительную переменную
				tr_rot[j] = tr_rot[j + 1];			// меняем местами
				tr_rot[j + 1] = ibuf;				// значения элементов
			}
		}
	}

	printf("\n   Modes ");
	for (i = Nvib; i < 3 * num_atom; i++)
		if (i == 3 * num_atom - 1) printf("%d ", tr_rot[i] + 1);
		else printf("%d, ", tr_rot[i] + 1);
	printf("was taken as translational-rotational modes.\n   Please check the program choice.\n\n\n   Mode number");
	for (i = 0; i < 3 * num_atom; i++)
		printf("%10d", i + 1);
	printf("\n   Frequencies    ");
	for (i = 0; i < 3 * num_atom; i++)
		printf("%10.2f", eigene[i]);
	printf("\n   Reduced masses ");

	for (j = 0; j < 3 * num_atom; j++) {
		red_m[j] = 0.;
		for (i = 0; i < 3 * num_atom; i++)
			red_m[j] += eigenh[3 * num_atom * j + i] * eigenh[3 * num_atom * j + i] / mass[i / 3];

		red_m[j] = 1. / red_m[j];
	}
	for (i = 0; i < 3 * num_atom; i++)
		printf("%10.4f", red_m[i]);
	printf("\n\n");

	for (i = 0; i < 3 * num_atom; i++)
		for (j = 0; j < 3 * num_atom; j++) {
			eigenh[3 * num_atom * j + i] /= sqrt(mass[i / 3]);
			if (j == 0) {
				if (i % 3 == 0) printf("  %2d   %3s    X   ", i / 3 + 1, name_particle[i / 3]);
				if (i % 3 == 1) printf("              Y   ");
				if (i % 3 == 2) printf("              Z   ");
			}
			if (j == 3 * num_atom - 1)
				printf("%10.6f%\n", eigenh[3 * num_atom * j + i]);
			else
				printf("%10.6f%", eigenh[3 * num_atom * j + i]);
		}
	printf("\n   TRLDP*         ");
	for (i = 0; i < 3 * num_atom; i++)
		printf("%10.6f", Pa[i]);
	printf("\n\n\n   *TRLDP --- translational-rotational linear dependece parameter\n   The least values related to translational and/or rotational modes.\n   Warning! This parameter is not Sayvetz translational/rotational conditions\n");
	printf("   Trans.-rot. analysis based on vibrational analysis in \"Gaussian\"(c) by Joseph W. Ochterski.\n   You can find more information on site https://gaussian.com\n\n");

	dbuf = 0.;
	for (i = 0; i < 3 * num_atom; i++)
		if (eigene[i] > 0)
			dbuf += eigene[i];
	dbuf *= 0.5;

	printf("Zero point energy\n   %11.6f eV      %20.6f cm^-1\n   %11.6f kcal/mol%20.6f kJ/mol\n\n\n", dbuf * 1.239945E-4, dbuf, dbuf * 2.859170E-3, dbuf * 1.196277E-2);

	free(Pa);
	free(tr_rot);
	//free(mass);
	free(hess1);

	return 0;
}

int k_A_BC_CTST(double* r0, double* T_list, int N_T) {
	init_U();

	int i, j, k, state;
	double grad[9], hess[81], eigene[9], eigenv[81], B0[3], mu[9], sigma, E_atom, we, m_AB, Ea, *Q_TS, *Q_R, M = 0, theta_BC, K_CTST;
	double alpha_tr, alpha_r, alpha_v;
	double a_AHO = 0.72, b_AHO = 0.15;

	state = locate_min(r0, 30, 1e-6);

	if(state != 0){
		printf("     Cannot locate TS.\n     Give another initial geometry and restart.\n");
		return 1;
	}

	E_atom = U_grad_hess_cart(r0, grad, hess);
	Ea = E_atom - Ev_BC[0];
	E_atom = D_lim - E_atom;
	we = 2. * (Ev_AB[0] - Te_AB);
	m_AB = mass[0]*mass[1]/(mass[0]+mass[1]);
	for(i=0; i<N_ATOM; i++)
		M += mass[i];
	theta_BC = 24*(mass[1]+mass[2])/(mass[1]*mass[2]) / sqr(re_BC);


	alpha_tr = pow(2.0 * M_PI * amu_kg * kb_c, 1.5) * pow(h_c, -3) * 1.0e-6; 		// constant to calculate translational statistical sum per 1 cm^-3
	alpha_r = bohr_m / h1_c * sqrt(2.0 * amu_kg * kb_c);					// constant to calculate rotational statistical sum per 1 cm^-3
	alpha_v = h_c * c_c * 1.0e2 / kb_c;										// constant to convert cm^-1 to K for calculation vibrational stat.sum

	freq_analysis(N_ATOM, (double*)Ze, r0, hess, eigene, eigenv, mu, B0, &sigma);

	for(i=0; i<3*N_ATOM; i++)
		if(eigene[i]>1e-2)
			Ea += 0.5*eigene[i]*1.239945E-4;  // in eV
	
	Q_TS = (double*)calloc(N_T, sizeof(double));
	Q_R = (double*)calloc(N_T, sizeof(double));

	printf("\n\n         Stat.sum data for located TS\nT, K\tQ_transl, cm^-3\tQ_rot\tQ_vib(HO)\tQ_vib(AHO)\tQ_vib(AHO)/Q_vib(HO)\tQ_total, cm^-3\n");

	for(i=0; i<N_T; i++){
		double T = T_list[i];
		double Q1 = 1., Q2 = 1., Q3, Emax;
		// Translational stat.sum 
		Q_TS[i] = alpha_tr * pow(M * T, 1.5);
		printf("%.2f\t%20.10e\t", T, Q_TS[i]);

		// Rotational stat.sum
		if (N_ATOM == 1) {
			Q1 = 1.;
			Q_TS[i] *= Q1;
			printf("%.10e\t", Q1);
		}
		else {
			if ((B0[0] < 1e-7) && (B0[1] < 1e-7) && (B0[0] < 1e-7)) {
				printf("\nIt seems that program deal with one atom.\nPlease check all data and restart\n");
				return 1;
			}
			else
				if (B0[0] < 1e-4) Q1 = sqrt(T / B0[1]) * sqrt(T / B0[2]) / sigma;
				else Q1 = sqrt(M_PI) * sqrt(T / B0[0]) * sqrt(T / B0[1]) * sqrt(T / B0[2]) / sigma;

			Q_TS[i] *= Q1;
			printf("%.10e\t", Q1);
		}

		Q1 = 1.;
		// Vibrational stat.sum
		// HO-approximation

		for (j = 0; j < 3 * N_ATOM; j++)
			if (eigene[j] > 1.0e-2) Q1 /= (1.0 - exp(-alpha_v  * eigene[j] / T));

		// AHO-approximation
		for (j = 0; j < 3 * N_ATOM; j++)
			if (eigene[j] > 1.0e-2) {
				Q3 = 0.;
				Emax = E_atom * 8.064877E+3 / (N_ATOM - 1) * pow(eigene[j] / 8.064877E+3 / we, a_AHO) * pow(mu[j] / m_AB, b_AHO);		// 8.064877E+3 --- from eV to cm^-1
				//				Emax=1e6;
				for (k = 0; k <= 2 * Emax / eigene[j] - 0.5; k++)
					Q3 += exp(-alpha_v / T * (eigene[j] * (k + 0.5) - 0.25 * sqr(eigene[j]) / Emax * sqr(k + 0.5)));
				Q3 /= exp(-alpha_v / T * (eigene[j] * 0.5 - 0.25 * 0.25 * sqr(eigene[j]) / Emax));
				Q2 *= Q3;
			}

		printf("%.10e\t%.10e\t%.10e\t", Q1, Q2, Q2 / Q1);

		Q_TS[i] *= Q2;
		printf("%.10e\n", Q_TS[i]);

	}
		


	printf("\n\n         Stat.sum data for reactant A+BC\nT, K\tQ_transl, cm^-3\tQ_rot\tQ_vib\tQ_total, cm^-3\n");

	for (i = 0; i < N_T; i++) {
		double T = T_list[i];
		double Q1 = 1., Q2 = 1., Q3, Emax;
		// Translational stat.sum 
		Q_R[i] = alpha_tr * pow(mass[0] * T, 1.5) * alpha_tr * pow((mass[1] + mass[2]) * T, 1.5);
		printf("%.2f\t%20.10e\t", T, Q_R[i]);

		// Rotational stat.sum
		

		if((Ze[1]==Ze[2])&&(mass[1]==mass[2]))
			Q1 = 0.5*T/theta_BC;
		else
			Q1 = T/theta_BC;
		Q_R[i] *= Q1;
		printf("%.10e\t", Q1);

		Q1 = 0.;
		// Vibrational stat.sum
		for(j=0; j<dim_Ev_BC; j++){
			Q1 += exp(-Ev_BC[j] * 1.160363E+4 / T);			// 1.160363E+4 --- from eV to K
		}
		Q1 /= exp(-Ev_BC[0] * 1.160363E+4 / T);

		printf("%.10e\t", Q1);

		Q_R[i] *= Q1;
		printf("%.10e\n", Q_R[i]);

	}

	printf("\n\nT, K\tk_CTST, cm^3/(mol s)\tGamma\n");

	for(i=0; i<N_T; i++){
		double gamma;
		gamma = 1. + 1. / 24. * sqr(eigene[0] * 1.438786 / T_list[i]);
		K_CTST = kb_c*T_list[i] / h_c * Q_TS[i] / Q_R[i] * exp(-Ea * 1.160363E+4 / T_list[i]) * Na;
		printf("%.2f\t%.10e\t%.10e\n", T_list[i], K_CTST, gamma);
	}


	free(Q_R);
	free(Q_TS);

	deinit_U();

	return 0;
}

int locate_CP(double* r1, int NSTEP, double OPTTOL, bool CP) { // CP --- false for local minimum, true for saddle point
	
	double* grad, * hess0, * eigenr, * eigenvec, * dr, trmax = 0.0002, norm_dr;
	int i, j, k, N;

	N = N_ATOM;

	grad = (double*)calloc(9 * N * (2 * N + 1), sizeof(double));

	hess0 = grad + 3 * N;
	eigenr = hess0 + 9 * N * N;
	eigenvec = eigenr + 3 * N;
	dr = eigenvec + 9 * N * N;

	if ((OPTTOL <= 0.) || (OPTTOL > 1e-1))
		OPTTOL = 5.2e-3;

	for (k = 0; k < NSTEP; k++) {
		double U;
		//		U_grad_hess_LEPS_dist(r, grad, hess0);
		//U_grad_hess_net_cart(r1, grad, hess0);
		U = U_grad_hess_cart(r1, grad, hess0);
		double buf = 0., grms;

		grms = 0.;
		for (i = 0; i < 3 * N; i++) {
			if (buf < fabs(grad[i])) buf = fabs(grad[i]);
			grms += grad[i] * grad[i];
			dr[i] = 0.;
		}

		grms = sqrt(grms / (3 * N));

		printf("NSTEP=\t%d\n", k + 1);
		printf("Energy=\t%.10f\teV\n", U);
		printf("Max.gradient=\t%7.3e\tRMS_gradient=\t%7.3e\n\n", buf, grms);

		if ((buf < OPTTOL) && (grms < OPTTOL / 3.))
			break;

		//		invert_matrix(hess0, a1, DIMD);

		//rs(3 * N, hess0, eigenr, true, eigenvec);
		hess_eigensystem(r1, hess0, eigenr, eigenvec);

		if(CP){

			double min=1000.;
			int count = 0;
			for(i=0; i<3*N; i++)
				if((min>eigenr[i])&&(eigenr[i]!=0.)){
					min = eigenr[i];
					count = i;
				}
			
			for(i=0; i<3*N; i++)
				if(i!=count)
					eigenr[i] = fabs(eigenr[i]);
				else
					eigenr[i] = -fabs(eigenr[i]);
		}
		else {
			for(i=0; i<3*N; i++)
				eigenr[i] = fabs(eigenr[i]);
		}


		for (i = 0; i < 3 * N; i++)
			for (j = 0; j < 3 * N; j++)
				if (fabs(eigenr[i]) > 1e-2) 
					dr[j] += eigenvec[i * 3 * N + j] * grad[j] * eigenvec[i * 3 * N + j] / eigenr[i];

		norm_dr = 0.;
		for(i=0; i<3*N; i++)
			norm_dr += dr[i]*dr[i];
		norm_dr = sqrt(norm_dr);

		if(norm_dr>trmax)
			for(i=0; i<3*N; i++)
				dr[i] *= trmax/norm_dr;

		for (i = 0; i < 3 * N; i++)
			r1[i] -= dr[i];

		printf("\tCOORDIANTS OF ATOMS (IN ANGSTROM)\n\n");
		for (i = 0; i < N; i++)
			printf("\tATOM%d\t%17.10f\t%17.10f\t%17.10f\n", i + 1, r1[3 * i], r1[3 * i + 1], r1[3 * i + 2]);

		printf("\n");
	}

	if (k >= NSTEP) {
		printf("\n\tGEOMETRY WAS NOT CONVERGED\n\n\tPLEASE TRY WITH ANOTHER POINT OR INCREASE ITERATION STEPS\n\n");
		free(grad);
		return 1;
	}

	free(grad);

	return 0;
}

int hess_eigensystem(double* coord, double* hess, double* eigene, double* eigenh) {

	char PG[10];
	int Nvib = 3, i, j, k, l, num_atom = 3, * tr_rot, d, d1;
	double sigma, * hess1, * D, * Pa, M, Rc[3], II[9], Ixx[3], Ivec[9], dbuf, a1[6];

	point_group_calc(num_atom, (double*)Ze, coord, PG, 0.02);
	sigma = (double)group_order(PG);

	if ((strcmp("Cinfv", PG) == 0) || (strcmp("Dinfh", PG) == 0))
		Nvib++;

	tr_rot = (int*)calloc(3 * num_atom, sizeof(int));
	hess1 = (double*)calloc(9 * num_atom * num_atom, sizeof(double));
	D = (double*)calloc(3 * num_atom * 6, sizeof(double));
	Pa = (double*)calloc(3 * num_atom, sizeof(double));

	M = 0.;
	for (i = 0; i < num_atom; i++)
		M += mass[i];

	for (j = 0; j < 3; j++) {
		Rc[j] = 0.;
		for (i = 0; i < num_atom; i++)
			Rc[j] += mass[i] * coord[3 * i + j];

		Rc[j] = Rc[j] / (M * bohr_A);
	}

	for (i = 0; i < 9; i++)
		II[i] = 0.;

	for (i = 0; i < num_atom; i++) {
		II[0] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));	// Ixx
		II[4] += mass[i] * ((coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));		// Iyy
		II[8] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]));		// Izz

		II[1] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixy
		II[2] -= mass[i] * (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixz
		II[5] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 2] / bohr_A - Rc[2]);	// Iyz
	}

	II[3] = II[1];		// Iyx
	II[6] = II[2];		// Izx
	II[7] = II[5];		// Izy

	rs(3, II, Ixx, true, Ivec);

	for (i = 0; i < num_atom; i++)
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				if (k == j) D[3 * num_atom * j + 3 * i + k] = sqrt(mass[i]);
				else D[3 * num_atom * j + 3 * i + k] = 0.;

	for (i = 3; i < 6; i++)
		for (j = 0; j < num_atom; j++)
			for (k = 0; k < 3; k++)
				for (d = 0; d < 3; d++)
					for (d1 = 0; d1 < 3; d1++)
						D[3 * num_atom * i + 3 * j + k] += sqrt(mass[j]) * (e_ijk(k, d, d1) * Ivec[3 * (i - 3) + d] * (coord[3 * j + d1] / bohr_A - Rc[d1]));

	for (j = 0; j < 6; j++) {
		dbuf = 0.;
		for (i = 0; i < 3 * num_atom; i++)
			dbuf += D[3 * num_atom * j + i] * D[3 * num_atom * j + i];

		for (i = 0; i < 3 * num_atom; i++)
			if (dbuf >= 1e-3) D[3 * num_atom * j + i] /= sqrt(dbuf);
	}

	for (i = 0; i < 3 * num_atom; i++)
		for (j = 0; j < 3 * num_atom; j++)
			hess1[i * 3 * num_atom + j] = hess[i * 3 * num_atom + j];

	rs(3 * num_atom, hess1, eigene, true, eigenh);

	for (i = 0; i < 3 * num_atom; i++) {
		for (k = 0; k < 6; k++)
			a1[k] = 0.;
		for (j = 0; j < 3 * num_atom; j++)
			for (k = 0; k < 6; k++)
				a1[k] += eigenh[3 * num_atom * i + j] * D[3 * num_atom * k + j];


		for (j = 0; j < 3 * num_atom; j++) {
			dbuf = 0.;
			for (k = 0; k < 6; k++)
				dbuf += a1[k] * D[3 * num_atom * k + j];
			Pa[i] += sqr(eigenh[3 * num_atom * i + j] - dbuf);
		}

		Pa[i] = sqrt(Pa[i]);
	}

	for (i = 0; i < 3 * num_atom; i++)
		tr_rot[i] = i;

	for (i = 0; i < 3 * num_atom; i++) {
		for (j = 0; j < 3 * num_atom - 1; j++) {
			if (fabs(Pa[tr_rot[j]]) < fabs(Pa[tr_rot[j + 1]])) {
				int ibuf = tr_rot[j];				// создали дополнительную переменную
				tr_rot[j] = tr_rot[j + 1];			// меняем местами
				tr_rot[j + 1] = ibuf;				// значения элементов
			}
		}
	}

	for (i = Nvib; i < 3 * num_atom; i++)
		eigene[tr_rot[i]] = 0.;

	free(D);
	free(Pa);
	free(tr_rot);
	free(hess1);
	return 0;
}

int hess_mw_eigensystem(double* coord, double* hess, double* eigene, double* eigenh) {
	
	char PG[10];
	int Nvib = 3, i, j, k, l, num_atom = 3, *tr_rot, d, d1;
	double sigma, *hess1, *D, *Pa, M, Rc[3], II[9], Ixx[3], Ivec[9], dbuf, a1[6];

	point_group_calc(num_atom, (double*)Ze, coord, PG, 0.02);
	sigma = (double)group_order(PG);

	if ((strcmp("Cinfv", PG) == 0) || (strcmp("Dinfh", PG) == 0))
		Nvib++;

	tr_rot = (int*)calloc(3 * num_atom, sizeof(int));
	hess1 = (double*)calloc(9 * num_atom * num_atom, sizeof(double));
	D = (double*)calloc(3 * num_atom * 6, sizeof(double));
	Pa = (double*)calloc(3 * num_atom, sizeof(double));

	M = 0.;
	for (i = 0; i < num_atom; i++)
		M += mass[i];

	for (j = 0; j < 3; j++) {
		Rc[j] = 0.;
		for (i = 0; i < num_atom; i++)
			Rc[j] += mass[i] * coord[3 * i + j];

		Rc[j] = Rc[j] / (M * bohr_A);
	}

	for (i = 0; i < 9; i++)
		II[i] = 0.;

	for (i = 0; i < num_atom; i++) {
		II[0] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));	// Ixx
		II[4] += mass[i] * ((coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]) + (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i + 2] / bohr_A - Rc[2]));		// Iyy
		II[8] += mass[i] * ((coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 1] / bohr_A - Rc[1]) + (coord[3 * i] / bohr_A - Rc[0]) * (coord[3 * i] / bohr_A - Rc[0]));		// Izz

		II[1] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixy
		II[2] -= mass[i] * (coord[3 * i + 2] / bohr_A - Rc[2]) * (coord[3 * i] / bohr_A - Rc[0]);	// Ixz
		II[5] -= mass[i] * (coord[3 * i + 1] / bohr_A - Rc[1]) * (coord[3 * i + 2] / bohr_A - Rc[2]);	// Iyz
	}

	II[3] = II[1];		// Iyx
	II[6] = II[2];		// Izx
	II[7] = II[5];		// Izy

	rs(3, II, Ixx, true, Ivec);

	for (i = 0; i < num_atom; i++)
		for (k = 0; k < 3; k++)
			for (j = 0; j < 3; j++)
				if (k == j) D[3 * num_atom * j + 3 * i + k] = sqrt(mass[i]);
				else D[3 * num_atom * j + 3 * i + k] = 0.;

	for (i = 3; i < 6; i++)
		for (j = 0; j < num_atom; j++)
			for (k = 0; k < 3; k++)
				for (d = 0; d < 3; d++)
					for (d1 = 0; d1 < 3; d1++)
						D[3 * num_atom * i + 3 * j + k] += sqrt(mass[j]) * (e_ijk(k, d, d1) * Ivec[3 * (i - 3) + d] * (coord[3 * j + d1] / bohr_A - Rc[d1]));

	for (j = 0; j < 6; j++) {
		dbuf = 0.;
		for (i = 0; i < 3 * num_atom; i++)
			dbuf += D[3 * num_atom * j + i] * D[3 * num_atom * j + i];

		for (i = 0; i < 3 * num_atom; i++)
			if (dbuf >= 1e-3) D[3 * num_atom * j + i] /= sqrt(dbuf);
	}

	for (i = 0; i < 3 * num_atom; i++)
		for (j = 0; j < 3 * num_atom; j++)
			hess1[i * 3 * num_atom + j] = hess[i * 3 * num_atom + j] / (sqrt(mass[i / 3]) * sqrt(mass[j / 3]));

	rs(3 * num_atom, hess1, eigene, true, eigenh);

	for (i = 0; i < 3 * num_atom; i++) {
		for (k = 0; k < 6; k++)
			a1[k] = 0.;
		for (j = 0; j < 3 * num_atom; j++)
			for (k = 0; k < 6; k++)
				a1[k] += eigenh[3 * num_atom * i + j] * D[3 * num_atom * k + j];


		for (j = 0; j < 3 * num_atom; j++) {
			dbuf = 0.;
			for (k = 0; k < 6; k++)
				dbuf += a1[k] * D[3 * num_atom * k + j];
			Pa[i] += sqr(eigenh[3 * num_atom * i + j] - dbuf);
		}

		Pa[i] = sqrt(Pa[i]);
	}

	for (i = 0; i < 3 * num_atom; i++)
		tr_rot[i] = i;

	for (i = 0; i < 3 * num_atom; i++) {
		for (j = 0; j < 3 * num_atom - 1; j++) {
			if (fabs(Pa[tr_rot[j]]) < fabs(Pa[tr_rot[j + 1]])) {
				int ibuf = tr_rot[j];				// создали дополнительную переменную
				tr_rot[j] = tr_rot[j + 1];			// меняем местами
				tr_rot[j + 1] = ibuf;				// значения элементов
			}
		}
	}

	for (i = Nvib; i < 3 * num_atom; i++)
		eigene[tr_rot[i]] = 0.;

	free(D);
	free(Pa);
	free(tr_rot);
	free(hess1);
	return 0;
}

double Fi(double *r2, double *r1, double *grad1, double *hess, double lam2, double stride, double *dr ) {
	int i, j;
	double grad[9], hess1[81], inv_hess[81], F2;

	for (i = 0; i < 9; i++)
		for (j = 0; j < 9; j++)
			hess1[i * 9 + j] = hess[i * 9 + j] / (sqrt(mass[i / 3]) * sqrt(mass[j / 3]));

	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
			grad[3*i+j] = grad1[3*i+j]/sqrt(mass[i]);
	
	for(i=0; i<9; i++)
		hess1[9*i+i] -= lam2;

	invert_matrix(hess1, inv_hess, 9);

	for(i=0; i<9; i++){
		dr[i] = 0.;
		for(j=0; j<9; j++)
			dr[i] -= inv_hess[9*i+j]*( grad[j] - lam2*(r2[j] - r1[j]) );
	}
	F2 = 0;
	for(i=0; i<9; i++)
		F2 += sqr( r2[i] - r1[i] + dr[i]);

	F2 -= 0.25 * sqr(stride);

	//printf("dr\n");
	//for (i = 0; i < 3; i++) {
	//	printf("%.1f", charge[i]);
	//	for (j = 0; j < 3; j++)
	//		printf("\t%.10f", dr[3 * i + j] / sqrt(mass[i]));
	//	printf("\n");
	//}
	//printf("\n");

	return F2;
}

int MEP_GS2(double* r, double stride, int maxiter, bool TS, bool dir, double opttol) {

	int i, j, microiter, iter;
	double q[9], r1[9], r2[9], dr[9], grad[9], hess[81], U, eigenv[9], eigenh[81], norm_g, lam1, lam2, dlam, irc, F1, F2;
	double dr_max, max_g, U1, U2;
	FILE* xyz;

	for(i=0; i<9; i++)
		q[i] = r[i];

	for(i=0; i<3; i++)
		for(j=0; j<3; j++)
			r[3*i+j] *= sqrt(mass[i]);

	irc = 0.;

	U2 = U_grad_hess_cart(q, grad, hess);
	U1 = U2;

	xyz = fopen("irc.xyz", "w");

	printf("irc, ae\tU, eV\n0.0000\t%.10f\n", U2);
	fprintf(xyz, "irc=\t0.0000\tU=\t%.10f\n", U2);
	for(i=0; i<3; i++){
		fprintf(xyz, "%.1f", Ze[i]);
		for(j=0; j<3; j++)
			fprintf(xyz, "\t%.10f", q[3*i+j]);
		fprintf(xyz, "\n");
	}
	fprintf(xyz, "\n");

	if(TS){
		if(dir){
			hess_mw_eigensystem(q, hess, eigenv, eigenh);
			for(i=0; i<9; i++){
				r1[i] = r[i] + 0.5*stride*eigenh[i];
				r2[i] = r[i] + stride*eigenh[i];
			}
		}
		else {
			hess_mw_eigensystem(q, hess, eigenv, eigenh);
			for(i=0; i<9; i++){
				r1[i] = r[i] - 0.5*stride*eigenh[i];
				r2[i] = r[i] - stride*eigenh[i];
			}
		}
	}
	else {
		norm_g = 0;
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				norm_g += grad[3*i+j]*grad[3*i+j]/mass[i];
		norm_g = sqrt(norm_g);

		for(i=0; i<9; i++){
			r1[i] = r[i] - 0.5*stride*grad[i]/sqrt(mass[i/3])/norm_g;
			r2[i] = r[i] - stride*grad[i]/sqrt(mass[i/3])/norm_g;
		}
	}
	
	for(iter=0; iter<maxiter; iter++){
		irc += stride;
		for(microiter = 0; microiter <20 ; microiter++){

			for(i=0; i<3; i++)
				for(j=0; j<3; j++)
					q[3*i+j] = r2[3*i+j]/sqrt(mass[i]);

			U = U_grad_hess_cart(q, grad, hess);
			hess_mw_eigensystem(q, hess, eigenv, eigenh);
	
			lam1 = eigenv[0] - 0.5;
			lam2 = eigenv[0]-1.5;

			F1 = Fi(r2, r1, grad, hess, lam1, stride, dr);
			
			for(i=0; ; i++){
				F2 = Fi(r2, r1, grad, hess, lam2, stride, dr);
						
				if (fabs(F2) < 1e-10)
					break;

				dlam = F2 * (lam2 - lam1) / (F2 - F1);

				F1 = F2;
				lam1 = lam2;
				lam2 = lam2 - dlam;
			}

			dr_max = 0;
			for(i=0; i<9; i++)
				if(dr_max<fabs(dr[i]))
					dr_max = fabs(dr[i]);

			if(dr_max<1e-9)
				break;

			for(i=0; i<9; i++)
				r2[i] += dr[i];
		}
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				q[3*i+j] = r2[3*i+j]/sqrt(mass[i]);
		U = U_grad_cart(q, grad);

		U1 = U2;
		U2 = U;

		printf("%.5f\t%.10f\n", irc, U);
		fprintf(xyz, "irc=\t%.5f\tU=\t%.10f\n", irc, U);
		for (i = 0; i < 3; i++) {
			fprintf(xyz, "%.1f", Ze[i]);
			for (j = 0; j < 3; j++)
				fprintf(xyz, "\t%.10f", q[3*i+j]);
			fprintf(xyz, "\n");
		}
		fprintf(xyz, "\n");

		norm_g = 0;
		max_g = 0;
		for(i=0; i<3; i++)
			for(j=0; j<3; j++){
				norm_g += grad[3*i+j]*grad[3*i+j]/mass[i];
				if(max_g<fabs(grad[3*i+j]))
					max_g = fabs(grad[3*i+j]);
			}
				
		norm_g = sqrt(norm_g);
		//printf("%e\t%e\n", norm_g, max_g);

		if(max_g<opttol)
			break;

		if(U1<U2+1e-10)
			break;

		for(i=0; i<9; i++){
			r[i] = r2[i];
			r1[i] = r[i] - 0.5*stride*grad[i]/sqrt(mass[i/3])/norm_g;
			r2[i] = r[i] - stride*grad[i]/sqrt(mass[i/3])/norm_g;
		}
	}

	fclose(xyz);

	return 0;
}

int group_order(char* group) {
	if (strcmp("C1", group) == 0) return 1;
	if (strcmp("Ci", group) == 0) return 1;
	if (strcmp("Cs", group) == 0) return 1;
	if (strcmp("Kh", group) == 0) return 1;

	if (strcmp("Cinfv", group) == 0) return 1;
	if (strcmp("Dinfh", group) == 0) return 2;

	if (strcmp("S4", group) == 0) return 2;
	if (strcmp("S6", group) == 0) return 3;
	if (strcmp("S8", group) == 0) return 4;

	if (strcmp("C2", group) == 0) return 2;
	if (strcmp("C3", group) == 0) return 3;
	if (strcmp("C4", group) == 0) return 4;
	if (strcmp("C5", group) == 0) return 5;
	if (strcmp("C6", group) == 0) return 6;
	if (strcmp("C7", group) == 0) return 7;
	if (strcmp("C8", group) == 0) return 8;

	if (strcmp("C2v", group) == 0) return 2;
	if (strcmp("C3v", group) == 0) return 3;
	if (strcmp("C4v", group) == 0) return 4;
	if (strcmp("C5v", group) == 0) return 5;
	if (strcmp("C6v", group) == 0) return 6;
	if (strcmp("C7v", group) == 0) return 7;
	if (strcmp("C8v", group) == 0) return 8;

	if (strcmp("C2h", group) == 0) return 2;
	if (strcmp("C3h", group) == 0) return 3;
	if (strcmp("C4h", group) == 0) return 4;
	if (strcmp("C5h", group) == 0) return 5;
	if (strcmp("C6h", group) == 0) return 6;
	if (strcmp("C7h", group) == 0) return 7;
	if (strcmp("C8h", group) == 0) return 8;

	if (strcmp("D2", group) == 0) return 2;
	if (strcmp("D3", group) == 0) return 3;
	if (strcmp("D4", group) == 0) return 4;
	if (strcmp("D5", group) == 0) return 5;
	if (strcmp("D6", group) == 0) return 6;
	if (strcmp("D7", group) == 0) return 7;
	if (strcmp("D8", group) == 0) return 8;

	if (strcmp("D2", group) == 0) return 4;
	if (strcmp("D3", group) == 0) return 6;
	if (strcmp("D4", group) == 0) return 8;
	if (strcmp("D5", group) == 0) return 10;
	if (strcmp("D6", group) == 0) return 12;
	if (strcmp("D7", group) == 0) return 14;
	if (strcmp("D8", group) == 0) return 16;

	if (strcmp("D2d", group) == 0) return 4;
	if (strcmp("D3d", group) == 0) return 6;
	if (strcmp("D4d", group) == 0) return 8;
	if (strcmp("D5d", group) == 0) return 10;
	if (strcmp("D6d", group) == 0) return 12;
	if (strcmp("D7d", group) == 0) return 14;
	if (strcmp("D8d", group) == 0) return 16;

	if (strcmp("D2h", group) == 0) return 4;
	if (strcmp("D3h", group) == 0) return 6;
	if (strcmp("D4h", group) == 0) return 8;
	if (strcmp("D5h", group) == 0) return 10;
	if (strcmp("D6h", group) == 0) return 12;
	if (strcmp("D7h", group) == 0) return 14;
	if (strcmp("D8h", group) == 0) return 16;

	if (strcmp("T", group) == 0) return 12;
	if (strcmp("Td", group) == 0) return 12;

	if (strcmp("O", group) == 0) return 24;
	if (strcmp("Oh", group) == 0) return 24;

	if (strcmp("I", group) == 0) return 60;
	if (strcmp("Ih", group) == 0) return 60;

	return 1;
}