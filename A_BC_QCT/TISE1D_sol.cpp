#include "TISE1D_sol.h"

//double V(double x) {
//
//	//return 0.5*x * x;
//	//return 37259.7*(1.-exp(-2.286*(x-0.9706)))*(1.-exp(-2.286*(x-0.9706)));
//	return 4.39*(1.-exp(-2.286*(x-0.9706)))*(1.-exp(-2.286*(x-0.9706)));
//}

int TISE1D_sol(double *V, double x1, double xN, int N, int order, double Emin, int iElim) {

	double h, gamma, *a, *b, *w, *c, Elim;
	int i, j, k;
	lapack_int *ipiv;
	unsigned long int t_buf;

	h = (xN - x1) / N;
	//J = 0;
	//gamma = 2.;
	gamma = 478.450668 * 0.503915;
	Elim = V[iElim];

	a = (double*)calloc((N + 1) * (N + 1), sizeof(double));
	b = (double*)calloc((N + 1) * (N + 1), sizeof(double));
	c = (double*)calloc((N + 1) * (N + 1), sizeof(double));
	w = (double*)calloc((N + 1), sizeof(double));
	ipiv = (lapack_int*)calloc((N + 1), sizeof(lapack_int));

	t_buf = time(NULL);                              // Calculating of current time
	switch(order) {
		case 4:
		// Numerov method for TISE1D 4-th order
			for(i=0; i<N+1; i++){
				a[i*(N+1)+i] = -2./h/h/gamma;
				if(i<N) a[i*(N+1)+i+1] = 1./h/h/gamma;
				if(i>0) a[i*(N+1)+i-1] = 1./h/h/gamma;
			}

			for(i=0; i<N+1; i++){
				b[i*(N+1)+i] = 10./12.;
				if(i<N) b[i*(N+1)+i+1] = 1./12;
				if(i>0) b[i*(N+1)+i-1] = 1./12;
			}
			break;
		case 6:
			// Numerov method for TISE1D 6-th order
				for(i=0; i<N+1; i++){
					a[i*(N+1)+i] = -2.5/h/h/gamma;
					if(i<N) a[i*(N+1)+i+1] = 4./3./h/h/gamma;
					if(i>0) a[i*(N+1)+i-1] = 4./3./h/h/gamma;
					if(i>1) a[i*(N+1)+i-2] = -1./12./h/h/gamma;
					if(i<N-1) a[i*(N+1)+i+2] = -1./12./h/h/gamma;
				}

				for(i=0; i<N+1; i++){
					b[i*(N+1)+i] = 14./15.;
					if(i<N) b[i*(N+1)+i+1] = 2./45.;
					if(i>0) b[i*(N+1)+i-1] = 2./45.;
					if(i>1) b[i*(N+1)+i-2] = -1./90.;
					if(i<N-1) b[i*(N+1)+i+2] = -1./90.;
				}
			break;
		case 8:
			// Numerov method for TISE1D 8-th order
			for(i=0; i<N+1; i++){
				a[i*(N+1)+i] = -49./18./h/h/gamma;
				if(i<N) a[i*(N+1)+i+1] = 1.5/h/h/gamma;
				if(i>0) a[i*(N+1)+i-1] = 1.5/h/h/gamma;
				if(i>1) a[i*(N+1)+i-2] = -3./20./h/h/gamma;
				if(i<N-1) a[i*(N+1)+i+2] = -3./20./h/h/gamma;
				if(i>2) a[i*(N+1)+i-3] = 1./90./h/h/gamma;
				if(i<N-2) a[i*(N+1)+i+3] = 1./90./h/h/gamma;
			}

			for(i=0; i<N+1; i++){
				b[i*(N+1)+i] = 54./56.;
				if(i<N) b[i*(N+1)+i+1] = 15./560.;
				if(i>0) b[i*(N+1)+i-1] = 15./560.;
				if(i>1) b[i*(N+1)+i-2] = -6./560.;
				if(i<N-1) b[i*(N+1)+i+2] = -6./560.;
				if(i>2) b[i*(N+1)+i-3] = 1./560.;
				if(i<N-2) b[i*(N+1)+i+3] = 1./560.;
			}
			break;
		case 10:
			// Numerov method for TISE1D 10-th order
			for(i=0; i<N+1; i++){
				a[i*(N+1)+i] = -14350./5040./h/h/gamma;
				if(i<N) a[i*(N+1)+i+1] = 8064./5040./h/h/gamma;
				if(i>0) a[i*(N+1)+i-1] = 8064./5040./h/h/gamma;
				if(i>1) a[i*(N+1)+i-2] = -1008./5040./h/h/gamma;
				if(i<N-1) a[i*(N+1)+i+2] = -1008./5040./h/h/gamma;
				if(i>2) a[i*(N+1)+i-3] = 128./5040./h/h/gamma;
				if(i<N-2) a[i*(N+1)+i+3] = 128./5040./h/h/gamma;
				if(i>3) a[i*(N+1)+i-4] = -9./5040./h/h/gamma;
				if(i<N-3) a[i*(N+1)+i+4] = -9./5040./h/h/gamma;
			}

			for(i=0; i<N+1; i++){
				b[i*(N+1)+i] = 1. - 70./3150;
				if(i<N) b[i*(N+1)+i+1] = 56./3150.;
				if(i>0) b[i*(N+1)+i-1] = 56./3150.;
				if(i>1) b[i*(N+1)+i-2] = -28./3150.;
				if(i<N-1) b[i*(N+1)+i+2] = -28./3150.;
				if(i>2) b[i*(N+1)+i-3] = 8./3150.;
				if(i<N-2) b[i*(N+1)+i+3] = 8./3150.;
				if(i>3) b[i*(N+1)+i-4] = -1./3150.;
				if(i<N-3) b[i*(N+1)+i+4] = -1./3150.;
			}
			break;
		default:
			// Numerov method for TISE1D 4-th order
			for(i=0; i<N+1; i++){
				a[i*(N+1)+i] = -2./h/h/gamma;
				if(i<N) a[i*(N+1)+i+1] = 1./h/h/gamma;
				if(i>0) a[i*(N+1)+i-1] = 1./h/h/gamma;
			}

			for(i=0; i<N+1; i++){
				b[i*(N+1)+i] = 10./12.;
				if(i<N) b[i*(N+1)+i+1] = 1./12;
				if(i>0) b[i*(N+1)+i-1] = 1./12;
			}
			break;
	}

	for(i=0; i<N+1; i++){
		//c[i*(N+1)+i] = V(x1) + 1 / (478.450668 * 0.948204436)*J*(J+1);
		//x1 += h;
		c[i*(N+1)+i] = V[i];
	}

	printf("\n\tMatrixes A and B were filled for %d seconds.\n\n", time(NULL) - t_buf);

	t_buf = time(NULL);                              // Calculating of current time

	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N+1, N+1, b, N+1, ipiv);

	LAPACKE_dgetri(LAPACK_ROW_MAJOR, N+1, b, N+1, ipiv);
	
	printf("\n\tMatrix was inverted for %d seconds.\n\n", time(NULL) - t_buf);

	printf("\n\tMatrix multiplication.\n\n");
	t_buf = time(NULL);                              // Calculating of current time

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N+1, N+1, N+1, -1.0, b, N+1, a, N+1, 1.0, c, N+1);
	
	printf("\n\tIt takes %d seconds.\n\n", time(NULL) - t_buf);

	t_buf = time(NULL);
	printf("\n\tCalculating Eigenvalues and Eigenvectors.\n\n");

	//LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', iElim, c, N+1, w);
	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'N', 'U', N+1, c, N+1, w);

	printf("\n\tEigenvalues were found for %d seconds.\n\n", time(NULL)-t_buf);
	
	for(i=0; i<N+1; i++){
		//if((w[i]>Emin)&&(w[i]<Elim))
		if(w[i]<V[N])
			printf("v=%d\t%.10f\tError=%.4f\n", i, w[i]*8.064877E+3, fabs(w[i]-i-0.5)/w[i]);
	}
	//double x2 = x1;
	//for(i=0; i<N+1; i++){
	//	printf("%f", x2);
	//	for(j=0; j<N+1; j++)
	//		//if((w[j]>Emin)&&(w[j]<Elim))
	//		if(w[j]<V[N])
	//			printf("\t%e", c[i*(N+1)+j]);
	//	printf("\n");
	//	x2 += h;
	//}

	free(ipiv);
	free(w);
	free(a);
	free(b);
	free(c);

	return 0;
}
