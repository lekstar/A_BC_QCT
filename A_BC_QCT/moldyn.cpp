#include "moldyn.h"
#include "PES.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

int H_eq_A_BC(double* Q, double* P, double* k, double* K, double* V) {

    double R[3] = { 0., 0., 0. }, grad[3];
    int i;

    for (i = 0; i < 3; i++) {
        R[0] += (mass[2] / (mass[1] + mass[2]) * Q[i] + Q[3 + i]) * (mass[2] / (mass[1] + mass[2]) * Q[i] + Q[3 + i]);
        R[2] += Q[i] * Q[i];
        R[1] += (mass[1] / (mass[1] + mass[2]) * Q[i] - Q[3 + i]) * (mass[1] / (mass[1] + mass[2]) * Q[i] - Q[3 + i]);
    }

    for (i = 0; i < 3; i++) {
        R[i] = sqrt(R[i]);
        //		printf("%e\n", R[i]);
    }
    //	printf("%e\n", U_LEPS_dist(R, dim));

    for (i = 0; i < 3; i++)
        k[i] = P[i] * (mass[1] + mass[2]) / (mass[1] * mass[2]);

    for (i = 3; i < 6; i++)
        k[i] = P[i] * (mass[0] + mass[1] + mass[2]) / (mass[0] * (mass[1] + mass[2]));

    for (i = 6; i < 9; i++)
        k[i] = P[i] / (mass[0] + mass[1] + mass[2]);

    *V = U_grad_dist(R, grad);

    for (i = 0; i < 3; i++) {
        k[i + 9] = -1. / R[0] * mass[2] / (mass[1] + mass[2]) * (mass[2] / (mass[1] + mass[2]) * Q[i] + Q[3 + i]) * grad[0] * eVtoE - Q[i] / R[2] * grad[2] * eVtoE -
            1. / R[1] * mass[1] / (mass[1] + mass[2]) * (mass[1] / (mass[1] + mass[2]) * Q[i] - Q[3 + i]) * grad[1] * eVtoE;

        k[i + 12] = -1. / R[0] * (mass[2] / (mass[1] + mass[2]) * Q[i] + Q[3 + i]) * grad[0] * eVtoE + 1. / R[1] * (mass[1] / (mass[1] + mass[2]) * Q[i] - Q[3 + i]) * grad[1] * eVtoE;

        k[15 + i] = 0.;
    }

    *K = 0.;

    for (i = 0; i < 3; i++)
        *K += P[i] * P[i] * (mass[1] + mass[2]) / (mass[1] * mass[2]) + P[i + 3] * P[i + 3] * (mass[0] + mass[1] + mass[2]) / (mass[0] * (mass[1] + mass[2])) + P[i + 6] * P[i + 6] / (mass[0] * (mass[1] + mass[2]));

    *K = *K * 0.5 / eVtoE;

    return 0;
}

int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps) {
	double time, R1, r[18], k1[18], k2[18], k3[18], k4[18], k5[18], bias[18], Q[9], P[9], K, V, H, Ri[3], * data, ro, Ev[3], Q0[9], P0[9];
	int i, status = -1, nstep, dim, vJ[2];

	dim = 1 + 2 * 3 * N_ATOM + 2;

	time = 0.;
	nstep = 1;
	ro = 6.0;

	for (i = 0; i < 3 * N_ATOM; i++) {
		Q0[i] = Q_start[i];
		P0[i] = P_start[i];
	}

	data = (double*)calloc(nstep * dim, sizeof(double)); 		// array of time, coordinates and momentum

	//	step = 0.01;

	for (i = 0; i < 3; i++) {
		// Coordinates
		r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
		r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
		r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

		// Momentum
		r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
		r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
		r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
	}

	H_eq_A_BC(Q0, P0, k1, &K, &V);
	H = K + V;

	data[0] = time;
	for (i = 1; i <= 18; i++)
		data[i] = r[i - 1];
	data[19] = K;
	data[20] = V;

	for (; time <= t_fin; ) {

		H_eq_A_BC(Q0, P0, k1, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k1[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + k1[i];
			P[i] = P0[i] + k1[i + 3 * N_ATOM];
		}
		H_eq_A_BC(Q, P, k2, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k2[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.5 * (k1[i] + k2[i]);
			P[i] = P0[i] + 0.5 * (k1[i + 3 * N_ATOM] + k2[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k3, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k3[i] *= step;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.375 * (k1[i] + k3[i]);
			P[i] = P0[i] + 0.375 * (k1[i + 3 * N_ATOM] + k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k4, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k4[i] = k1[i] + 1.333333333333333 * step * k4[i];

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 1.5 * (k4[i] - k3[i]);
			P[i] = P0[i] + 1.5 * (k4[i + 3 * N_ATOM] - k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k5, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k5[i] *= step * 0.333333333333333333;

		for (i = 0; i < 6 * N_ATOM; i++) {
			bias[i] = 0.2 * k4[i] - 0.3 * k3[i] - 0.1 * k5[i];
		}

		R1 = bias[0];

		for (i = 1; i < 6 * N_ATOM; i++)
			if (fabs(R1) < fabs(bias[i]))
				R1 = bias[i];

		H_eq_A_BC(Q0, P0, k1, &K, &V);

		if ((fabs(R1) > eps) || (fabs(H - K - V) > 1e-5)) {
			step /= 1.05;
		}
		else {
			time += step;
			nstep++;
			for (i = 0; i < 3 * N_ATOM; i++) {
				Q0[i] += 0.5 * (k4[i] + k5[i]);
				P0[i] += 0.5 * (k4[i + 3 * N_ATOM] + k5[i + 3 * N_ATOM]);
			}

			H_eq_A_BC(Q0, P0, k1, &K, &V);

			//            printf("Process %d\t%f\t%f\t%f\t%f\n", omp_get_thread_num(), time, step, K, V);

			for (i = 0; i < 3; i++) {
				// Coordinates
				r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
				r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
				r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

				// Momentum
				r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
				r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
				r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
			}

			data = (double*)realloc((void*)data, nstep * dim * sizeof(double));

			data[(nstep - 1) * dim] = time;
			for (i = 1; i <= 18; i++)
				data[(nstep - 1) * dim + i] = r[i - 1];
			data[(nstep - 1) * dim + 19] = K;
			data[(nstep - 1) * dim + 20] = V;

			for (i = 0; i < 3; i++)
				Ri[i] = 0;

			for (i = 0; i < 3; i++) {
				Ri[0] += (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]) * (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]);
				Ri[2] += Q0[i] * Q0[i];
				Ri[1] += (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]) * (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]);
			}

			if ((Ri[0] + Ri[1] + Ri[2] > 3. * (re_AB + re_AC + re_BC))) {

				calc_Ev(data, dim, nstep, Ev, Ri);

				if ((Ev[0] > D_lim) && (Ev[1] > D_lim) && (Ev[2] > D_lim)) {
					status = 4;									// A+B+C products
					if (time > 70)
						break;
				}
				else
					if ((Ev[0] > D_lim) && (Ev[1] > D_lim)) {
						status = 3;						// A+BC products, no reaction

						if ((Ri[0] > ro) && (Ri[1] > ro))
							if (time > 70)
								break;
					}
					else
						if ((Ev[0] > D_lim) && (Ev[2] > D_lim)) {
							status = 2;							// AC+B products

							if ((Ri[0] > ro) && (Ri[2] > ro))
								if (time > 70)
									break;
						}
						else
							if ((Ev[1] > D_lim) && (Ev[2] > D_lim)) {
								status = 1;							// AB+C products

								if ((Ri[1] > ro) && (Ri[2] > ro))
									if (time > 70)
										break;
							}
							else {
								status = 5;						// ABC* product
							}
				//if(time>200){
				//	printf("%f\t%f\t%f\t%f\t%f\t%f\n", Ri[0], Ri[1], Ri[2], Ev[0], Ev[1], Ev[2]);
				//	//getch();
				//}

			}

			if (fabs(R1) < 32. * eps)
				if (step < 0.1)
					step *= 1.05;
		}
	}

	switch (status) {
	case 1: {
		calc_vJ(0, data, dim, nstep, vJ);
		break;
	}
	case 2: {
		calc_vJ(1, data, dim, nstep, vJ);
		break;
	}
	case 3: {
		calc_vJ(2, data, dim, nstep, vJ);
		break;
	}
	}

	free(data);

	//	if (status == -1)
			//printf("%e\t%e\t%e\n", Ev[0], Ev[1], Ev[2]);
	//		getch();

	return status;
}

int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps, int* react) {
    double time, R1, r[18], k1[18], k2[18], k3[18], k4[18], k5[18], bias[18], Q[9], P[9], K, V, H, Ri[3], *data, ro, Ev[3], Q0[9], P0[9], t_lim1;
    int i, status = -1, nstep, dim, vJ[2];

    dim = 1 + 2 * 3 * N_ATOM + 2;

    time = 0.;
    nstep = 1;
    ro = 6.0;
	t_lim1 = 2.*ro/( P_start[5]/( mass[0]*(mass[1]*mass[2])/(mass[0] + mass[1] + mass[2]) ) );

    for (i=0; i<3*N_ATOM; i++) {
        Q0[i] = Q_start[i];
        P0[i] = P_start[i];
    }

    data = (double*)calloc(nstep * dim, sizeof(double)); 		// array of time, coordinates and momentum

    //	step = 0.01;

    for (i = 0; i < 3; i++) {
        // Coordinates
        r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
        r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
        r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

        // Momentum
        r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
        r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
        r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
    }

    H_eq_A_BC(Q0, P0, k1, &K, &V);
    H = K + V;

    data[0] = time;
    for (i = 1; i <= 18; i++)
        data[i] = r[i - 1];
    data[19] = K;
    data[20] = V;

    for (; time <= t_fin; ) {

        H_eq_A_BC(Q0, P0, k1, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k1[i] *= step * 0.333333333333333333;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + k1[i];
            P[i] = P0[i] + k1[i + 3 * N_ATOM];
        }
        H_eq_A_BC(Q, P, k2, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k2[i] *= step * 0.333333333333333333;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 0.5 * (k1[i] + k2[i]);
            P[i] = P0[i] + 0.5 * (k1[i + 3 * N_ATOM] + k2[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k3, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k3[i] *= step;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 0.375 * (k1[i] + k3[i]);
            P[i] = P0[i] + 0.375 * (k1[i + 3 * N_ATOM] + k3[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k4, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k4[i] = k1[i] + 1.333333333333333 * step * k4[i];

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 1.5 * (k4[i] - k3[i]);
            P[i] = P0[i] + 1.5 * (k4[i + 3 * N_ATOM] - k3[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k5, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k5[i] *= step * 0.333333333333333333;

        for (i = 0; i < 6 * N_ATOM; i++) {
            bias[i] = 0.2 * k4[i] - 0.3 * k3[i] - 0.1 * k5[i];
        }

        R1 = bias[0];

        for (i = 1; i < 6 * N_ATOM; i++)
            if (fabs(R1) < fabs(bias[i]))
                R1 = bias[i];

        H_eq_A_BC(Q0, P0, k1, &K, &V);

        if ((fabs(R1) > eps) || (fabs(H - K - V) > 1e-5)) {
            step /= 1.05;
        }
        else {
            time += step;
            nstep++;
            for (i = 0; i < 3 * N_ATOM; i++) {
                Q0[i] += 0.5 * (k4[i] + k5[i]);
                P0[i] += 0.5 * (k4[i + 3 * N_ATOM] + k5[i + 3 * N_ATOM]);
            }

            H_eq_A_BC(Q0, P0, k1, &K, &V);

            //            printf("Process %d\t%f\t%f\t%f\t%f\n", omp_get_thread_num(), time, step, K, V);

            for (i = 0; i < 3; i++) {
                // Coordinates
                r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
                r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
                r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

                // Momentum
                r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
                r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
                r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
            }

            data = (double*)realloc((void*)data, nstep * dim * sizeof(double));

            data[(nstep - 1) * dim] = time;
            for (i = 1; i <= 18; i++)
                data[(nstep - 1) * dim + i] = r[i - 1];
            data[(nstep - 1) * dim + 19] = K;
            data[(nstep - 1) * dim + 20] = V;

			for(i=0; i<3; i++)
				Ri[i] = 0;

			for (i = 0; i < 3; i++) {
				Ri[0] += (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]) * (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]);
				Ri[2] += Q0[i] * Q0[i];
				Ri[1] += (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]) * (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]);
			}

            if ((Ri[0] + Ri[1] + Ri[2] > 3. * (re_AB + re_AC + re_BC))) {

				calc_Ev(data, dim, nstep, Ev, Ri);

                if ((Ev[0] > D_lim) && (Ev[1] > D_lim) && (Ev[2] > D_lim)) {
                    status = 4;									// A+B+C products
					if(time > t_lim1)
						break;
                }
                else
                    if ((Ev[0] > D_lim) && (Ev[1] > D_lim)) {
                        status = 3;						// A+BC products, no reaction

                        if ((Ri[0] > ro) && (Ri[1] > ro))
							if (time > t_lim1)
								break;
                    }
                    else
                        if ((Ev[0] > D_lim) && (Ev[2] > D_lim)) {
                            status = 2;							// AC+B products

                            if ((Ri[0] > ro) && (Ri[2] > ro))
								if (time > t_lim1)
									break;
                        }
                        else
                            if ((Ev[1] > D_lim) && (Ev[2] > D_lim)) {
                                status = 1;							// AB+C products

                                if ((Ri[1] > ro) && (Ri[2] > ro))
									if (time > t_lim1)
										break;
                            }
                            else {
                                status = 5;						// ABC* product
                            }
				//if(time>200){
				//	printf("%f\t%f\t%f\t%f\t%f\t%f\n", Ri[0], Ri[1], Ri[2], Ev[0], Ev[1], Ev[2]);
				//	//getch();
				//}

            }

            if (fabs(R1) < 32. * eps)
                if (step < 0.1)
                    step *= 1.05;
        }
    }

	switch(status){
		case 1: {
				calc_vJ(0, data, dim, nstep, vJ);
				break;
			}
		case 2: {
				calc_vJ(1, data, dim, nstep, vJ);
				break;
			}
		case 3: {
				calc_vJ(2, data, dim, nstep, vJ);
				break;
			}
	}

    free(data);

//	if (status == -1)
		//printf("%e\t%e\t%e\n", Ev[0], Ev[1], Ev[2]);
//		getch();

#pragma omp critical
    {
        if ((status >= 1) && (status <= 3))
            react[status - 1]++;
        else
            if(status>1)
                react[status - 1]++;
    }

    return status;
}

int traject_A_BC(double* Q_start, double* P_start, double step, double t_fin, double eps, int** react) {
    double time, R1, r[18], k1[18], k2[18], k3[18], k4[18], k5[18], bias[18], Q[9], P[9], K, V, H, Ri[3], * data, ro, Ev[3], Q0[9], P0[9], t_lim1;
    int i, status = -1, nstep, dim, vJ[2];

    dim = 1 + 2 * 3 * N_ATOM + 2;

    time = 0.;
    nstep = 1;
    ro = 6.0;
	t_lim1 = 2.*ro/( P_start[5]/( mass[0]*(mass[1]*mass[2])/(mass[0] + mass[1] + mass[2]) ) );

    for (i = 0; i < 3 * N_ATOM; i++) {
        Q0[i] = Q_start[i];
        P0[i] = P_start[i];
    }

    data = (double*)calloc(nstep * dim, sizeof(double)); 		// array of time, coordinates and momentum

    //	step = 0.01;

    for (i = 0; i < 3; i++) {
        // Coordinates
        r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
        r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
        r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

        // Momentum
        r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
        r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
        r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
    }

    H_eq_A_BC(Q0, P0, k1, &K, &V);
    H = K + V;

    data[0] = time;
    for (i = 1; i <= 18; i++)
        data[i] = r[i - 1];
    data[19] = K;
    data[20] = V;

    for (; time <= t_fin; ) {

        H_eq_A_BC(Q0, P0, k1, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k1[i] *= step * 0.333333333333333333;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + k1[i];
            P[i] = P0[i] + k1[i + 3 * N_ATOM];
        }
        H_eq_A_BC(Q, P, k2, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k2[i] *= step * 0.333333333333333333;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 0.5 * (k1[i] + k2[i]);
            P[i] = P0[i] + 0.5 * (k1[i + 3 * N_ATOM] + k2[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k3, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k3[i] *= step;

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 0.375 * (k1[i] + k3[i]);
            P[i] = P0[i] + 0.375 * (k1[i + 3 * N_ATOM] + k3[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k4, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k4[i] = k1[i] + 1.333333333333333 * step * k4[i];

        for (i = 0; i < 3 * N_ATOM; i++) {
            Q[i] = Q0[i] + 1.5 * (k4[i] - k3[i]);
            P[i] = P0[i] + 1.5 * (k4[i + 3 * N_ATOM] - k3[i + 3 * N_ATOM]);
        }
        H_eq_A_BC(Q, P, k5, &K, &V);
        for (i = 0; i < 6 * N_ATOM; i++)
            k5[i] *= step * 0.333333333333333333;

        for (i = 0; i < 6 * N_ATOM; i++) {
            bias[i] = 0.2 * k4[i] - 0.3 * k3[i] - 0.1 * k5[i];
        }

        R1 = bias[0];

        for (i = 1; i < 6 * N_ATOM; i++)
            if (fabs(R1) < fabs(bias[i]))
                R1 = bias[i];

        H_eq_A_BC(Q0, P0, k1, &K, &V);

        if ((fabs(R1) > eps) || (fabs(H - K - V) > 1e-5)) {
            step /= 1.05;
        }
        else {
            time += step;
            nstep++;
            for (i = 0; i < 3 * N_ATOM; i++) {
                Q0[i] += 0.5 * (k4[i] + k5[i]);
                P0[i] += 0.5 * (k4[i + 3 * N_ATOM] + k5[i + 3 * N_ATOM]);
            }

            H_eq_A_BC(Q0, P0, k1, &K, &V);

            //            printf("Process %d\t%f\t%f\t%f\t%f\n", omp_get_thread_num(), time, step, K, V);

            for (i = 0; i < 3; i++) {
                // Coordinates
                r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
                r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
                r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

                // Momentum
                r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
                r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
                r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
            }

            data = (double*)realloc((void*)data, nstep * dim * sizeof(double));

            data[(nstep - 1) * dim] = time;
            for (i = 1; i <= 18; i++)
                data[(nstep - 1) * dim + i] = r[i - 1];
            data[(nstep - 1) * dim + 19] = K;
            data[(nstep - 1) * dim + 20] = V;

            for (i = 0; i < 3; i++)
                Ri[i] = 0;

            for (i = 0; i < 3; i++) {
                Ri[0] += (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]) * (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]);
                Ri[2] += Q0[i] * Q0[i];
                Ri[1] += (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]) * (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]);
            }

            if ((Ri[0] + Ri[1] + Ri[2] > 3. * (re_AB + re_AC + re_BC))) {

                calc_Ev(data, dim, nstep, Ev, Ri);

                if ((Ev[0] > D_lim) && (Ev[1] > D_lim) && (Ev[2] > D_lim)) {
                    status = 4;									// A+B+C products
                    if (time > t_lim1)
                        break;
                }
                else
                    if ((Ev[0] > D_lim) && (Ev[1] > D_lim)) {
                        status = 3;						// A+BC products, no reaction

                        if ((Ri[0] > ro) && (Ri[1] > ro))
                            if (time > t_lim1)
                                break;
                    }
                    else
                        if ((Ev[0] > D_lim) && (Ev[2] > D_lim)) {
                            status = 2;							// AC+B products

                            if ((Ri[0] > ro) && (Ri[2] > ro))
                                if (time > t_lim1)
                                    break;
                        }
                        else
                            if ((Ev[1] > D_lim) && (Ev[2] > D_lim)) {
                                status = 1;							// AB+C products

                                if ((Ri[1] > ro) && (Ri[2] > ro))
                                    if (time > t_lim1)
                                        break;
                            }
                            else {
                                status = 5;						// ABC* product
                            }
                //if(time>200){
                //	printf("%f\t%f\t%f\t%f\t%f\t%f\n", Ri[0], Ri[1], Ri[2], Ev[0], Ev[1], Ev[2]);
                //	//getch();
                //}

            }

            if (fabs(R1) < 32. * eps)
                if (step < 0.1)
                    step *= 1.05;
        }
    }

    switch (status) {
    case 1: {
        calc_vJ(0, data, dim, nstep, vJ);
        break;
    }
    case 2: {
        calc_vJ(1, data, dim, nstep, vJ);
        break;
    }
    case 3: {
        calc_vJ(2, data, dim, nstep, vJ);
        break;
    }
    }

    free(data);

    //	if (status == -1)
            //printf("%e\t%e\t%e\n", Ev[0], Ev[1], Ev[2]);
    //		getch();
#pragma omp critical
    {
        if ((status >= 1) && (status <= 3))
            react[status - 1][vJ[0]]++;
        else
            if (status > 1)
                react[status - 1][0]++;
    }
    return status;
}

int traject_A_BC_print(double* Q_start, double* P_start, double step, double t_fin, double eps) {
	double time, R1, r[18], k1[18], k2[18], k3[18], k4[18], k5[18], bias[18], Q[9], P[9], K, V, H, Ri[3], * data, ro, Ev[3], Q0[9], P0[9], t_lim1;
	int i, status = -1, nstep, dim, vJ[2];
	FILE* drc;

	dim = 1 + 2 * 3 * N_ATOM + 2;

	time = 0.;
	nstep = 1;
	ro = 8.0;
	t_lim1 = 2.*ro/( P_start[5]/( mass[0]*(mass[1]*mass[2])/(mass[0] + mass[1] + mass[2]) ) );

	for (i = 0; i < 3 * N_ATOM; i++) {
		Q0[i] = Q_start[i];
		P0[i] = P_start[i];
	}

	data = (double*)calloc(nstep * dim, sizeof(double)); 		// array of time, coordinates and momentum

	//	step = 0.01;

	for (i = 0; i < 3; i++) {
		// Coordinates
		r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
		r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
		r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

		// Momentum
		r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
		r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
		r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
	}

	H_eq_A_BC(Q0, P0, k1, &K, &V);
	H = K + V;

	data[0] = time;
	for (i = 1; i <= 18; i++)
		data[i] = r[i - 1];
	data[19] = K;
	data[20] = V;

	for (; time <= t_fin; ) {

		H_eq_A_BC(Q0, P0, k1, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k1[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + k1[i];
			P[i] = P0[i] + k1[i + 3 * N_ATOM];
		}
		H_eq_A_BC(Q, P, k2, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k2[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.5 * (k1[i] + k2[i]);
			P[i] = P0[i] + 0.5 * (k1[i + 3 * N_ATOM] + k2[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k3, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k3[i] *= step;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.375 * (k1[i] + k3[i]);
			P[i] = P0[i] + 0.375 * (k1[i + 3 * N_ATOM] + k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k4, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k4[i] = k1[i] + 1.333333333333333 * step * k4[i];

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 1.5 * (k4[i] - k3[i]);
			P[i] = P0[i] + 1.5 * (k4[i + 3 * N_ATOM] - k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k5, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k5[i] *= step * 0.333333333333333333;

		for (i = 0; i < 6 * N_ATOM; i++) {
			bias[i] = 0.2 * k4[i] - 0.3 * k3[i] - 0.1 * k5[i];
		}

		R1 = bias[0];

		for (i = 1; i < 6 * N_ATOM; i++)
			if (fabs(R1) < fabs(bias[i]))
				R1 = bias[i];

		H_eq_A_BC(Q0, P0, k1, &K, &V);

		if ((fabs(R1) > eps) || (fabs(H - K - V) > 1e-5)) {
			step /= 1.05;
		}
		else {
			time += step;
			nstep++;
			for (i = 0; i < 3 * N_ATOM; i++) {
				Q0[i] += 0.5 * (k4[i] + k5[i]);
				P0[i] += 0.5 * (k4[i + 3 * N_ATOM] + k5[i + 3 * N_ATOM]);
			}

			H_eq_A_BC(Q0, P0, k1, &K, &V);

			//            printf("Process %d\t%f\t%f\t%f\t%f\n", omp_get_thread_num(), time, step, K, V);

			for (i = 0; i < 3; i++) {
				// Coordinates
				r[i] = (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
				r[i + 3] = -mass[2] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];
				r[i + 6] = mass[1] / (mass[1] + mass[2]) * Q0[i] - mass[0] / (mass[0] + mass[1] + mass[2]) * Q0[i + 3] + Q0[i + 6];

				// Momentum
				r[i + 9] = P0[i + 3] + mass[0] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
				r[i + 12] = -P0[i] - mass[1] / (mass[1] + mass[2]) * P0[i + 3] + mass[1] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
				r[i + 15] = P0[i] - mass[2] / (mass[1] + mass[2]) * P0[i + 3] + mass[2] / (mass[0] + mass[1] + mass[2]) * P0[i + 6];
			}

			data = (double*)realloc((void*)data, nstep * dim * sizeof(double));

			data[(nstep - 1) * dim] = time;
			for (i = 1; i <= 18; i++)
				data[(nstep - 1) * dim + i] = r[i - 1];
			data[(nstep - 1) * dim + 19] = K;
			data[(nstep - 1) * dim + 20] = V;

			for (i = 0; i < 3; i++)
				Ri[i] = 0;

			for (i = 0; i < 3; i++) {
				Ri[0] += (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]) * (mass[2] / (mass[1] + mass[2]) * Q0[i] + Q0[3 + i]);
				Ri[2] += Q0[i] * Q0[i];
				Ri[1] += (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]) * (mass[1] / (mass[1] + mass[2]) * Q0[i] - Q0[3 + i]);
			}

			if ((Ri[0] + Ri[1] + Ri[2] > 3. * (re_AB + re_AC + re_BC))) {

				calc_Ev(data, dim, nstep, Ev, Ri);

				if ((Ev[0] > D_lim) && (Ev[1] > D_lim) && (Ev[2] > D_lim)) {
					status = 4;									// A+B+C products
					if (time > t_lim1)
						break;
				}
				else
					if ((Ev[0] > D_lim) && (Ev[1] > D_lim)) {
						status = 3;						// A+BC products, no reaction

						if ((Ri[0] > ro) && (Ri[1] > ro))
							if (time > t_lim1)
								break;
					}
					else
						if ((Ev[0] > D_lim) && (Ev[2] > D_lim)) {
							status = 2;							// AC+B products

							if ((Ri[0] > ro) && (Ri[2] > ro))
								if (time > t_lim1)
									break;
						}
						else
							if ((Ev[1] > D_lim) && (Ev[2] > D_lim)) {
								status = 1;							// AB+C products

								if ((Ri[1] > ro) && (Ri[2] > ro))
									if (time > t_lim1)
										break;
							}
							else {
								status = 5;						// ABC* product
							}
				//if(time>200){
				//	printf("%f\t%f\t%f\t%f\t%f\t%f\n", Ri[0], Ri[1], Ri[2], Ev[0], Ev[1], Ev[2]);
				//	//getch();
				//}

			}

			if (fabs(R1) < 32. * eps)
				if (step < 0.1)
					step *= 1.05;
		}
	}

	switch (status) {
	case 1: {
		calc_vJ(0, data, dim, nstep, vJ);
		break;
	}
	case 2: {
		calc_vJ(1, data, dim, nstep, vJ);
		break;
	}
	case 3: {
		calc_vJ(2, data, dim, nstep, vJ);
		break;
	}
	}

	drc = fopen("traject.drc", "w");

	fprintf(drc, "### PREAMBLE AND COMMENTS ###\n### EPILOGUE ###\n### ELEMENTAL COMPOSITION ###\n");
	for(i=0; i<3; i++)
		fprintf(drc, "AN%d:\t%d\n", i+1, (int)Ze[i]);
	fprintf(drc, "###   TABLE OF RESULTS    ###\n#	t, fs	Epot, eV\n\tX, A	Y, A	Z, A	Vx, A/fs	Vy, A/fs	Vz, A/fs\n\n");

	for(i=0; i<nstep; i++){
		fprintf(drc, "%d\t%.8f\t%.8f\n", i, data[i * dim], data[i * dim + 20]);
		for(int j=0; j<3; j++){
			fprintf(drc, "\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t%.9E\t\n", data[i*dim+1+3*j], data[i*dim+2+3*j], data[i*dim+3+3*j], data[i*dim+10+3*j]/mass[j], data[i*dim+11+3*j]/mass[j], data[i*dim+12+3*j]/mass[j]);
		}
		fprintf(drc, "\n");
	}

	fclose(drc);

	free(data);

	//	if (status == -1)
			//printf("%e\t%e\t%e\n", Ev[0], Ev[1], Ev[2]);
	//		getch();

	return status;
}


int create_init_coord(double *Q, double *P, double *x, double *Rmp0, double b, double Tv, double Vr, int J) {

	double phi, eta, theta, Rmp, P0;

	phi = 2. * M_PI * x[0];
	eta = 2. * M_PI * x[1];
	theta = acos(2. * x[2] - 1.);

	if (x[3] < 0.5)
		Rmp = Rmp0[0];
	else
		Rmp = Rmp0[1];

	Q[3] = 0.;
	Q[4] = b*sqrt(x[5]);
	Q[5] = -sqrt(100. - Q[4] * Q[4]) - 0.5 * Tv * Vr * x[4];

	P[3] = 0;
	P[4] = 0;
	P[5] = mass[0] * (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]) * (Vr);		// reduced mass for A, BC and Vr in A/fs, usual less 0.2 A/fs

	Q[0] = Rmp * sin(theta) * cos(phi);
	Q[1] = Rmp * sin(theta) * sin(phi);
	Q[2] = Rmp * (2. * x[2] - 1.);

	P0 = sqrt(J * (J + 1)) * 6.35077993E-03 / Rmp;

	P[0] = -P0 * (sin(phi) * cos(eta) + cos(phi) * cos(theta) * sin(eta));
	P[1] = P0 * (cos(phi) * cos(eta) - sin(phi) * cos(theta) * sin(eta));
	P[2] = P0 * sin(theta) * sin(eta);

	return 0;
}

int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps) {

	int state;
	double P[9], Q[9];

	create_init_coord(Q, P, x, Rmp, b, Tv, Vr, J);

	state = traject_A_BC(Q, P, 0.04, 3000., 1e-9);

	return state;
}

int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps, int* react) {

	int state;
	double P[9], Q[9];

	create_init_coord(Q, P, x, Rmp, b, Tv, Vr, J);

	state = traject_A_BC(Q, P, 0.04, 3000., 1e-9, react);

	return state;
}

int traject_A_BC(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps, int** react) {

	int state;
	double P[9], Q[9];

	create_init_coord(Q, P, x, Rmp, b, Tv, Vr, J);

	state = traject_A_BC(Q, P, 0.04, 3000., 1e-9, react);

	return state;
}

int traject_A_BC_print(double* x, double* Rmp, double b, double Tv, double Vr, int J, double step, double t_fin, double eps) {

	int state;
	double P[9], Q[9];

	create_init_coord(Q, P, x, Rmp, b, Tv, Vr, J);

	state = traject_A_BC_print(Q, P, 0.04, 3000., 1e-9);

	return state;
}

int calc_Ev(double* data, int dim, int nstep, double* Ev, double* R) {
	double scalRP, modR, modP, modM, J, mu, V, R1[3];
	int i, j;

	i = nstep - 1;

	// Inernal energy of BC molecule
	scalRP = 0.;
	modR = 0.;
	modP = 0.;
	mu = mass[1] * mass[2] / (mass[1] + mass[2]);
	for (j = 0; j < 3; j++) {
		scalRP += (data[i * dim + 7 + j] - data[i * dim + 4 + j]) * (-mass[2] / (mass[1] + mass[2]) * data[i * dim + 13 + j] + mass[1] / (mass[1] + mass[2]) * data[i * dim + 16 + j]);
		modR += sqr(data[i * dim + 7 + j] - data[i * dim + 4 + j]);
		modP += sqr(-mass[2] / (mass[1] + mass[2]) * data[i * dim + 13 + j] + mass[1] / (mass[1] + mass[2]) * data[i * dim + 16 + j]);
	}

	R1[2] = sqrt(modR);
	R1[0] = 20.;
	R1[1] = 20.;

	R[2] = R1[2];

	modM = modR * modP - sqr(scalRP);
	J = sqrt(0.25 + modM / sqr(0.00635077993)) - 0.5;

	V = U_dist(R1);

	Ev[2] = 0.5 * scalRP * scalRP / modR / mu / 0.00964853326 + V +0.5 * J * (J + 1) * 0.004180159286 / mu / modR;

	// Inernal energy of AB molecule
	scalRP = 0.;
	modR = 0.;
	modP = 0.;
	mu = mass[0] * mass[1] / (mass[0] + mass[1]);
	for (j = 0; j < 3; j++) {
		scalRP += (data[i * dim + 4 + j] - data[i * dim + 1 + j]) * (-mass[1] / (mass[0] + mass[1]) * data[i * dim + 10 + j] + mass[0] / (mass[0] + mass[1]) * data[i * dim + 13 + j]);
		modR += sqr(data[i * dim + 4 + j] - data[i * dim + 1 + j]);
		modP += sqr(-mass[1] / (mass[0] + mass[1]) * data[i * dim + 10 + j] + mass[0] / (mass[0] + mass[1]) * data[i * dim + 13 + j]);
	}

	R1[0] = sqrt(modR);
	R1[1] = 20.;
	R1[2] = 20.;

	R[0] = R1[0];

	modM = modR * modP - sqr(scalRP);
	J = sqrt(0.25 + modM / sqr(0.00635077993)) - 0.5;

	V = U_dist(R1);

	Ev[0] = 0.5 * scalRP * scalRP / modR / mu / 0.00964853326 + V +0.5 * J * (J + 1) * 0.004180159286 / mu / modR;

	// Inernal energy of AC molecule
	scalRP = 0.;
	modR = 0.;
	modP = 0.;
	mu = mass[0] * mass[2] / (mass[0] + mass[2]);
	for (j = 0; j < 3; j++) {
		scalRP += (data[i * dim + 7 + j] - data[i * dim + 1 + j]) * (-mass[2] / (mass[0] + mass[2]) * data[i * dim + 10 + j] + mass[0] / (mass[0] + mass[2]) * data[i * dim + 16 + j]);
		modR += sqr(data[i * dim + 7 + j] - data[i * dim + 1 + j]);
		modP += sqr(-mass[2] / (mass[0] + mass[2]) * data[i * dim + 10 + j] + mass[0] / (mass[0] + mass[2]) * data[i * dim + 16 + j]);
	}

	R1[1] = sqrt(modR);
	R1[0] = 20.;
	R1[2] = 20.;

	R[1] = R1[1];

	modM = modR * modP - sqr(scalRP);
	J = sqrt(0.25 + modM / sqr(0.00635077993)) - 0.5;

	V = U_dist(R1);

	Ev[1] = 0.5 * scalRP * scalRP / modR / mu / 0.00964853326 + V +0.5 * J * (J + 1) * 0.004180159286 / mu / modR;

//		printf("%.8f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", data[i*dim], Ev[0], Ev[1], Ev[2], R[0], R[1], R[2]);

	return 0;
}

int calc_vJ(int n, double* data, int dim, int nstep, int* vJ) {

	int i, j, mi, mj, qi, pi, qj, pj;
	double scalRP, modR, modP, modM, mu, Ek[3], J = 0, Te;

	Ek[0] = 1e10;
	Ek[1] = 1e10;
	Ek[2] = 1e10;

	switch (n) {
	case 0: {
		qi = 1;
		qj = 4;
		pi = 10;
		pj = 13;
		mi = 0;
		mj = 1;
		Te = Te_AB;
		break;
	}
	case 1: {
		qi = 1;
		qj = 7;
		pi = 10;
		pj = 16;
		mi = 0;
		mj = 2;
		Te = Te_AC;
		break;
	}
	case 2: {
		qi = 4;
		qj = 7;
		pi = 13;
		pj = 16;
		mi = 1;
		mj = 2;
		Te = Te_BC;
		break;
	}
	default: {
		return -1;
		break;
	}
	}

	for (i = nstep - 1; i >= 0; i--) {
		// Inernal energy of AC molecule
		scalRP = 0.;
		modR = 0.;
		modP = 0.;
		mu = mass[mi] * mass[mj] / (mass[mi] + mass[mj]);
		for (j = 0; j < 3; j++) {
			scalRP += (data[i * dim + qj + j] - data[i * dim + qi + j]) * (-mass[mj] / (mass[mi] + mass[mj]) * data[i * dim + pi + j] + mass[mi] / (mass[mi] + mass[mj]) * data[i * dim + pj + j]);
			modR += sqr(data[i * dim + qj + j] - data[i * dim + qi + j]);
			modP += sqr(-mass[mj] / (mass[mi] + mass[mj]) * data[i * dim + pi + j] + mass[mi] / (mass[mi] + mass[mj]) * data[i * dim + pj + j]);
		}

		modM = modR * modP - sqr(scalRP);
		if (i >= nstep - 2)
			J = sqrt(0.25 + modM / sqr(0.00635077993)) - 0.5;

		Ek[0] = Ek[1];
		Ek[1] = Ek[2];
		Ek[2] = 0.5 * scalRP * scalRP / modR / mu / 0.00964853326;

		if ((Ek[0] < Ek[1]) && (Ek[1] > Ek[2]))
			break;
	}

	Ek[1] += Te;

	switch (n) {
	case 0: {
		if (Ek[1] < Ev_AB[0]) {
			vJ[0] = 0;
			break;
		}

		for (i = 0; i < dim_Ev_AB - 1; i++)
			if ((Ek[1] >= Ev_AB[i]) && (Ek[1] < Ev_AB[i + 1]))
				break;

		if (Ek[1] - Ev_AB[i] >= Ev_AB[i + 1] - Ek[1])
			vJ[0] = i + 1;
		else
			vJ[0] = i;

		if (Ek[1] > Ev_AB[dim_Ev_AB - 1])
			vJ[0] = dim_Ev_AB - 1;

		break;
	}
	case 1: {
		if (Ek[1] < Ev_AC[0]) {
			vJ[0] = 0;
			break;
		}

		for (i = 0; i < dim_Ev_AC - 1; i++)
			if ((Ek[1] >= Ev_AC[i]) && (Ek[1] < Ev_AC[i + 1]))
				break;

		if (Ek[1] - Ev_AC[i] >= Ev_AC[i + 1] - Ek[1])
			vJ[0] = i + 1;
		else
			vJ[0] = i;

		if (Ek[1] > Ev_AC[dim_Ev_AC - 1])
			vJ[0] = dim_Ev_AC - 1;

		break;
	}
	case 2: {
		if (Ek[1] < Ev_BC[0]) {
			vJ[0] = 0;
			break;
		}

		for (i = 0; i < dim_Ev_BC - 1; i++)
			if ((Ek[1] >= Ev_BC[i]) && (Ek[1] < Ev_BC[i + 1]))
				break;

		if (Ek[1] - Ev_BC[i] >= Ev_BC[i + 1] - Ek[1])
			vJ[0] = i + 1;
		else
			vJ[0] = i;

		if (Ek[1] > Ev_BC[dim_Ev_BC - 1])
			vJ[0] = dim_Ev_BC - 1;

		break;
	}
	default: {
		return -1;
		break;
	}
	}

	vJ[1] = (int)(J + 0.500000000001);    // Rounding to closest integer
	//	printf("v=%d\tJ=%d\n", vJ[0], vJ[1]);

	return 0;
}

int find_Rmp_BC(int J, double E, double* Rmp) {
	int i, error = 0;

	double ri[3], grad[3], dr, hess[9], re;

	//	E *= 1.239945E-4; // from cm^-1 to eV

	re = re_BC;

	ri[0] = 20.;
	ri[2] = re;
	ri[1] = ri[0] + ri[2];

	for (i = 0; i < 30; i++) {
		U_grad_hess_dist(ri, grad, hess);
		if (fabs(grad[2] - J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2] / ri[2]) < 1e-5) {
			//			printf("re=%.12f\n", ri[2]);
			break;
		}

		ri[2] -= (grad[2] - J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2] / ri[2]) / (hess[8] + 3 * J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2] / ri[2] / ri[2]);
		ri[1] = ri[0] + ri[2];
	}

	re = ri[2];
	ri[2] -= 0.1;
	ri[1] = ri[0] + ri[2];

	for (i = 0; i < 30; i++) {
		double U;
		U = U_grad_dist(ri, grad) + 0.5 * J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2];
		//			printf("%d\t%.9e\t%.9e\n", i+1, U, input[0]);
		dr = (U - E) / (grad[2] - J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2] / ri[2]);
		if (fabs(dr) > 0.2) {
			if (dr < 0)
				dr = -0.2;
			else
				dr = 0.2;
		}
		ri[2] -= dr;
		ri[1] = ri[0] + ri[2];

		if (fabs(U - E) < 1e-6) {
			Rmp[0] = ri[2];
			//			printf("%d\t%.9e\t%.9e\n", i+1, U, ri[2]);
			break;
		}
	}
	if (i >= 20) {
		Rmp[0] = ri[2];
		error++;
	}

	ri[0] = 20.;
	ri[2] = Rmp[0] + fabs(re - Rmp[0]) + 0.1;
	ri[1] = ri[0] + ri[2];

	for (i = 0; i < 30; i++) {
		double U;
		U = U_grad_dist(ri, grad) + 0.5 * J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2];
		//			printf("%d\t%.9e\t%.9e\n", i+1, U, input[0]);
		dr = (U - E) / (grad[2] - J * (J + 1) * 4.180159286E-03 * (mass[1] + mass[2]) / (mass[1] * mass[2]) / ri[2] / ri[2] / ri[2]);
		if (fabs(dr) > 0.2) {
			if (dr < 0)
				dr = -0.2;
			else
				dr = 0.2;
		}
		ri[2] -= dr;
		ri[1] = ri[0] + ri[2];

		if (fabs(U - E) < 1e-6) {
			Rmp[1] = ri[2];
			//printf("%d\t%.9e\t%.9e\n", i+1, U, ri[2]);
			break;
		}
	}

	if (i >= 20) {
		Rmp[1] = ri[2];
		error++;
	}

	return error;
}

double find_Tv_BC(double Rp, int J) {
	double time, R1, k1[18], k2[18], k3[18], k4[18], k5[18], bias[18], Q[9], P[9], K, V, H, max[3], eps = 1e-8, step = 0.05, Tv = 0.;
	int i;

	double P0[9], Q0[9], x[3], phi, eta, theta, Rmp, P1, R_BC;

#pragma omp critical
	{
		for (i = 0; i < 3; i++)
			x[i] = myrand();
	}

	phi = 2. * M_PI * x[0];
	eta = 2. * M_PI * x[1];
	theta = acos(2. * x[2] - 1.);

	Rmp = Rp;

	Q0[3] = 0.;
	Q0[4] = 0.;
	Q0[5] = -sqrt(100. - Q0[4] * Q0[4]);

	P0[3] = 0;
	P0[4] = 0;
	P0[5] = 0;

	Q0[0] = Rmp * sin(theta) * cos(phi);
	Q0[1] = Rmp * sin(theta) * sin(phi);
	Q0[2] = Rmp * (2. * x[2] - 1.);

	P1 = sqrt(J * (J + 1)) * 6.35077993E-03 / Rmp;

	P0[0] = -P1 * (sin(phi) * cos(eta) + cos(phi) * cos(theta) * sin(eta));
	P0[1] = P1 * (cos(phi) * cos(eta) - sin(phi) * cos(theta) * sin(eta));
	P0[2] = P1 * sin(theta) * sin(eta);

	time = 0.;

	H_eq_A_BC(Q0, P0, k1, &K, &V);
	H = K + V;

	R_BC = 0.;
	for (i = 0; i < 3; i++)
		R_BC += sqr(Q0[i]);
	R_BC = sqrt(R_BC);

	for (i = 0; i < 3; i++)
		max[i] = R_BC;

	for (; time <= 100.; ) {

		H_eq_A_BC(Q0, P0, k1, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k1[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + k1[i];
			P[i] = P0[i] + k1[i + 3 * N_ATOM];
		}
		H_eq_A_BC(Q, P, k2, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k2[i] *= step * 0.333333333333333333;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.5 * (k1[i] + k2[i]);
			P[i] = P0[i] + 0.5 * (k1[i + 3 * N_ATOM] + k2[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k3, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k3[i] *= step;

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 0.375 * (k1[i] + k3[i]);
			P[i] = P0[i] + 0.375 * (k1[i + 3 * N_ATOM] + k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k4, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k4[i] = k1[i] + 1.333333333333333 * step * k4[i];

		for (i = 0; i < 3 * N_ATOM; i++) {
			Q[i] = Q0[i] + 1.5 * (k4[i] - k3[i]);
			P[i] = P0[i] + 1.5 * (k4[i + 3 * N_ATOM] - k3[i + 3 * N_ATOM]);
		}
		H_eq_A_BC(Q, P, k5, &K, &V);
		for (i = 0; i < 6 * N_ATOM; i++)
			k5[i] *= step * 0.333333333333333333;

		for (i = 0; i < 6 * N_ATOM; i++) {
			bias[i] = 0.2 * k4[i] - 0.3 * k3[i] - 0.1 * k5[i];
		}

		R1 = bias[0];

		for (i = 1; i < 6 * N_ATOM; i++)
			if (fabs(R1) < fabs(bias[i]))
				R1 = bias[i];

		H_eq_A_BC(Q0, P0, k1, &K, &V);

		if ((fabs(R1) > eps) || (fabs(H - K - V) > 3e-6)) {
			step /= 1.05;
		}
		else {
			time += step;
			for (i = 0; i < 3 * N_ATOM; i++) {
				Q0[i] += 0.5 * (k4[i] + k5[i]);
				P0[i] += 0.5 * (k4[i + 3 * N_ATOM] + k5[i + 3 * N_ATOM]);
			}

			H_eq_A_BC(Q0, P0, k1, &K, &V);

			R_BC = 0.;
			for (i = 0; i < 3; i++)
				R_BC += sqr(Q0[i]);
			R_BC = sqrt(R_BC);

			max[0] = max[1];
			max[1] = max[2];
			max[2] = R_BC;

			if ((max[0] < max[1]) && (max[1] > max[2])) {
				Tv = time - step;
				break;
			}

			if (fabs(R1) < 32. * eps)
				if (step < 0.1)
					step *= 1.05;
		}
	}

	return Tv;
}

double find_Tv_BC(int N, double Rp, int J) {
	int i;
	double sum = 0.;
	if (N <= 0)
		return 0.;
	for (i = 0; i < N; i++) {
		sum += find_Tv_BC(Rp, J);
	}

	return sum / N;
}

double find_Tv_BC_par(int N, double Rp, int J, int num_th) {

	long int i;
	double sum = 0., Tv;
	if (N <= 0)
		return 0.;

#pragma omp parallel for num_threads(num_th) private(Tv, i) schedule(dynamic,3) reduction(+: sum)
	{
		for (i=0; i<N; i++) {

			Tv = find_Tv_BC(Rp, J);
			sum += Tv;
		}
	}

	return sum / N;

}

int P_r_b_fix(int N, double* Rmp, int J, double b, double Vr, double Tv) {

	int i, j, state, **react;

	double x[6];

	react = new int* [5];

	react[0] = new int[dim_Ev_AB];
	react[1] = new int[dim_Ev_AC];
	react[2] = new int[dim_Ev_BC];
	react[3] = new int[1];
	react[4] = new int[1];

	for (i = 0; i < dim_Ev_AB; i++)
		react[0][i] = 0;
	for (i = 0; i < dim_Ev_AC; i++)
		react[1][i] = 0;
	for (i = 0; i < dim_Ev_BC; i++)
		react[2][i] = 0;
	react[3][0] = 0;
	react[4][0] = 0;

	for (j = 0; j < N; j++) {
		for (i = 0; i < 5; i++)
			x[i] = myrand();

		x[5] = 1.;

		state = traject_A_BC(x, Rmp, b, Tv, Vr, J, 0.04, 300., 1e-9, react);

		if (state < 1)
			continue;
	}

	//	printf("\nProducts\n\tA+BC\tAB+C\tAC+B\tA+B+C\tABC*\n");
	//	for(i=0; i<5; i++)
	//		printf("\t%d", react[i]);
	//	printf("\n");

	printf("Number of trajectories:\n\t%d\n", N);

	printf("\n\nReaction A+BC(v\',J\')-->C+AB(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AB; i++)
		printf("%d\t%d\n", i, react[0][i]);

	printf("\n\nReaction A+BC(v\',J\')-->B+AC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AC; i++) {
		printf("%d\t%d\n", i, react[1][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+BC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_BC; i++) {
		printf("%d\t%d\n", i, react[2][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+B+C\n\t%d\n", react[3][0]);

	printf("\n\nReaction A+BC(v\',J\')-->ABC*\n\t%d\n", react[4][0]);

	for (i = 0; i < 5; i++)
		delete[] react[i];
	delete[] react;

	return 0;
}

int P_r_b_fix(int N, double* Rmp, int J, double b, double Vr, double Tv, int* react) {

	int i, j;

	double *x;

	for (i=0; i<5; i++)
		react[i] = 0;

	x = (double*)malloc(N*6 * sizeof(double));
	for(i=0; i<N; i++){
		for(j=0; j<5; j++)
			x[i*6+j] = myrand();
		x[i*6+5] = 1.;
	}

	for (j=0; j<N; j++) {
		traject_A_BC(&x[j*6], Rmp, b, Tv, Vr, J, 0.04, 300., 1e-9, react);
	}

		printf("\nProducts\n\tA+BC\tAB+C\tAC+B\tA+B+C\tABC*\n");
		for(i=0; i<5; i++)
			printf("\t%d", react[i]);
		printf("\n");

	free(x);

	return 0;
}

int P_r_b_fix_par(int N, double* Rmp, int J, double b, double Vr, double Tv, int* react, int num_threads) {

	int i, j, state;

	double *x;

	for(i=0; i<5; i++)
		react[i] = 0;

	x = (double*)malloc(N * 6 * sizeof(double));
	
	for(i=0; i<N; i++) {
		for(j=0; j<5; j++)
			x[i*6+j] = myrand();
		x[i*6+5] = 1.;
	}
#pragma omp parallel for num_threads(num_threads) private(j, state) schedule(dynamic,1)
	for (j = 0; j < N; j++) {
		state = traject_A_BC(&x[j*6], Rmp, b, Tv, Vr, J, 0.04, 3000., 1e-9, react);

		if (state < 1) {
			//printf("Decline trajectory\n");
			continue;
		}
	}

	//printf("\nProducts\n\tA+BC\tAB+C\tAC+B\tA+B+C\tABC*\n");
	//for (i = 0; i < 5; i++)
	//	printf("\t%d", react[i]);
	//printf("\n");

	free(x);

	return 0;
}

double find_b_max_A_BC(int N, double* Rmp, int J, double b_start, double Vr, double Tv) {
	int i, j, state;
	double b, *x, step;

	x = (double*)malloc(N * 6 * sizeof(double));

	step = 0.2;

	for(b=b_start; b<5.0; b+=step){
		
		for(i=0; i<N; i++){
			for(j=0; j<5; j++)
				x[i*6+j] = myrand();
			x[i*6+5] = 1.;
		}

		for(i=0; i<N; i++){

			state = traject_A_BC(&x[i*6], Rmp, b, Tv, Vr, J, 0.04, 3000., 1e-9);

			if(state!=3)
				break;
		}

		if(state==3){
			b -= step;
			step *= 0.5;
		}
		if(step < 0.02)
			break;
	}

	free(x);

	return b;
}

double find_b_max_A_BC_par(int N, double* Rmp, int J, double b_start, double Vr, double Tv, int num_threads) {
	int i, k, react[5];
	double b, step;

	k = N / (3*num_threads);
	N = (k+1) * 3*num_threads;

	step = 0.2;

	for (b = b_start; b < 5.0; b += step) {

		for (i = 0; i < k+1; i++) {
			
			P_r_b_fix_par( 3*num_threads, Rmp, J, b, Vr, Tv, react, num_threads );

			if(react[2] != 3*num_threads)
				break;
		}

		if (react[2] == 3*num_threads) {
			if(b>=step)
				b -= step;
			step *= 0.5;
		}
		if (step < 0.02)
			break;
	}

	if(b>1e-3)
		b += step;

	return b;
}

int P_r_b_max(int N, double* Rmp, int J, double b, double Vr, double Tv) {

	int i, j, state, ** react;

	double x[6];

	react = new int* [5];

	react[0] = new int[dim_Ev_AB];
	react[1] = new int[dim_Ev_AC];
	react[2] = new int[dim_Ev_BC];
	react[3] = new int[1];
	react[4] = new int[1];

	for (i = 0; i < dim_Ev_AB; i++)
		react[0][i] = 0;
	for (i = 0; i < dim_Ev_AC; i++)
		react[1][i] = 0;
	for (i = 0; i < dim_Ev_BC; i++)
		react[2][i] = 0;
	react[3][0] = 0;
	react[4][0] = 0;

	for (j = 0; j < N; j++) {
		for (i = 0; i < 6; i++)
			x[i] = myrand();

		state = traject_A_BC(x, Rmp, b, Tv, Vr, J, 0.04, 300., 1e-9, react);
	}

	printf("Number of trajectories:\n\t%d\n", N);

	printf("\n\nReaction A+BC(v\',J\')-->C+AB(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AB; i++)
		printf("%d\t%d\n", i, react[0][i]);

	printf("\n\nReaction A+BC(v\',J\')-->B+AC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AC; i++) {
		printf("%d\t%d\n", i, react[1][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+BC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_BC; i++) {
		printf("%d\t%d\n", i, react[2][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+B+C\n\t%d\n", react[3][0]);

	printf("\n\nReaction A+BC(v\',J\')-->ABC*\n\t%d\n", react[4][0]);

	for (i = 0; i < 5; i++)
		delete[] react[i];
	delete[] react;

	return 0;
}

int P_r_b_max(int N, double* Rmp, int J, double b, double Vr, double Tv, int **react) {

	int i, j, state;

	double x[6];

	for (i = 0; i < dim_Ev_AB; i++)
		react[0][i] = 0;
	for (i = 0; i < dim_Ev_AC; i++)
		react[1][i] = 0;
	for (i = 0; i < dim_Ev_BC; i++)
		react[2][i] = 0;
	react[3][0] = 0;
	react[4][0] = 0;

	for (j = 0; j < N; j++) {
		for (i = 0; i < 6; i++)
			x[i] = myrand();

		state = traject_A_BC(x, Rmp, b, Tv, Vr, J, 0.04, 300., 1e-9, react);
	}

	printf("Number of trajectories:\n\t%d\n", N);

	printf("\n\nReaction A+BC(v\',J\')-->C+AB(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AB; i++)
		printf("%d\t%d\n", i, react[0][i]);

	printf("\n\nReaction A+BC(v\',J\')-->B+AC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AC; i++) {
		printf("%d\t%d\n", i, react[1][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+BC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_BC; i++) {
		printf("%d\t%d\n", i, react[2][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+B+C\n\t%d\n", react[3][0]);

	printf("\n\nReaction A+BC(v\',J\')-->ABC*\n\t%d\n", react[4][0]);

	return 0;
}

int P_r_b_max_par(int N, double* Rmp, int J, double b, double Vr, double Tv, int num_threads) {

	int i, j, state, ** react;

	double *x;

	react = new int* [5];

	react[0] = new int[dim_Ev_AB];
	react[1] = new int[dim_Ev_AC];
	react[2] = new int[dim_Ev_BC];
	react[3] = new int[1];
	react[4] = new int[1];

	for (i = 0; i < dim_Ev_AB; i++)
		react[0][i] = 0;
	for (i = 0; i < dim_Ev_AC; i++)
		react[1][i] = 0;
	for (i = 0; i < dim_Ev_BC; i++)
		react[2][i] = 0;
	react[3][0] = 0;
	react[4][0] = 0;

	x = (double*)malloc(N * 6 * sizeof(double));

	for(i=0; i<N; i++)
		for(j=0; j<6; j++)
			x[i*6+j] = myrand();

#pragma omp parallel for num_threads(num_threads) private(j) shared(react) schedule(dynamic,1)
	for (j = 0; j < N; j++) {
		state = traject_A_BC(&x[j*6], Rmp, b, Tv, Vr, J, 0.04, 300., 1e-9, react);
	}

	printf("Number of trajectories:\n\t%d\n", N);

	printf("\n\nReaction A+BC(v\',J\')-->C+AB(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AB; i++)
		printf("%d\t%d\n", i, react[0][i]);

	printf("\n\nReaction A+BC(v\',J\')-->B+AC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_AC; i++) {
		printf("%d\t%d\n", i, react[1][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+BC(v\'\',J\'\')\nv\\J");
	printf("\n");
	for (i = 0; i < dim_Ev_BC; i++) {
		printf("%d\t%d\n", i, react[2][i]);
	}

	printf("\n\nReaction A+BC(v\',J\')-->A+B+C\n\t%d\n", react[3][0]);

	printf("\n\nReaction A+BC(v\',J\')-->ABC*\n\t%d\n", react[4][0]);

	free(x);

	for (i = 0; i < 5; i++)
		delete[] react[i];
	delete[] react;

	return 0;
}

int P_r_b_max_par(int N, double* Rmp, int J, double b_max, double Vr, double Tv, int** react, int num_threads) {

	int i, j;

	double *x;

	for (i = 0; i < dim_Ev_AB; i++)
		react[0][i] = 0;
	for (i = 0; i < dim_Ev_AC; i++)
		react[1][i] = 0;
	for (i = 0; i < dim_Ev_BC; i++)
		react[2][i] = 0;
	react[3][0] = 0;
	react[4][0] = 0;

	x = (double*)malloc(N * 6 * sizeof(double));

	for(i=0; i<N; i++)
		for(j=0; j<6; j++)
			x[i*6+j] = myrand();

#pragma omp parallel for num_threads(num_threads) private(j) shared(react) schedule(dynamic,1)
	for (j = 0; j < N; j++) {
		traject_A_BC(&x[j * 6], Rmp, b_max, Tv, Vr, J, 0.04, 400., 1e-9, react);
	}

	free(x);

	return 0;
}

int k_A_BC_QCT(int Ntr, int v1, int J1, double Ej, double* Etran, int dim_mesh, int num_threads) {

	int i, j, k, *** react, err;

	double dEa, Ev, mu, b_start, ki, kerr, Etr, T, rmp[2], Tv, *Vr, *bmax1;

	Vr = (double*)calloc(dim_mesh, sizeof(double));
	bmax1 = (double*)calloc(dim_mesh, sizeof(double));

	mu = mass[0] * (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]);		// Reduced mass A_BC

	for (i = 0; i < dim_mesh; i++)
		Vr[i] = sqrt(Etran[i] * 0.00964853326 * 2 / mu);			// from Etran in eV to Vr if A/fs

	react = new int** [dim_mesh];

	for (i = 0; i < dim_mesh; i++) {
		react[i] = new int* [5];

		react[i][0] = new int[dim_Ev_AB];
		for (j = 0; j < dim_Ev_AB; j++)
			react[i][0][j] = 0;

		react[i][1] = new int[dim_Ev_AC];
		for (j = 0; j < dim_Ev_AC; j++)
			react[i][1][j] = 0;

		react[i][2] = new int[dim_Ev_BC];
		for (j = 0; j < dim_Ev_BC; j++)
			react[i][2][j] = 0;

		react[i][3] = new int[1];
		react[i][3][0] = 0;

		react[i][4] = new int[1];
		react[i][4][0] = 0;
	}

	init_U();

	err = find_Rmp_BC(J1, Ej, rmp);

	printf("      Reaction A+BC(v\'=%d; J\'=%d)-->products\n\n", v1, J1);
	printf("      Coordinates (in A) of turning points at energy E(v\', J\')=%.10f\teV\n", Ej);
	printf("      R1=\t%.10f\tR2=\t%.10f\n", rmp[0], rmp[1]);

	Tv = find_Tv_BC_par(10000, rmp[1], J1, num_threads);

	printf("      Oscillation period at this energy level Tv=\t%.10f\tfs\n\n", Tv);

	printf("      Finding B_max parameter\n\nVr, A/fs\tEtr, eV\tbmax, A\n");

	b_start = 0.;
	for (i = 0; i < dim_mesh; i++) {
		bmax1[i] = find_b_max_A_BC_par(500, rmp, J1, b_start, Vr[i], Tv, num_threads);
		printf("%.10f\t%.10f\t%.4f\n", Vr[i], Etran[i], bmax1[i]);
		b_start = (double)((int)(9 * bmax1[i])) * 0.1;
	}

	for (i = 0; i < dim_mesh; i++) {
		if (bmax1[i] > 1e-4) {
			//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
			//            P_r_b_max(300, 6.547100242E-01, 9.079290792E-01, 0, bmax1[i], Vr[i], 8.003627381, m, order);

			//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
			printf("      Start of iteration %d/%d\n", i + 1, dim_mesh);
			P_r_b_max_par(Ntr, rmp, J1, bmax1[i], Vr[i], Tv, react[i], num_threads);
			printf("      End of iteration %d/%d\n\n", i + 1, dim_mesh);
		}
		else {
			printf("      Skip iteration %d/%d\n\n", i + 1, dim_mesh);
		}
	}

	Etr = mu * 0.5 / 0.00964853326 * 1.160363E+4; // in K/Vr/Vr

	// Reaction Reaction A+BC-->C+AB

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_AB; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_AB; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][0][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][0][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][0][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][0][j]) * sqrt((Ntr - react[i][0][j]) / (double)Ntr / (double)(react[i][0][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][0][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][0][j]) * sqrt((Ntr - react[i + 1][0][j]) / (double)Ntr / (double)(react[i + 1][0][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_AB; j++) {
				react1[0] += react[i][0][j];
				react1[1] += react[i + 1][0][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	// Reaction Reaction A+BC-->B+AC

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_AC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_AC; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][1][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][1][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][1][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][1][j]) * sqrt((Ntr - react[i][1][j]) / (double)Ntr / (double)(react[i][1][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][1][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][1][j]) * sqrt((Ntr - react[i + 1][1][j]) / (double)Ntr / (double)(react[i + 1][1][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_AC; j++) {
				react1[0] += react[i][1][j];
				react1[1] += react[i + 1][1][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	// Reaction Reaction A+BC-->A+BC

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_BC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_BC; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][2][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][2][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][2][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][2][j]) * sqrt((Ntr - react[i][2][j]) / (double)Ntr / (double)(react[i][2][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][2][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][2][j]) * sqrt((Ntr - react[i + 1][2][j]) / (double)Ntr / (double)(react[i + 1][2][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_BC; j++) {
				react1[0] += react[i][2][j];
				react1[1] += react[i + 1][2][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}


	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+B+C\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;
			react1[0] = react[i][3][0];
			react1[1] = react[i + 1][3][0];

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->ABC*\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;
			react1[0] = react[i][4][0];
			react1[1] = react[i + 1][4][0];

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	deinit_U();

	for (int i = 0; i < dim_mesh; i++) {
		for (int j = 0; j < 5; j++)
			delete[] react[i][j];
		delete[] react[i];
	}
	delete[] react;

	free(bmax1);
	free(Vr);

	return 0;

}

int k_A_BC_QCT(int Ntr, int v1, int J1, double Ej, double* Vr, double* bmax1, int dim_mesh, int num_threads) {

	int i, j, k, *** react, err;

	double dEa, Ev, mu, b_start, ki, kerr, Etr, T, rmp[2], Tv;

	mu = mass[0] * (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]);		// Reduced mass A_BC

	react = new int** [dim_mesh];

	for (i = 0; i < dim_mesh; i++) {
		react[i] = new int* [5];

		react[i][0] = new int[dim_Ev_AB];
		for (j = 0; j < dim_Ev_AB; j++)
			react[i][0][j] = 0;

		react[i][1] = new int[dim_Ev_AC];
		for (j = 0; j < dim_Ev_AC; j++)
			react[i][1][j] = 0;

		react[i][2] = new int[dim_Ev_BC];
		for (j = 0; j < dim_Ev_BC; j++)
			react[i][2][j] = 0;

		react[i][3] = new int[1];
		react[i][3][0] = 0;

		react[i][4] = new int[1];
		react[i][4][0] = 0;
	}

	init_U();

	err = find_Rmp_BC(J1, Ej, rmp);

	printf("      Reaction A+BC(v\'=%d; J\'=%d)-->products\n\n", v1, J1);
	printf("      Coordinates (in A) of turning points at energy E(v\', J\')=%.10f\teV\n", Ej);
	printf("      R1=\t%.10f\tR2=\t%.10f\n", rmp[0], rmp[1]);

	Tv = find_Tv_BC_par(1000, rmp[1], J1, num_threads);

	printf("      Oscillation period at this energy level Tv=\t%.10f\tfs\n\n", Tv);

	for (i = 0; i < dim_mesh; i++) {
		if (bmax1[i] > 1e-4) {
			//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
			//            P_r_b_max(300, 6.547100242E-01, 9.079290792E-01, 0, bmax1[i], Vr[i], 8.003627381, m, order);

			//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
			printf("      Start of iteration %d/%d\n", i + 1, dim_mesh);
			P_r_b_max_par(Ntr, rmp, J1, bmax1[i], Vr[i], Tv, react[i], num_threads);
			printf("      End of iteration %d/%d\n\n", i + 1, dim_mesh);
		}
		else {
			printf("      Skip iteration %d/%d\n\n", i + 1, dim_mesh);
		}
	}

	Etr = mu * 0.5 / 0.00964853326 * 1.160363E+4; // in K/Vr/Vr

	// Reaction Reaction A+BC-->C+AB

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_AB; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_AB; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][0][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][0][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][0][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][0][j]) * sqrt((Ntr - react[i][0][j]) / (double)Ntr / (double)(react[i][0][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][0][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][0][j]) * sqrt((Ntr - react[i + 1][0][j]) / (double)Ntr / (double)(react[i + 1][0][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_AB; j++) {
				react1[0] += react[i][0][j];
				react1[1] += react[i + 1][0][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	// Reaction Reaction A+BC-->B+AC

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_AC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_AC; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][1][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][1][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][1][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][1][j]) * sqrt((Ntr - react[i][1][j]) / (double)Ntr / (double)(react[i][1][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][1][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][1][j]) * sqrt((Ntr - react[i + 1][1][j]) / (double)Ntr / (double)(react[i + 1][1][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_AC; j++) {
				react1[0] += react[i][1][j];
				react1[1] += react[i + 1][1][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	// Reaction Reaction A+BC-->A+BC

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC(v\'\')\nT, K\t10^3/T\t", v1, J1);

	for (j = 0; j < dim_Ev_BC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("\n");
	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		for (j = 0; j < dim_Ev_BC; j++) {
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double ksum1, ksum2, ksum3, ksum4;
				ksum1 = bmax1[i] * bmax1[i] * (react[i][2][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][2][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react[i][2][j] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react[i][2][j]) * sqrt((Ntr - react[i][2][j]) / (double)Ntr / (double)(react[i][2][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react[i + 1][2][j] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[i + 1][2][j]) * sqrt((Ntr - react[i + 1][2][j]) / (double)Ntr / (double)(react[i + 1][2][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			printf("%.8e\t%.8e\t", ki, kerr);
		}
		printf("\n");
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;

			for (j = 0; j < dim_Ev_BC; j++) {
				react1[0] += react[i][2][j];
				react1[1] += react[i + 1][2][j];
			}

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}


	printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+B+C\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;
			react1[0] = react[i][3][0];
			react1[1] = react[i + 1][3][0];

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	printf("Reaction A+BC(v\'=%d; J\'=%d)-->ABC*\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

	for (T = 300; T < 4001; T += 100) {
		printf("%.2f\t%.7f\t", T, 1000. / T);
		ki = 0.;
		kerr = 0.;
		for (i = 0; i < dim_mesh - 1; i++) {
			double react1[2] = { 0., 0. };
			double ksum1, ksum2, ksum3, ksum4;
			react1[0] = react[i][4][0];
			react1[1] = react[i + 1][4][0];

			ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			if (react1[0] <= 0)
				ksum3 = 0.;
			else
				ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
			if (react1[1] <= 0)
				ksum4 = 0.;
			else
				ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
			ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
		}
		ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
		kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

		printf("%.8e\t%.8e\n", ki, kerr);
	}

	deinit_U();

	for (int i = 0; i < dim_mesh; i++) {
		for (int j = 0; j < 5; j++)
			delete[] react[i][j];
		delete[] react[i];
	}
	delete[] react;


	return 0;
}

int k_A_BC_QCT(int Ntr, int v1, double *Ej, int dimJ, double* Etran, int dim_mesh, int num_threads) { //int J1,

	int i, j, k, ***react, err, J1, dim_T=39;

	double dEa, Ev, mu, b_start, ki, kerr, Etr, T, rmp[2], Tv, *Vr, *bmax1, *wj, sum, *data1, ***data2;

	Vr = (double*)calloc(dim_mesh, sizeof(double));
	bmax1 = (double*)calloc(dim_mesh, sizeof(double));
	wj = (double*)calloc(dimJ, sizeof(double));
	data1 = (double*)calloc(dimJ * (dim_mesh + 5), sizeof(double));



	mu = mass[0] * (mass[1] + mass[2]) / (mass[0] + mass[1] + mass[2]);		// Reduced mass A_BC

	for (i = 0; i < dim_mesh; i++)
		Vr[i] = sqrt(Etran[i] * 0.00964853326 * 2 / mu);			// from Etran in eV to Vr if A/fs

	react = new int** [dim_mesh*dimJ];

	for (i = 0; i < dim_mesh * dimJ; i++) {
		react[i] = new int* [5];

		react[i][0] = new int[dim_Ev_AB];
		for (j = 0; j < dim_Ev_AB; j++)
			react[i][0][j] = 0;

		react[i][1] = new int[dim_Ev_AC];
		for (j = 0; j < dim_Ev_AC; j++)
			react[i][1][j] = 0;

		react[i][2] = new int[dim_Ev_BC];
		for (j = 0; j < dim_Ev_BC; j++)
			react[i][2][j] = 0;

		react[i][3] = new int[1];
		react[i][3][0] = 0;

		react[i][4] = new int[1];
		react[i][4][0] = 0;
	}

	data2 = new double** [dimJ];

	for (i = 0; i < dimJ; i++) {
		data2[i] = new double* [5];

		data2[i][0] = new double[(2*(dim_Ev_AB+1)+1)*dim_T];
		for (j = 0; j < (2*(dim_Ev_AB+1)+1)*dim_T; j++)
			data2[i][0][j] = 0;

		data2[i][1] = new double[(2*(dim_Ev_AC+1)+1)*dim_T];
		for (j = 0; j < (2*(dim_Ev_AC+1)+1)*dim_T; j++)
			data2[i][1][j] = 0;

		data2[i][2] = new double[(2*(dim_Ev_BC+1)+1)*dim_T];
		for (j = 0; j < (2*(dim_Ev_BC+1)+1)*dim_T; j++)
			data2[i][2][j] = 0;
		
		data2[i][3] = new double[3*dim_T];
		for (j = 0; j < 3*dim_T; j++)
			data2[i][3][j] = 0;

		data2[i][4] = new double[3*dim_T];
		for (j = 0; j < 3*dim_T; j++)
			data2[i][4][j] = 0;
	}

	init_U();

	for (J1=0; J1<dimJ; J1++) {

		err = find_Rmp_BC(J1, Ej[J1], rmp);

		data1[J1*(dim_mesh+5)] = J1;
		data1[J1*(dim_mesh+5)+1] = Ej[J1];
		data1[J1*(dim_mesh+5)+2] = rmp[0];
		data1[J1*(dim_mesh+5)+3] = rmp[1];

		printf("      Reaction A+BC(v\'=%d; J\'=%d)-->products\n\n", v1, J1);
		printf("      Coordinates (in A) of turning points at energy E(v\', J\')=%.10f\teV\n", Ej[J1]);
		printf("      R1=\t%.10f\tR2=\t%.10f\n", rmp[0], rmp[1]);

		Tv = find_Tv_BC_par(10000, rmp[1], J1, num_threads);
		
		data1[J1*(dim_mesh+5)+4] = Tv;

		printf("      Oscillation period at this energy level Tv=\t%.10f\tfs\n\n", Tv);

		printf("      Finding B_max parameter\n\nVr, A/fs\tEtr, eV\tbmax, A\n");

		b_start = 0.;
		for (i = 0; i < dim_mesh; i++) {
			bmax1[i] = find_b_max_A_BC_par(700, rmp, J1, b_start, Vr[i], Tv, num_threads);
			printf("%.10f\t%.10f\t%.4f\n", Vr[i], Etran[i], bmax1[i]);
			if(J1>0)
				b_start = data1[(J1-1)*(dim_mesh+5)+5+i];
			else
				b_start = (double)((int)(9 * bmax1[i])) * 0.1;
		}

		for(i=0; i<dim_mesh; i++)
			data1[J1*(dim_mesh+5)+5+i] = bmax1[i];

		for (i = 0; i < dim_mesh; i++) {
			if (bmax1[i] > 1e-4) {
				//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
				//            P_r_b_max(300, 6.547100242E-01, 9.079290792E-01, 0, bmax1[i], Vr[i], 8.003627381, m, order);

				//            printf("\nBmax=%.5f\tVr=%.7f\n", bmax1[i], Vr[i]);
				printf("      Start of iteration %d/%d\n", i + 1, dim_mesh);
				P_r_b_max_par(Ntr, rmp, J1, bmax1[i], Vr[i], Tv, react[J1*dim_mesh+i], num_threads);
				printf("      End of iteration %d/%d\n\n", i + 1, dim_mesh);
			}
			else {
				printf("      Skip iteration %d/%d\n\n", i + 1, dim_mesh);
			}
		}

		Etr = mu * 0.5 / 0.00964853326 * 1.160363E+4; // in K/Vr/Vr

		// Reaction Reaction A+BC-->C+AB

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB(v\'\')\nT, K\t10^3/T\t", v1, J1);

		for (j = 0; j < dim_Ev_AB; j++)
			printf("k%d\tk_err%d\t", j, j);
		printf("\n");
		k = 0;
		for (T = 200; T < 4001; T += 100) {

			data2[J1][0][k*(2*(dim_Ev_AB+1)+1)] = T;

			printf("%.2f\t%.7f\t", T, 1000. / T);
			for (j = 0; j < dim_Ev_AB; j++) {
				ki = 0.;
				kerr = 0.;
				for (i = 0; i < dim_mesh - 1; i++) {
					double ksum1, ksum2, ksum3, ksum4;
					ksum1 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][0][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][0][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					if (react[J1*dim_mesh+i][0][j] <= 0)
						ksum3 = 0.;
					else
						ksum3 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][0][j]) * sqrt((Ntr - react[J1*dim_mesh+i][0][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i][0][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					if (react[J1*dim_mesh+i+1][0][j] <= 0)
						ksum4 = 0.;
					else
						ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][0][j]) * sqrt((Ntr - react[J1*dim_mesh+i+1][0][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i+1][0][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
					kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				}
				ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
				kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

				data2[J1][0][k*(2*(dim_Ev_AB+1)+1)+1+2*j] = ki;
				data2[J1][0][k*(2*(dim_Ev_AB+1)+1)+2+2*j] = kerr;

				printf("%.8e\t%.8e\t", ki, kerr);
			}
			printf("\n");
			k++;
		}
		printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB\nEtr, eV\ts_all\tS_err_all\n", v1, J1);

		for(i=0; i<dim_mesh; i++){
			double react1 = 0.;
			double sigma, dsigma;
			for (j=0; j<dim_Ev_AB; j++)
				react1 += react[J1*dim_mesh+i][0][j];
			sigma = M_PI * sqr(bmax1[i]) * (double)react1 / Ntr;
			if (react1 == 0)
				dsigma = 0.;
			else
				dsigma = sigma * sqrt((Ntr - react1) / (double)Ntr / (double)(react1));
			printf("%.5f\t%.10e\t%.10e\n", Etran[i], sigma, dsigma);
		}

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->C+AB\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

		k = 0;
		for (T = 200; T < 4001; T += 100) {
			printf("%.2f\t%.7f\t", T, 1000. / T);
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double react1[2] = { 0., 0. };
				double ksum1, ksum2, ksum3, ksum4;

				for (j = 0; j < dim_Ev_AB; j++) {
					react1[0] += react[J1*dim_mesh+i][0][j];
					react1[1] += react[J1*dim_mesh+i+1][0][j];
				}

				ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react1[0] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react1[1] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			data2[J1][0][k*(2*(dim_Ev_AB+1)+1)+1+2*dim_Ev_AB] = ki;
			data2[J1][0][k*(2*(dim_Ev_AB+1)+1)+2+2*dim_Ev_AB] = kerr;

			printf("%.8e\t%.8e\n", ki, kerr);
			k++;
		}

		// Reaction Reaction A+BC-->B+AC

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC(v\'\')\nT, K\t10^3/T\t", v1, J1);

		for (j = 0; j < dim_Ev_AC; j++)
			printf("k%d\tk_err%d\t", j, j);
		printf("\n");
		k = 0;
		for (T = 200; T < 4001; T += 100) {

			data2[J1][1][k*(2*(dim_Ev_AC+1)+1)] = T;

			printf("%.2f\t%.7f\t", T, 1000. / T);
			for (j = 0; j < dim_Ev_AC; j++) {
				ki = 0.;
				kerr = 0.;
				for (i = 0; i < dim_mesh - 1; i++) {
					double ksum1, ksum2, ksum3, ksum4;
					ksum1 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][1][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][1][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					if (react[J1*dim_mesh+i][1][j] <= 0)
						ksum3 = 0.;
					else
						ksum3 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][1][j]) * sqrt((Ntr - react[J1*dim_mesh+i][1][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i][1][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					if (react[J1*dim_mesh+i+1][1][j] <= 0)
						ksum4 = 0.;
					else
						ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][1][j]) * sqrt((Ntr - react[J1*dim_mesh+i+1][1][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i+1][1][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
					kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				}
				ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
				kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

				data2[J1][1][k*(2*(dim_Ev_AC+1)+1)+1+2*j] = ki;
				data2[J1][1][k*(2*(dim_Ev_AC+1)+1)+2+2*j] = kerr;

				printf("%.8e\t%.8e\t", ki, kerr);
			}
			printf("\n");
			k++;
		}

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC\nEtr, eV\ts_all\tS_err_all\n", v1, J1);

		for (i = 0; i < dim_mesh; i++) {
			double react1 = 0.;
			double sigma, dsigma;
			for (j = 0; j < dim_Ev_AC; j++)
				react1 += react[J1 * dim_mesh + i][0][j];
			sigma = M_PI * sqr(bmax1[i]) * (double)react1 / Ntr;
			if (react1 == 0)
				dsigma = 0.;
			else
				dsigma = sigma * sqrt((Ntr - react1) / (double)Ntr / (double)(react1));
			printf("%.5f\t%.10e\t%.10e\n", Etran[i], sigma, dsigma);
		}

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->B+AC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

		k = 0;
		for (T = 200; T < 4001; T += 100) {
			printf("%.2f\t%.7f\t", T, 1000. / T);
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double react1[2] = { 0., 0. };
				double ksum1, ksum2, ksum3, ksum4;

				for (j = 0; j < dim_Ev_AC; j++) {
					react1[0] += react[J1*dim_mesh+i][1][j];
					react1[1] += react[J1*dim_mesh+i+1][1][j];
				}

				ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react1[0] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react1[1] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			data2[J1][1][k*(2*(dim_Ev_AC+1)+1)+1+2*dim_Ev_AC] = ki;
			data2[J1][1][k*(2*(dim_Ev_AC+1)+1)+2+2*dim_Ev_AC] = kerr;

			printf("%.8e\t%.8e\n", ki, kerr);
			k++;
		}

		// Reaction Reaction A+BC-->A+BC

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC(v\'\')\nT, K\t10^3/T\t", v1, J1);

		for (j = 0; j < dim_Ev_BC; j++)
			printf("k%d\tk_err%d\t", j, j);
		printf("\n");
		k = 0;
		for (T = 200; T < 4001; T += 100) {

			data2[J1][2][k*(2*(dim_Ev_BC+1)+1)] = T;

			printf("%.2f\t%.7f\t", T, 1000. / T);
			for (j = 0; j < dim_Ev_BC; j++) {
				ki = 0.;
				kerr = 0.;
				for (i = 0; i < dim_mesh - 1; i++) {
					double ksum1, ksum2, ksum3, ksum4;
					ksum1 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][2][j]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][2][j]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					if (react[J1*dim_mesh+i][2][j] <= 0)
						ksum3 = 0.;
					else
						ksum3 = bmax1[i] * bmax1[i] * (react[J1*dim_mesh+i][2][j]) * sqrt((Ntr - react[J1*dim_mesh+i][2][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i][2][j])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
					if (react[J1*dim_mesh+i+1][2][j] <= 0)
						ksum4 = 0.;
					else
						ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react[J1*dim_mesh+i+1][2][j]) * sqrt((Ntr - react[J1*dim_mesh+i+1][2][j]) / (double)Ntr / (double)(react[J1*dim_mesh+i+1][2][j])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
					ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
					kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				}
				ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
				kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

				data2[J1][2][k*(2*(dim_Ev_BC+1)+1)+1+2*j] = ki;
				data2[J1][2][k*(2*(dim_Ev_BC+1)+1)+2+2*j] = kerr;

				printf("%.8e\t%.8e\t", ki, kerr);
			}
			printf("\n");
			k++;
		}

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+BC\nT, K\t10^3/T\tk_all\tk_err_all\n", v1, J1);

		k = 0;
		for (T = 200; T < 4001; T += 100) {
			printf("%.2f\t%.7f\t", T, 1000. / T);
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double react1[2] = { 0., 0. };
				double ksum1, ksum2, ksum3, ksum4;

				for (j = 0; j < dim_Ev_BC; j++) {
					react1[0] += react[J1*dim_mesh+i][2][j];
					react1[1] += react[J1*dim_mesh+i+1][2][j];
				}

				ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react1[0] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react1[1] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			data2[J1][2][k*(2*(dim_Ev_BC+1)+1)+1+2*dim_Ev_BC] = ki;
			data2[J1][2][k*(2*(dim_Ev_BC+1)+1)+2+2*dim_Ev_BC] = kerr;

			printf("%.8e\t%.8e\n", ki, kerr);
			k++;
		}


		printf("Reaction A+BC(v\'=%d; J\'=%d)-->A+B+C\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

		k = 0;
		for (T = 200; T < 4001; T += 100) {

			data2[J1][3][k*3] = T;

			printf("%.2f\t%.7f\t", T, 1000. / T);
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double react1[2] = { 0., 0. };
				double ksum1, ksum2, ksum3, ksum4;
				react1[0] = react[J1*dim_mesh+i][3][0];
				react1[1] = react[J1*dim_mesh+i+1][3][0];

				ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react1[0] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react1[1] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			data2[J1][3][k*3+1] = ki;
			data2[J1][3][k*3+2] = kerr;

			printf("%.8e\t%.8e\n", ki, kerr);
			k++;
		}

		printf("Reaction A+BC(v\'=%d; J\'=%d)-->ABC*\nT, K\t10^3/T\tk\tk_err\n", v1, J1);

		k = 0;
		for (T = 200; T < 4001; T += 100) {

			data2[J1][4][k*3] = T;

			printf("%.2f\t%.7f\t", T, 1000. / T);
			ki = 0.;
			kerr = 0.;
			for (i = 0; i < dim_mesh - 1; i++) {
				double react1[2] = { 0., 0. };
				double ksum1, ksum2, ksum3, ksum4;
				react1[0] = react[J1*dim_mesh+i][4][0];
				react1[1] = react[J1*dim_mesh+i+1][4][0];

				ksum1 = bmax1[i] * bmax1[i] * (react1[0]) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				ksum2 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				if (react1[0] <= 0)
					ksum3 = 0.;
				else
					ksum3 = bmax1[i] * bmax1[i] * (react1[0]) * sqrt((Ntr - react1[0]) / (double)Ntr / (double)(react1[0])) * exp(-Etr * Vr[i] * Vr[i] / T) * Vr[i] * Vr[i];
				if (react1[1] <= 0)
					ksum4 = 0.;
				else
					ksum4 = bmax1[i + 1] * bmax1[i + 1] * (react1[1]) * sqrt((Ntr - react1[1]) / (double)Ntr / (double)(react1[1])) * exp(-Etr * Vr[i + 1] * Vr[i + 1] / T) * Vr[i + 1] * Vr[i + 1];
				ki += (ksum1 + ksum2) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
				kerr += (ksum3 + ksum4) * (Vr[i + 1] * Vr[i + 1] - Vr[i] * Vr[i]);
			}
			ki *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;
			kerr *= 0.5 * Etr * Etr / T * 8.762702193E+11 / sqrt(mu * T) / Ntr * M_PI;

			data2[J1][4][k*3+1] = ki;
			data2[J1][4][k*3+2] = kerr;

			printf("%.8e\t%.8e\n", ki, kerr);
			k++;
		}
	}

	printf("\n\n    Summary of Reaction A+BC(v\'=%d; J\'=%d)-->products\n\n", v1, J1);

	printf("\t\t\t\t\t");
	for(i=0; i<dim_mesh/2; i++)
		printf("\t");
	printf("E_transl, eV\n");
	printf("J\tE(J), eV\tRm, A\tRp, A\tTv, fs");
	for(i=0; i<dim_mesh; i++)
		printf("\tbmax%d, A", i);
	printf("\n");

	for(i=0; i<dimJ; i++){
		printf("%.0f", data1[i*(dim_mesh+5)]);
		for(j=1; j<5; j++)
			printf("\t%.10f", data1[i*(dim_mesh+5)+j]);
		for(j=5; j<dim_mesh+5; j++)
			printf("\t%.3f", data1[i*(dim_mesh+5)+j]);
		printf("\n");
	}
	printf("\n\n     Averaged over J number Rate constants for Reaction A+BC(v\')\n\n");


	printf("     Reaction A+BC(v\'=%d)-->C+AB(v\'\')\nT, K\t10^3/T\t", v1);
	for (j = 0; j < dim_Ev_AB; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("k_all\tk_err_all");
	printf("\n");
	for(i=0; i<dim_T; i++){
		
		T = data2[0][0][i*(2*(dim_Ev_AB+1)+1)];

		printf("%.2f\t%.10f", T, 1000 / T);

		sum = 0.;
		for(k=0; k<dimJ; k++)
			sum += (2*k+1)*exp(-Ej[k]*1.160363E+4/T);
		for(k=0; k<dimJ; k++)
			wj[k] = (2*k+1)*exp(-Ej[k]*1.160363E+4/T)/sum;
		
		for(k=0; k<dim_Ev_AB+1; k++){
			ki = 0;
			kerr = 0;

			for(j=0; j<dimJ; j++){
				ki += data2[j][0][i*(2*(dim_Ev_AB+1)+1)+1+2*k] * wj[j];
				kerr += data2[j][0][i*(2*(dim_Ev_AB+1)+1)+2+2*k] * wj[j];
			}

			printf("\t%.10e\t%.10e", ki, kerr);
		}
		printf("\n");
	}
	printf("\n\n");

	printf("     Reaction A+BC(v\'=%d)-->B+AC(v\'\')\nT, K\t10^3/T\t", v1);
	for (j = 0; j < dim_Ev_AC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("k_all\tk_err_all");
	printf("\n");
	for(i=0; i<dim_T; i++){
		
		T = data2[0][1][i*(2*(dim_Ev_AC+1)+1)];

		printf("%.2f\t%.10f", T, 1000 / T);

		sum = 0.;
		for(k=0; k<dimJ; k++)
			sum += (2*k+1)*exp(-Ej[k]*1.160363E+4/T);
		for(k=0; k<dimJ; k++)
			wj[k] = (2*k+1)*exp(-Ej[k]*1.160363E+4/T)/sum;
		
		for(k=0; k<dim_Ev_AC+1; k++){
			ki = 0;
			kerr = 0;

			for(j=0; j<dimJ; j++){
				ki += data2[j][1][i*(2*(dim_Ev_AC+1)+1)+1+2*k] * wj[j];
				kerr += data2[j][1][i*(2*(dim_Ev_AC+1)+1)+2+2*k] * wj[j];
			}

			printf("\t%.10e\t%.10e", ki, kerr);
		}
		printf("\n");
	}
	printf("\n\n");

	printf("     Reaction A+BC(v\'=%d)-->A+BC(v\'\')\nT, K\t10^3/T\t", v1);
	for (j = 0; j < dim_Ev_BC; j++)
		printf("k%d\tk_err%d\t", j, j);
	printf("k_all\tk_err_all");
	printf("\n");
	for(i=0; i<dim_T; i++){
		
		T = data2[0][2][i*(2*(dim_Ev_BC+1)+1)];

		printf("%.2f\t%.10f", T, 1000 / T);

		sum = 0.;
		for(k=0; k<dimJ; k++)
			sum += (2*k+1)*exp(-Ej[k]*1.160363E+4/T);
		for(k=0; k<dimJ; k++)
			wj[k] = (2*k+1)*exp(-Ej[k]*1.160363E+4/T)/sum;
		
		for(k=0; k<dim_Ev_BC+1; k++){
			ki = 0;
			kerr = 0;

			for(j=0; j<dimJ; j++){
				ki += data2[j][2][i*(2*(dim_Ev_BC+1)+1)+1+2*k] * wj[j];
				kerr += data2[j][2][i*(2*(dim_Ev_BC+1)+1)+2+2*k] * wj[j];
			}

			printf("\t%.10e\t%.10e", ki, kerr);
		}
		printf("\n");
	}
	printf("\n\n");


	printf("Reaction A+BC(v\'=%d)-->A+B+C\nT, K\t10^3/T\tk\tk_err\n", v1);
	for(i=0; i<dim_T; i++){
		
		T = data2[0][3][i*3];

		printf("%.2f\t%.10f", T, 1000 / T);

		sum = 0.;
		for(k=0; k<dimJ; k++)
			sum += (2*k+1)*exp(-Ej[k]*1.160363E+4/T);
		for(k=0; k<dimJ; k++)
			wj[k] = (2*k+1)*exp(-Ej[k]*1.160363E+4/T)/sum;
		
		ki = 0;
		kerr = 0;

		for(j=0; j<dimJ; j++){
			ki += data2[j][3][i*3+1] * wj[j];
			kerr += data2[j][3][i*3+2] * wj[j];
		}

		printf("\t%.10e\t%.10e\n", ki, kerr);
		
	}
	printf("\n\n");

	printf("Reaction A+BC(v\'=%d)-->ABC*\nT, K\t10^3/T\tk\tk_err\n", v1);
	for(i=0; i<dim_T; i++){
		
		T = data2[0][4][i*3];

		printf("%.2f\t%.10f", T, 1000 / T);

		sum = 0.;
		for(k=0; k<dimJ; k++)
			sum += (2*k+1)*exp(-Ej[k]*1.160363E+4/T);
		for(k=0; k<dimJ; k++)
			wj[k] = (2*k+1)*exp(-Ej[k]*1.160363E+4/T)/sum;
		
		ki = 0;
		kerr = 0;

		for(j=0; j<dimJ; j++){
			ki += data2[j][4][i*3+1] * wj[j];
			kerr += data2[j][4][i*3+2] * wj[j];
		}

		printf("\t%.10e\t%.10e\n", ki, kerr);
		
	}
	printf("\n\n");

	deinit_U();

	for (int i = 0; i < dimJ; i++) {
		for (int j = 0; j < 5; j++)
			delete[] data2[i][j];
		delete[] data2[i];
	}
	delete[] data2;

	for (int i = 0; i < dim_mesh*dimJ; i++) {
		for (int j = 0; j < 5; j++)
			delete[] react[i][j];
		delete[] react[i];
	}
	delete[] react;

	free(data1);
	free(wj);
	free(bmax1);
	free(Vr);

	return 0;

}
