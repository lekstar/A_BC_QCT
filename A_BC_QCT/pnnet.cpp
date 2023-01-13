#include "pnnet.h"

void NNet::init_laylink(const int n, const int *a){
	layers = n;
	
	link = new int[n];
	for(int i=0; i<n; i++)
		link[i] = a[i];
	
	laylink_state = 1;
}
		
int NNet::init_weight(const double *w){
	if(laylink_state==0)
		return 1;
	
	int i, j, k, l;
	
	weight = new double*[layers];

	l = 0;
	
	for(i=0; i<layers-1; i++){
		weight[i] = new double[link[i+1]*link[i]];
		for(j=0; j<link[i+1]; j++){
			for(k=0; k<link[i]; k++){
				weight[i][j*link[i]+k] = w[l];
				l++;
			}
		}
	}
	weight_state = 1;
	
	return 0;
}

int NNet::init_bias(const double *b){
	if(laylink_state==0)
		return 1;
	
	int i, j, l = 0;
	
	bias = new double*[layers-1];

	for(i=0; i<layers-1; i++){
		bias[i] = new double[link[i+1]];
		for(j=0; j<link[i+1]; j++){
			bias[i][j] = b[l];
			l++;
		}	
	}
	
	bias_state = 1;
	
	return 0;
}
		
int NNet::init_Sij(){
	if(laylink_state==0)
		return 1;
	
//	int i, j, k, l;
	
//	Sij = new double*[layers];
//
//	for(i=0; i<layers; i++){
//		Sij[i] = new double[link[i]];
//		for(j=0; j<link[i]; j++)
//			Sij[i][j] = 0.;	
//	}
//	
//	Gij = new double**[link[0]];
//
//	for(i=0; i<link[0]; i++){
//		Gij[i] = new double*[layers];
//		for(j=0; j<layers; j++){
//			Gij[i][j] = new double[link[j]];
//		}
//	}
//	
//	Hij = new double**[link[0]*link[0]];
//
//	for(i=0; i<link[0]*link[0]; i++){
//		Hij[i] = new double*[layers];
//		for(j=0; j<layers; j++){
//			Hij[i][j] = new double[link[j]];
//		}
//	}
//	
//	Sij_state = 1;
	
	return 0;
}
		
void NNet::init_net(const int n, const int *a, const double *w, const double *b){
	if(net_state==0){
		init_laylink(n, a);
		init_weight(w);
		init_bias(b);
//		init_Sij();
	}
	
	net_state = 1;
}

int NNet::is_init(){
	
	return net_state;
}

void NNet::delete_net(){
	if(net_state!=0){
		
		for(int i=0; i<layers-1; i++){
			delete [] weight[i];
		}
		delete [] weight;
		
		for(int i=0; i<layers-1; i++)
			delete [] bias[i];
		delete [] bias;
		
//		for(int i=0; i<layers; i++)
//			delete [] Sij[i];
//		delete [] Sij;
//		
//		for(int i=0; i<link[0]; i++){
//			for(int j=0; j<layers; j++)
//				delete [] Gij[i][j];
//			delete [] Gij[i];
//		}
//		delete [] Gij;
//		
//		for(int i=0; i<link[0]*link[0]; i++){
//			for(int j=0; j<layers; j++)
//				delete [] Hij[i][j];
//			delete [] Hij[i];
//		}
//		delete [] Hij;
						
		delete [] link;
		
		laylink_state = 0;
		weight_state = 0;
		bias_state = 0;
//		Sij_state = 0;
		net_state = 0;
	}
	
	if(weight_state!=0){
		for(int i=0; i<layers-1; i++){
			delete [] weight[i];
		}
		delete [] weight;
		
		weight_state = 0;
	}
	
	if(bias_state!=0){
		for(int i=0; i<layers-1; i++)
			delete [] bias[i];
		delete [] bias;
		
		bias_state = 0;
	}
	
//	if(Sij_state!=0){
//		for(int i=0; i<layers; i++)
//			delete [] Sij[i];
//		delete [] Sij;
//		
//		for(int i=0; i<link[0]; i++){
//			for(int j=0; j<layers; j++)
//				delete [] Gij[i][j];
//			delete [] Gij[i];
//		}
//		delete [] Gij;
//		
//		for(int i=0; i<link[0]*link[0]; i++){
//			for(int j=0; j<layers; j++)
//				delete [] Hij[i][j];
//			delete [] Hij[i];
//		}
//		delete [] Hij;
//		
//		Sij_state = 0;
//	}
	
	if(laylink_state!=0){
		delete [] link;
				
		laylink_state = 0;
	}
}
		
int NNet::out(double *input, double *output){
	if(net_state==0)
		return 1;
	
	int i, j, k;
	double **Sij, sum;
	
	Sij = new double*[layers];

	for(i=0; i<layers; i++){
		//printf("link=%d\n", link[i]);
		Sij[i] = new double[link[i]];
		for(j=0; j<link[i]; j++)
			Sij[i][j] = 0.;	
	}
	
	for(i=0; i<link[0]; i++)
		Sij[0][i] = input[i];
	
	for(i=0; i<layers-1; i++){
		
		for(j=0; j<link[i+1]; j++)
			Sij[i+1][j] = bias[i][j];

		cblas_dgemv(CblasRowMajor, CblasNoTrans, link[i+1], link[i], 1.0, weight[i], link[i], Sij[i], 1, 1.0, Sij[i+1], 1);
		
		/*for(j=0; j<link[i+1]; j++){
			sum = 0.;
			#pragma omp simd reduction(+:sum)
				for(k=0; k<link[i]; k++)
					sum += Sij[i][k]*weight[i][j*link[i]+k];
				Sij[i+1][j] = sum + bias[i][j];
			}*/
		
		for(j=0; j<link[i+1]; j++)
			if(i!=layers-2)
				Sij[i+1][j] = tanh(Sij[i+1][j]);
	}

	for(i=0; i<link[layers-1]; i++)
		output[i] = Sij[layers-1][i];
	
	for(int i=0; i<layers; i++)
			delete [] Sij[i];
		delete [] Sij;
	
	return 0;
}
		
int NNet::grad_input(double *input, double *output, double *grad){
	if(net_state==0)
		return 1;
	
	int i, j, k, m;
	double **Sij, ***Gij, sum;
	
	Sij = new double*[layers];

    if(Sij==NULL){
        printf("\nCan not allocate memory for Sij\n");
        return 1;
    }

	for(i=0; i<layers; i++){
		Sij[i] = new double[link[i]];
        if(Sij[i]==NULL){
            printf("\nCan not allocate memory for Sij\n");
            return 1;
        }
		for(j=0; j<link[i]; j++)
			Sij[i][j] = 0.;	
	}
	
	Gij = new double**[link[0]];

	for(i=0; i<link[0]; i++){
		Gij[i] = new double*[layers];
		for(j=0; j<layers; j++){
			Gij[i][j] = new double[link[j]];
		}
	}
	
	for(i=0; i<link[0]; i++)
		Sij[0][i] = input[i];
	
	for(i=0; i<layers-1; i++){

		for(j=0; j<link[i+1]; j++)
			Sij[i+1][j] = bias[i][j];

		cblas_dgemv(CblasRowMajor, CblasNoTrans, link[i+1], link[i], 1.0, weight[i], link[i], Sij[i], 1, 1.0, Sij[i+1], 1);

		/*for(j=0; j<link[i+1]; j++){
				sum = 0.;
				#pragma omp simd reduction(+:sum)
				for(k=0; k<link[i]; k++)
					sum += Sij[i][k]*weight[i][j*link[i]+k];
				Sij[i+1][j] = sum + bias[i][j];
			}*/
		
		for(j=0; j<link[i+1]; j++)
			if(i!=layers-2)
				Sij[i+1][j] = tanh(Sij[i+1][j]);
	}

	for(i=0; i<link[layers-1]; i++)
		output[i] = Sij[layers-1][i];
	
	for(m=0; m<link[0]; m++){
		for(i=0; i<link[0]; i++)
			Gij[m][0][i] = 0.;
		Gij[m][0][m] = 1.;				
		
		for(i=0; i<layers-1; i++){

			for(j=0; j<link[i+1]; j++)
				Gij[m][i+1][j] = 0.;

			cblas_dgemv(CblasRowMajor, CblasNoTrans, link[i+1], link[i], 1.0, weight[i], link[i], Gij[m][i], 1, 1.0, Gij[m][i+1], 1);

			//for(j=0; j<link[i+1]; j++){
			//	sum = 0.;
			//	#pragma omp simd reduction(+:sum)
			//		for(k=0; k<link[i]; k++)
			//				sum+= Gij[m][i][k]*weight[i][j*link[i]+k];
			//	Gij[m][i+1][j] = sum;
			//}
				
			for(j=0; j<link[i+1]; j++)
				if(i!=layers-2)
					Gij[m][i+1][j] *= (1.-Sij[i+1][j]*Sij[i+1][j]);
		}
		
		grad[m] = Gij[m][layers-1][0];
	}
	
	for(int i=0; i<layers; i++)
			delete [] Sij[i];
		delete [] Sij;
	
	for(int i=0; i<link[0]; i++){
		for(int j=0; j<layers; j++)
			delete [] Gij[i][j];
		delete [] Gij[i];
	}
	delete [] Gij;
	
	return 0;
}
		
int NNet::hess_grad_input(double *input, double *output, double *grad, double *hess){
	if(net_state==0)
		return 1;
	
	int i, j, k, m, n;
	double sum1, sum2, sum3;
	double **Sij, ***Gij, ***Hij;
	
	Sij = new double*[layers];

	for(i=0; i<layers; i++){
		Sij[i] = new double[link[i]];
		for(j=0; j<link[i]; j++)
			Sij[i][j] = 0.;	
	}
	
	Gij = new double**[link[0]];

	for(i=0; i<link[0]; i++){
		Gij[i] = new double*[layers];
		for(j=0; j<layers; j++){
			Gij[i][j] = new double[link[j]];
		}
	}
	
	Hij = new double**[link[0]*link[0]];

	for(i=0; i<link[0]*link[0]; i++){
		Hij[i] = new double*[layers];
		for(j=0; j<layers; j++){
			Hij[i][j] = new double[link[j]];
		}
	}
	
	for(i=0; i<link[0]; i++)
		Sij[0][i] = input[i];
	
	for(i=0; i<layers-1; i++){
		for(j=0; j<link[i+1]; j++)
			Sij[i+1][j] = bias[i][j];

		cblas_dgemv(CblasRowMajor, CblasNoTrans, link[i+1], link[i], 1.0, weight[i], link[i], Sij[i], 1, 1.0, Sij[i+1], 1);
		
		/*for(j=0; j<link[i+1]; j++){
			sum = 0.;
			#pragma omp simd reduction(+:sum)
				for(k=0; k<link[i]; k++)
					sum += Sij[i][k]*weight[i][j*link[i]+k];
				Sij[i+1][j] = sum + bias[i][j];
			}*/
		
		for(j=0; j<link[i+1]; j++)
			if(i!=layers-2)
				Sij[i+1][j] = tanh(Sij[i+1][j]);
	}
		

	for(i=0; i<link[layers-1]; i++)
		output[i] = Sij[layers-1][i];
	
	for(m=0; m<link[0]; m++){
		for(i=0; i<link[0]; i++)
			if(i==m)
				Gij[m][0][i] = 1.;
			else
				Gij[m][0][i] = 0.;
		
		for(i=0; i<layers-1; i++){
			for(j=0; j<link[i+1]; j++)
				Gij[m][i+1][j] = 0.;

			cblas_dgemv(CblasRowMajor, CblasNoTrans, link[i+1], link[i], 1.0, weight[i], link[i], Gij[m][i], 1, 1.0, Gij[m][i+1], 1);

			//for(j=0; j<link[i+1]; j++){
			//	sum = 0.;
			//	#pragma omp simd reduction(+:sum)
			//		for(k=0; k<link[i]; k++)
			//				sum+= Gij[m][i][k]*weight[i][j*link[i]+k];
			//	Gij[m][i+1][j] = sum;
			//}
				
			for(j=0; j<link[i+1]; j++)
				if(i!=layers-2)
					Gij[m][i+1][j] *= (1.-Sij[i+1][j]*Sij[i+1][j]);
		}
		
		grad[m] = Gij[m][layers-1][0];
	}
	
	for(n=0; n<link[0]; n++)
		for(m=0; m<link[0]; m++){
			for(i=0; i<link[0]; i++)
				Hij[n*link[0]+m][0][i] = 0.;
			
			for(i=0; i<layers-1; i++){
				for(j=0; j<link[i+1]; j++){
					Hij[n*link[0]+m][i+1][j] = 0.;
					sum1 = 0.;
					sum2 = 0.;
					sum3 = 0.;
					for(k=0; k<link[i]; k++){
						sum1 += Hij[n*link[0]+m][i][k]*weight[i][j*link[i]+k];
						sum2 += Gij[m][i][k]*weight[i][j*link[i]+k];
						sum3 += Gij[n][i][k]*weight[i][j*link[i]+k];
					}
					if(i!=layers-2)
						Hij[n*link[0]+m][i+1][j] = (sum1 + sum2*sum3*(-2.*Sij[i+1][j]) ) * (1.-Sij[i+1][j]*Sij[i+1][j]);
					else
						Hij[n*link[0]+m][i+1][j] = sum1;
				}
			}
			
			hess[n*link[0]+m] = Hij[n*link[0]+m][layers-1][0];
		}
	
	for(int i=0; i<link[0]*link[0]; i++){
		for(int j=0; j<layers; j++)
			delete [] Hij[i][j];
		delete [] Hij[i];
	}
	delete [] Hij;
	
	for(int i=0; i<layers; i++)
			delete [] Sij[i];
		delete [] Sij;
	
	for(int i=0; i<link[0]; i++){
		for(int j=0; j<layers; j++)
			delete [] Gij[i][j];
		delete [] Gij[i];
	}
	delete [] Gij;
	
	return 0;
}
