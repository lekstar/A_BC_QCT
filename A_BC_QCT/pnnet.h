#ifndef PNNET_H
#define PNNET_H

#include <cstdlib>
#include <iostream>
#include <cmath>

#include <mkl.h>

#include <omp.h>

#ifndef sqr(x)
#define sqr(x)	((x)*(x))
#endif

class NNet{
	// Private section
	public:
		// Public Declarations
		
		NNet(){
			laylink_state = 0;
			weight_state = 0;
			bias_state = 0;
			net_state = 0;
		}
		
		~NNet(){
			delete_net();
		}
		
		void init_laylink(const int n, const int *a);
		
		int init_weight(const double *w);
		
		int init_bias(const double *b);
		
		int init_Sij();
		
		void init_net(const int n, const int *a, const double *w, const double *b);
		
		void delete_net();
		
		int out(double *input, double *output);
		
		int grad_input(double *input, double *output, double *grad);
		
		int hess_grad_input(double *input, double *output, double *grad, double *hess);
		
		int is_init();
		
	protected:
		// Protected Declarations
	
	private:
		// Private Declarations
		
		int layers;

		int* link;

		double** weight;

		double** bias;
		
//		double **Sij;
//		
//		double ***Gij;
//		
//		double ***Hij;
		
		int laylink_state;
		int weight_state;
		int bias_state;
//		int Sij_state;
		int net_state;
};

#endif
