#include "myrand.h"

unsigned long long int Xn = 1LL, Yn = 1LL;

int mysrand(){
	unsigned long long int N = time( NULL );
	
//	printf("%ld\n", time( NULL ));
	
	Xn = N % 11311LL;
	if(Xn==0LL) Xn = N / 11311LL;
	
	Yn = N % 43321LL;
	if(Yn==0LL) Yn = N / 43321LL;
	
	return 0;
}

double myrand(){
	
	unsigned long long int i;
	double buf;
	
//	unsigned long long int M = 1099511627776LL, K = 762939453125LL;

//	unsigned long long int M;
	
//	M = 1048576LL*1048576LL;
 	
// 	Yn = (Yn*23LL) % (100000001LL);
// 	
// 	for(i=0LL; i<Yn; i+=1LL)
// 		Xn = (762939453125LL*Xn) % 1099511627776LL;
// 	
//	buf = (double)Xn/(double)1099511627776LL;
//	i = Yn;
// (double)(rand()) / ((double)(RAND_MAX));

    Yn = (Yn*23LL+13) % (11LL);
	for(i=0LL; i<Yn+3; i+=1LL)
		Xn = (1103515245LL*Xn+12345) % 2147483648LL;
	
	buf = (double)Xn/(double)2147483648LL;

	return buf;
}
