#include "header.h"

int main( void ){
	/* Stellar parameters */

	double Rstar = Rsun ;
	double Mstar = Msun;
	double X = 0.7; /* H abundance */
	double mu = 4. / (3. + 5. * X) ;
	double n = 3.;
	double xstar = 0;
	
	size_t nx = 1024;
	double* theta = (double *) calloc (nx, sizeof(double));
	double* dtheta_dx = (double *) calloc (nx, sizeof(double));
	double* x = (double *) calloc (nx, sizeof(double));

	// Compute Pc - fonction (11)
	int tmp = (int) xstar;
	double Pc = G / (4 * pi) * Mstar * Mstar / (Rstar * Rstar * dtheta_dx[tmp] * dtheta_dx[tmp]);
	
	// Compute Rhoc - fonction (10)
	double rho_c = Mstar * xstar/ (Rstar * Rstar * Rstar * 4 * pi * dtheta_dx[tmp]);
	
	// Compute Delta - fonction (32)
	double delta = a/3
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*Pc*Pc*Pc
		/pow(rho_c, 3*(n+1)/n);
		
	// Compute Beta - Newton-Rapshon method on fonction (31)
	compute_beta(delta,1.0,1e-8);
	
	
	export();

	return EXIT_SUCCESS;
}







