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

	// Calcul de Pc - fonction (11)
	double Pc = G / (4 * pi) * Mstar * Mstar / (Rstar * Rstar * dtheta_dx[xstar] * dtheta_dx[xstar]);
	
	// Calcul de Rhoc - fonction (10)
	double Pc = Mstar * xstar/ (Rstar * Rstar * Rstar * 4 * pi * dtheta_dx[xstar]);
	free(&tmp)
	
	// Calcul de Delta - fonction (32)
	double delta = a/3
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*Pc*Pc*Pc
		/pow(rho_c, 3*(n+1)/n);
		
	// Calcul de Beta - Methode de Newton-Rapshon sur la fonction (31)
	compute_beta(delta,1.0,1e-8);
	
	
	export();

	return EXIT_SUCCESS;
}







