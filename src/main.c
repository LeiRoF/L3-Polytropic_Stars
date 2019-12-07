#include "header.h"

void unitary_tests();

int main( void ){
	//unitary_tests();
	
	/* Stellar parameters */

	double Rstar = Rsun ;
	double Mstar = Msun;
	double X = 0.7; /* H abundance */
	double mu = 4. / (3. + 5. * X) ;
	double n = 3.;
	double dx = 0.01;
	
	/* Variable declaration */
	double xstar = 0.;
	double* ptxstar = &xstar;
	double dtheta_dxstar;
	double* dtheta_dx = &dtheta_dxstar;
	
	size_t nx = 1024;
	double* theta = (double *) calloc (nx, sizeof(double));
	double* x = (double *) calloc (nx, sizeof(double));
	
	
	
	// Lane_emden -> compute theta[i] & x[i]
	lane_emden(&nx, &x, &theta, n, dx);
	// nx = i_star

	// compute xstar - function (28)
	boundary_values(nx, x, theta, dx, ptxstar, dtheta_dx);
	printf("\nXstar: %lf\n", xstar);
	
	// compute dtheta_dxstar - function (30)
	compute_dtheta_dxstar(nx, x, theta, dx, ptxstar, dtheta_dx);
	printf("\ndtheta_dxstar: %lf\n", dtheta_dxstar);

	// Compute Pc - fonction (11)
	double Pc = compute_Pc(Mstar, Rstar, dtheta_dxstar);
	printf("\nPc: %lf\n", Pc);
	
	// Compute Rhoc - fonction (10)
	double rho_c = compute_rho_c(Mstar, Rstar, xstar, dtheta_dxstar);
	printf("\nrho_c: %lf\n", rho_c);
	
	// Compute Delta - fonction (32)
	double delta = compute_delta(Pc, rho_c, n, mu);
	printf("\ndelta: %lf\n", delta);
		
	// Compute Beta - Newton-Rapshon method on fonction (31)
	double beta = compute_beta(delta,1.0,1e-8);
	printf("\nbeta: %lf\n", beta);	
	
	// Export f(beta)
	for (double i = 0; i<1000; i++){
		export_beta(i/100, f(i/100, delta), TRUE);
	}
	
	export();
	return EXIT_SUCCESS;
}

