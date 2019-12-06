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
	double xstar = 0.;
	double dx = 0.1;
	
	size_t nx = 1024;
	double* theta = (double *) calloc (nx, sizeof(double));
	double* dtheta_dx = (double *) calloc (nx, sizeof(double));
	double* x = (double *) calloc (nx, sizeof(double));

	// Compute Pc - fonction (11)
	double Pc = compute_Pc(Mstar, Rstar, dtheta_dx[(int) xstar]);
	
	// Compute Rhoc - fonction (10)
	double rho_c = compute_rho_c(Mstar, Rstar, xstar, dtheta_dx[(int) xstar]);
	
	// Compute Delta - fonction (32)
	double delta = compute_delta(Pc, rho_c, n, mu);
		
	// Compute Beta - Newton-Rapshon method on fonction (31)
	compute_beta(delta,1.0,1e-8);
	
	// Lane_emden -> calcul theta[i] et x[i]
	lane_emden(&nx, &x, &theta, n, dx);
	
	// calcul xstar
	//boundary_values(nx, &&x, &&theta, dx, &xstar, &&dtheta_dx);
	
	export();

	return EXIT_SUCCESS;
}


void unitary_tests(){
	printf("\n[DEBUG] *** TESTS UNITAIRES ***");
	printf("\n |");
	printf("\n | Test unitaire de: double compute_Pc(double Mstar, double Rstar, double dtheta_xstar)");
	printf("\n | Valeur: %lf", compute_Pc(1.9891e30, 6.95508e8, 1.0));
	printf("\n |");
	printf("\n | Test unitaire de: double compute_rho_c(double Mstar, double Rstar, double xstar double dtheta_xstar)");
	printf("\n | Valeur: %lf", compute_rho_c(1.9891e30, 6.95508e8, 100.0, 1.0));
	printf("\n |");
	printf("\n | Test unitaire de: double compute_delta(double Pc, double rho_c, double n, double mu)");
	printf("\n | Valeur: %lf", compute_rho_c(40e10, 47e3, 3.0, 0.5));
	printf("\n |");
	printf("\n | Test unitaire de: double compute_beta(double delta, double beta0, double epsilon)");
	printf("\n | Valeur: %lf", compute_beta(1.0,1.0,1e-8));
	printf("\n |");
	printf("\n | Test unitaire de: void export_theta_i(size_t i, double theta)");
	export_theta_i(0, 42, FALSE);
	export_theta_i(1, 666, TRUE);
	export_theta_i(2, 3.14159, TRUE);
	
	printf("\n |");
	printf("\n[DEBUG] *** FIN DES TESTS UNITAIRES ***");
}


