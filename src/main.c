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
	double dtheta_dxstar = 0.;
	double* dtheta_dx = &dtheta_dxstar;
	
	size_t nx = 1024;
	double* theta = (double *) calloc (nx, sizeof(double));
	double* x = (double *) calloc (nx, sizeof(double));
	double* rho = (double *) calloc (nx, sizeof(double));
	double* P = (double *) calloc (nx, sizeof(double));
	double* T = (double *) calloc (nx, sizeof(double));
	
	
	// Lane_emden -> compute theta[i] & x[i]
	printf("[INFO] Determination de theta(x) ...");
	lane_emden(&nx, &x, &theta, n, dx);
	printf(" Export dans theta_x.dat\n");
	// nx = i_star

	// compute xstar - function (28)
	printf("[INFO] Determination de x_star ...");
	boundary_values(nx, x, theta, dx, ptxstar, dtheta_dx);
	printf(" Xstar: %lf\n", xstar);
	
	// compute dtheta_dxstar - function (30)
	printf("[INFO] Determination de dtheta(x_star)/dx ...");
	compute_dtheta_dxstar(nx, x, theta, dx, xstar, dtheta_dx);
	printf(" dtheta(x_star)/dx: %lf\n", dtheta_dxstar);

	// Compute Pc - fonction (11)
	printf("[INFO] Determination de  P_c ...");
	double Pc = compute_Pc(Mstar, Rstar, dtheta_dxstar);
	printf(" P_c: %lf\n", Pc);
	
	// Compute Rhoc - fonction (10)
	printf("[INFO] Determination de rho_c ...");
	double rho_c = compute_rho_c(Mstar, Rstar, xstar, dtheta_dxstar);
	printf(" rho_c: %lf\n", rho_c);
	
	// Compute Delta - fonction (32)
	printf("[INFO] Determination de delta ...");
	double delta = compute_delta(Pc, rho_c, n, mu);
	printf(" delta: %lf\n", delta);
		
	// Compute Beta - Newton-Rapshon method on fonction (31)
	printf("[INFO] Determination de beta ...");
	double beta = compute_beta(delta,1.0,1e-8);
	printf(" beta: %lf\n", beta);
	
	// Compute rho(x) - function (12)
	printf("[INFO] Determination de rho(x) ...");
	compute_rho(nx, rho, rho_c, theta, x, n);
	printf(" Export dans rho_x.dat\n");
	
	// Compute P(x) - function (13)
	printf("[INFO] Determination de P(x) ...");
	compute_P(nx, P, Pc, theta, x, n);
	printf(" Export dans P_x.dat\n");
	
	// Compute T(x) - function (16)
	printf("[INFO] Determination de T(x) ...");
	compute_T(nx, T, mu, beta, Pc, rho_c, theta, x);
	printf(" Export dans T_x.dat\n");
	
	// Export f(beta)
	printf("[INFO] Export de Beta.dat ...");
	for (double i = 0; i<1000; i++){
		if((int) i%100 == 0)
			printf(".");
		export_for_grace(i/100, f(i/100, delta), "beta.dat", (int) i);
	}
	
	export();
	return EXIT_SUCCESS;
}

