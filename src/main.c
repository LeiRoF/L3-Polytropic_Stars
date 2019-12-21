#include "header.h"
int main(void){
	
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
	double* dtheta_dxs = &dtheta_dxstar;
	
	size_t nx = 1024;
	double* theta = (double *) calloc (nx, sizeof(double));
	double* dtheta_dx = (double *) calloc (nx, sizeof(double));
	double* x = (double *) calloc (nx, sizeof(double));
	double* rho = (double *) calloc (nx, sizeof(double));
	double* P = (double *) calloc (nx, sizeof(double));
	double* T = (double *) calloc (nx, sizeof(double));
	double* m = (double *) calloc (nx, sizeof(double));
	
	
	// Lane_emden -> compute theta[i] & x[i]
	printf("[INFO] Determination de theta(x) ...");
	lane_emden(&nx, &x, &theta, n, dx, &dtheta_dx);
	printf(" Export dans theta_x.dat\n");
	// nx = i_star

	rk4(x, n, nx);

	// compute xstar - function (28)
	printf("[INFO] Determination de x_star ...");
	boundary_values(nx, x, theta, dx, ptxstar, dtheta_dxs);
	printf(" Xstar: %e\n", xstar);
	
	// compute dtheta_dxstar - function (30)
	printf("[INFO] Determination de dtheta(x_star)/dx ...");
	compute_dtheta_dxstar(nx, x, theta, dx, xstar, dtheta_dxs);
	printf(" dtheta(x_star)/dx: %e\n", dtheta_dxstar);

	// Compute Pc - fonction (11)
	printf("[INFO] Determination de P_c ...");
	double Pc = compute_Pc(Mstar, Rstar, n, dtheta_dxstar);
	printf(" P_c: %e\n", Pc);
	
	// Compute Rhoc - fonction (10)
	printf("[INFO] Determination de rho_c ...");
	double rho_c = compute_rho_c(Mstar, Rstar, xstar, dtheta_dxstar);
	printf(" rho_c: %e\n", rho_c);
	
	// Compute Delta - fonction (32)
	printf("[INFO] Determination de delta ...");
	double delta = compute_delta(Pc, rho_c, n, mu);
	printf(" delta: %e\n", delta);
		
	// Compute Beta - Newton-Rapshon method on fonction (31)
	printf("[INFO] Determination de beta ...");
	double beta = compute_beta(delta,1.0,1e-8);
	printf(" beta: %e\n", beta);
	
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
	
	// Compute m(x) - function (14)
	printf("[INFO] Determination de m(x) ...");
	compute_m(nx, m, Mstar, x, xstar, dtheta_dx, dtheta_dxstar);
	printf(" Export dans m_x.dat\n");
	
	// Export f(beta)
	printf("[INFO] Export de Beta.dat ...");
	for (double i = 0; i<100; i++){
		if((int) i%100 == 0)
			printf(".");
		export_for_grace(i/100, f(i/100, delta), "beta.dat", (int) i);
	}
	printf(" Ok.\n");
	
	write_results(nx, mu, xstar, beta, Pc, rho_c, x, rho, P, T, m);

	printf("Mstar= %e\n", compute_Mstar(xstar, dtheta_dxstar));

	printf("\n");	
	return EXIT_SUCCESS;
}

