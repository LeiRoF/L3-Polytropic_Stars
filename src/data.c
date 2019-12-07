const double Msun = 1.9891e30;	// kg
const double Rsun = 6.95508e8;	// m
const double G = 6.673e-11;		// N*m^2*kg^(-2)
const double kb = 1.3806503e-23;	// J*K^(-1)
const double a = 7.565767e-16;	// J*m^{-3}*K^{-4}
const double amu = 1.66053873e-27;	// kg
const double pi = 3.14159265358979323846;

// Compute fonction (28) with (29) method
void boundary_values(size_t nx, double* const x, double* const theta, double dx, double* xstar, double* dtheta_dxstar){
	
	double xs[2], a, b, c;
	int tmp = (int) *xstar;
	
	size_t i = nx-1;
	
	// Attribution of a,b & c for solution: x+- = (-b +- sqrt(b^2 - 4ac)) / (2a)
	a= theta[i] - 2*theta[i-1] + theta[i-2];
	b= -theta[i] * (x[i-2] + x[i-1]) + 2*theta[i-1] * (x[i] + x[i-2]) - theta[i-2] * (x[i] + x[i-1]);
	c= theta[i] * x[i-2] * x[i-1] - 2*theta[i-1] * x[i] * x[i-2] + theta[i-2] * x[i] * x[i-1];
	
	xs[0] = (-b+sqrt(b*b-4*a*c))/(2*a);
	xs[1] = (-b-sqrt(b*b-4*a*c))/(2*a);

	printf("A: %lf   B: %lf   C: %lf   x[i-1]: %lf   x[i]: %lf   x+: %lf   x-: %lf\n", a, b, c, x[i-1], x[i], xs[0], xs[1]);
	
	if(xs[0] >= x[i-1] && xs[0] <= x[i])
		*xstar = xs[0];
	else if(xs[1] >= x[i-1] && xs[1] <= x[i])
		*xstar = xs[1];
	else{
		printf("[ERROR] Unable to compute Xstar. Assuming x[i-1] = Xstar\n");
		*xstar = x[i-1];
	}
		
	
}

// Comute functuon (30)
void compute_dtheta_dxstar(size_t nx, double* const x, double* const theta, double dx, double* xstar, double* dtheta_dxstar){
	size_t i = nx-1;
	*dtheta_dxstar = theta[i] * (2*(*xstar) * - x[i-2] - x[i-1])/(2*dx*dx)   -    theta[i-1] * (2*(*xstar) * - x[i] - x[i-2])/(2*dx*dx)   +   theta[i-2] * (2*(*xstar) * - x[i-1] - x[i])/(2*dx*dx);
}

// Compute Pc - fonction (11)
double compute_Pc(double Mstar, double Rstar, double dtheta_xstar){
	return G / (4 * pi) * Mstar * Mstar / (Rstar * Rstar * Rstar * Rstar * dtheta_xstar * dtheta_xstar);
}

// Compute Rhoc - fonction (10)
double compute_rho_c(double Mstar, double Rstar, double xstar, double dtheta_xstar){
	return Mstar * xstar/ (Rstar * Rstar * Rstar * 4 * pi * dtheta_xstar);
}

// Compute Delta - fonction (32)
double compute_delta(double Pc, double rho_c, double n, double mu){
	double tmp = (mu*amu/kb);
	tmp = tmp*tmp*tmp*tmp;
	return a/3*tmp*Pc*Pc*Pc
		/pow(rho_c, 3*(n+1)/n);
}


// Compute Beta - Newton-Rapshon method on fonction (31)
double f(double beta, double delta){
	return delta*beta*beta*beta*beta + beta - 1.0;
}

double df(double beta, double delta){
    return 4*delta*beta*beta*beta+1.0;
}

double compute_beta(double delta, double beta0, double epsilon){
	double beta = beta0;
	int loop = 0;
	
	do{
		beta = beta - f(beta, delta)/df(beta, delta);
		loop++;
	}while(f(beta, delta) > epsilon || f(beta, delta) < -epsilon);
	
	return beta ;
}
