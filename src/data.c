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
	
	// x +- = (-b +- sqrt(b^2 - 4ac)) / (2a)
	a= theta[nx] - 2*x[nx-1] + theta[nx-2];
	b= theta[nx] * (x[nx-2] + x[nx-1]) - 2*theta[nx-1] * (x[nx] + x[nx-2]) + theta[nx-2] * (x[nx] + x[nx-1]);
	c= theta[nx] * x[nx-2] * x[nx-1] - 2*theta[nx-1] * x[nx] * x[nx-2] + theta[nx-2] * x[nx] * x[nx-1];
	
	xs[0] = (-b+sqrt(b*b-4*a*c))/(2*a);
	xs[1] = (-b-sqrt(b*b-4*a*c))/(2*a);
	
	if(xs[0] >= x[nx-1] && xs[0] <= x[nx])
		*xstar = xs[0];
	else
		*xstar = xs[1];
	
}

// Comute functuon (30)
void compute_dtheta_dxstar(size_t nx, double* const x, double* const theta, double dx, double* xstar, double* dtheta_dxstar){
	*dtheta_dxstar = theta[nx] * (2*(*xstar) * - x[nx-2] - x[nx-1])/(2*dx*dx)   -    theta[nx-1] * (2*(*xstar) * - x[nx] - x[nx-2])/(2*dx*dx)   +   theta[nx-2] * (2*(*xstar) * - x[nx-1] - x[nx])/(2*dx*dx);
}

// Compute Pc - fonction (11)
double compute_Pc(double Mstar, double Rstar, double dtheta_xstar){
	return G / (4 * pi) * Mstar * Mstar / (Rstar * Rstar * dtheta_xstar * dtheta_xstar);
}

// Compute Rhoc - fonction (10)
double compute_rho_c(double Mstar, double Rstar, double xstar, double dtheta_xstar){
	return Mstar * xstar/ (Rstar * Rstar * Rstar * 4 * pi * dtheta_xstar);
}

// Compute Delta - fonction (32)
double compute_delta(double Pc, double rho_c, double n, double mu){
	return a/3
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*(mu*amu/kb)
		*Pc*Pc*Pc
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
