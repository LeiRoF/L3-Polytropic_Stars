const double Msun = 1.9891e30;	// kg
const double Rsun = 6.95508e8;	// m
const double G = 6.673e-11;		// N*m^2*kg^(-2)
const double kb = 1.3806503e-23;	// J*K^(-1)
const double a = 7.565767e-16;	// J*m^{-3}*K^{-4}
const double amu = 1.66053873e-27;	// kg
const double pi = 3.14159265358979323846;

void boundary_values(size_t nx, double* const x, double* const theta, double dx, double* xstar, double* dtheta_dxstar){
	
	double xs[2], a, b, c;
	int tmp = (int) *xstar;
	a= theta[tmp] - 2*x[tmp-1] + theta[tmp-2];
	b= theta[tmp] * (x[tmp-2] + x[tmp-1]) - 2*theta[tmp-1] * (x[tmp] + x[tmp-2]) + theta[tmp-2] * (x[tmp] + x[tmp-1]);
	c= theta[tmp] * x[tmp-2] * x[tmp-1] - 2*theta[tmp-1] * x[tmp] * x[tmp-2] + theta[tmp-2] * x[tmp] * x[tmp-1];
	xs[0] = (-b+sqrt(b*b-4*a*c))/(2*a);
	xs[0] = (-b-sqrt(b*b-4*a*c))/(2*a);
	
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


/*********/
/* THETA */
/*********/
/*
double compute_theta_i(double i, double dx, double x, double n){
	double theta = 1;
	double theta_old = 0;
	double j=0;

	do{
		theta += 1/(
				(1-dx/(2*x))*(1-dx/(2*x)))
			*(
				(1-dx/(2*x))*(1-dx/(2*x))
				*(theta - theta_old)
				-pow(theta,n)
				*x*x
			);
		j++;
	}while(j<i);
	return theta;
}

double compute_next_theta(double theta, double theta_old, double dx, double x, double n){
	theta += 1/(
			(1-dx/(2*x))*(1-dx/(2*x)))
		*(
			(1-dx/(2*x))*(1-dx/(2*x))
			*(theta - theta_old)
			-pow(theta,n)
			*x*x
		);
	return theta;
}*/
