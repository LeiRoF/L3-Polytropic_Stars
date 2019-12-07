

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

	//printf("A: %lf   B: %lf   C: %lf   x[i-1]: %lf   x[i]: %lf   x+: %lf   x-: %lf\n", a, b, c, x[i-1], x[i], xs[0], xs[1]);
	
	if(xs[0] >= x[i-1] && xs[0] <= x[i])
		*xstar = xs[0];
	else if(xs[1] >= x[i-1] && xs[1] <= x[i])
		*xstar = xs[1];
	else{
		printf("\n[ERROR] Unable to compute Xstar. Assuming x[i-1] = Xstar\n");
		*xstar = x[i-1];
	}
		
	
}

// Comute functuon (30)
void compute_dtheta_dxstar(size_t nx, double* const x, double* const theta, double dx, double xstar, double* dtheta_dxstar){
	size_t i = nx-1;
	
	//printf("theta[i] %lf   theta[i-1] %lf   theta[i-2] %lf   x[i] %lf   x[i-1] %lf   x[i-2] %lf   Xs %lf   dx %lf", theta[i] , theta[i-1] , theta[i-2] , x[i] , x[i-1] , x[i-2] , xstar , dx);
	*dtheta_dxstar = theta[i] * (2*xstar - x[i-2] - x[i-1])/(2*dx*dx)   -    theta[i-1] * (2*xstar - x[i] - x[i-2])/(2*dx*dx)   +   theta[i-2] * (2*xstar - x[i-1] - x[i])/(2*dx*dx);
}

// Compute Pc - fonction (11)
double compute_Pc(double Mstar, double Rstar, double dtheta_xstar){
	return G / (4 * pi) * Mstar * Mstar / (Rstar * Rstar * Rstar * Rstar * dtheta_xstar * dtheta_xstar);
}

// Compute Rhoc - fonction (10)
double compute_rho_c(double Mstar, double Rstar, double xstar, double dtheta_xstar){
	return Mstar * xstar/ (Rstar * Rstar * Rstar * 4 * pi * fabs(dtheta_xstar));
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

// Compute rho(x) - function (12)
void compute_rho(size_t nx, double* rho, double rho_c, double* theta, double* x, double n){
	for(size_t i = 0; i<nx; i++){
		if((int) i%100 == 0)
			printf(".");
		rho[i] = rho_c * pow(theta[i], n);
		export_for_grace(x[i]/x[nx-1], rho[i], "rho_x.dat", (int) i);
	}
}

// Compute P(x) - function (13)
void compute_P(size_t nx, double* P, double Pc, double* theta, double* x, double n){
	for(size_t i = 0; i<nx; i++){
		if((int) i%100 == 0)
			printf(".");
		P[i] = Pc * pow(theta[i], n+1);
		export_for_grace(x[i]/x[nx-1], P[i], "P_x.dat", (int) i);
	}
}

// Compute T(x) - function (16)
void compute_T(size_t nx, double* T, double mu, double beta, double Pc, double rho_c, double* theta, double* x){
	for(size_t i = 0; i<nx; i++){
		if((int) i%100 == 0)
			printf(".");
		T[i] = mu*amu*beta/kb * Pc/rho_c * theta[i];
		export_for_grace(x[i]/x[nx-1], T[i], "T_x.dat", (int) i);
	}
}

