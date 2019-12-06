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

/********/
/* BETA */
/********/

double f(double beta, double delta){
	return delta*beta*beta*beta*beta + beta - 1.0;
}

double df(double beta, double delta){
    return 4*delta*beta*beta*beta+1.0;
}

double compute_beta(double delta, double beta0, double epsilon){
	double beta = beta0;
	int loop = 0;

	/* DONE : code the Newton-Raphson method */
	do{
		beta = beta - f(beta, delta)/df(beta, delta);
		printf("[DEBUG] %d B %.15lf f(B) %le df(B) %le\n", loop, beta, f(beta, delta), df(beta, delta));
		loop++;
	}while(f(beta, delta) > epsilon || f(beta, delta) < -epsilon);
		printf("Valeur trouvée %lf\n", beta);
		return beta ;
}

