double g(double x, double y, double z, double n);
double s(double x, double y, double z, double n);

double g(double x, double y, double z, double n){
    return y;
}

double s(double x, double y, double z, double n){
	if(fabs(x) < 1e-5){
		return -1/3; //DL quand x -> 0
	}
    return -pow(y,n) - 2*z/x;
}


int rk4(double* x, double n, double nx){
	double k1, k2, k3, k4, h=1;
	int N=(int) nx, i=0;
	double y[N], z[N];

	y[0] = 1;
	z[0] = 0; // = dy/dx

	export_for_grace(x[0], y[0], "y_x.dat", 0); 

	for(i = 0; i < N; i++){
		k1 = s(x[i], y[i], z[i], n);
		k2 = s(x[i]+h/2, y[i]+h/2*z[i], z[i]*h/2*k1, n);
		k3 = s(x[i]+h/2,y[i]+h/2*z[i]+h*h/4*k1,z[i]+h/2*k2, n);
		k4 = s(x[i]+h,y[i]+h*z[i]+h*h/2*k2,z[i]+h*k3, n);

		y[i+1] = y[i] + h*z[i] + h*h/6*(k1+k2+k3);
		z[i+1] = z[i] + h/6*(k1+2*k2+2*k3+k4);

		export_for_grace(x[i+1], y[i+1], "y_x.dat", 1);
	}
	printf("\nRunge-Kutta: y[i+1]= %e, z[i+1]= %e\n\n", y[i+1], z[i+1]);
}

int rk4_bis(double* x, double n, double nx){
	double k11, k12, k13, k14, k21, k22, k23, k24, h=1;
	int N=(int) nx, i=0;
	double y[N], z[N];

	y[0] = 1;
	z[0] = 0; // = dy/dx

	export_for_grace(x[0], y[0], "yb_x.dat", 0); 

	for(i = 0; i < N; i++){
		k11 = h*g(x[i], y[i], z[i], n);
		k21 = h*s(x[i], y[i], z[i], n);
		
		k12 = h*g(x[i]+h/2, y[i]+k11/2, z[i]*k21/2, n);
		k22 = h*s(x[i]+h/2, y[i]+k11/2, z[i]*k21/2, n);
		
		k13 = h*g(x[i]+h/2, y[i]+k12/2, z[i]+k22/2, n);
		k23 = h*s(x[i]+h/2, y[i]+k12/2, z[i]+k22/2, n);
		
		k14 = h*g(x[i]+h, y[i]+k13,z[i]+k23, n);
		k24 = h*s(x[i]+h, y[i]+k13,z[i]+k23, n);

		y[i+1] = y[i] + (k11+2*k12+2*k13+k14)/6;
		z[i+1] = z[i] + (k21+2*k22+2*k23+k24)/6;

		export_for_grace(x[i+1], y[i+1], "yb_x.dat", 1);
	}
	printf("\nRunge-Kutta: y[i+1]= %e, z[i+1]= %e\n\n", y[i+1], z[i+1]);
}
