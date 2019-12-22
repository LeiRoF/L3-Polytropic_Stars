double g(double x, double y, double z, double n);

int rk4(double* x, double n, double nx){
	double k1, k2, k3, k4, h=1;
	int N=(int) nx, i=0;
	double y[N], z[N];

	y[0] = 1;
	z[0] = 0; // = dy/dx

	export_for_grace(x[0], y[0], "y_x.dat", 0); 

	for(i = 0; i < N; i++){
		k1 = g(x[i], y[i], z[i], n);
		k2 = g(x[i]+h/2, y[i]+h/2*z[i], z[i]*h/2*k1, n);
		k3 = g(x[i]+h/2,y[i]+h/2*z[i]+h*h/4*k1,z[i]+h/2*k2, n);
		k4 = g(x[i]+h,y[i]+h*z[i]+h*h/2*k2,z[i]+h*k3, n);

		y[i+1] = y[i] + h*z[i] + h*h/6*(k1+k2+k3);
		z[i+1] = z[i] + h/6*(k1+2*k2+2*k3+k4);

		export_for_grace(x[i+1], y[i+1], "y_x.dat", 1); 
	}
	printf("\nRunge-Kutta: y[i+1]= %e, z[i+1]= %e\n\n", y[i+1], z[i+1]);
}

double g(double x, double y, double z, double n){
    return (-x*pow(y,n)-2*z)/x;
}
























/*
int main(){
    //déclaration des variables
    double a, b, h;
    int k, N;
    double A1, A2, A3, A4;
    
    //conditions initiales
    a = 0;
    b = 5;
    h = 0.01/2;
    N = ((b - a) /  h)+1;
    double y[N], t[N];
    y[0] = 1;
    t[0] = a;
    
    //boucle sur k
    printf("%lf %lf\n", t[0], y[0]);
    for(k = 0; k<N-1; k++){
        A1 = f(y[k], t[k]);
        A2 = f(y[k] + h/2 * A1, t[k] + h/2);
        A3 = f(y[k] + h/2 * A2, t[k] + h/2);
        A4 = f(y[k] + h * A3, t[k] + h);
        
        t[k+1] = t[k] + h;
        y[k+1] = y[k] + h/6 * (A1 + 2*A2 + 2*A3 + A4);
        
        printf("%lf %lf\n", t[k+1], y[k+1]);
    }
    
    //printf("%lf %.15lf\n", h,  y[N-1]-exp(-b));
    
    return 0;
}
*/


