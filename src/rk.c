double f(double y, double t);

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

double f(double y, double t){
    double g = 0.1, RC = 1, q0 = 1;
    
    return -y - (g * y * y * q0) / RC;
}
