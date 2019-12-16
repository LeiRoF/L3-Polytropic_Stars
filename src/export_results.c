
void export_for_grace(double X, double Y, char name[], int write_after){
	FILE* file;
	if(write_after)
		file = fopen(name, "aw");	
	else
		file = fopen(name, "w");
	
	fprintf(file, "%lf %lf\n", X, Y);
	
	fclose(file);
}


double write_results(size_t nx, double mu, double xstar, double beta, double Pc, double rhoc, double* x, double* density, double* pressure, double* temperature, double* mass){

 /*************************************/
  /********** Writing results **********/
  /*************************************/

  /* opening file */
 
  FILE* file;
  file = fopen("pstar.txt", "w");

  /* printing some results to stdout */
 
  printf("mu    = %.*e\n", DBL_DIG, mu);
  printf("xstar = %.*e\n", DBL_DIG, xstar);
  printf("beta  = %.*e\n", DBL_DIG, beta);
  printf("P_c   = %.*e [Pa]\n", DBL_DIG, Pc);
  printf("rho_c = %.*e [kg/m^3]\n", DBL_DIG, rhoc);

  /* normalizing x */

  for(size_t ix = 0; ix < nx; ++ix) {
    x[ix] /= xstar;
  }

  for(size_t ix = 0; ix < nx - 1; ++ix) {
    /* We write cell wall values */
   
    fprintf(file, "%.*e %.*e %.*e %.*e %.*e\n",
        DBL_DIG, .5 * (x[ix] + x[ix + 1]),
        DBL_DIG, .5 * (density[ix] + density[ix + 1]),
        DBL_DIG, .5 * (pressure[ix] + pressure[ix+1]),
        DBL_DIG, .5 * (temperature[ix] + temperature[ix+1]),
        DBL_DIG, mass[ix]
        );
  }
}

