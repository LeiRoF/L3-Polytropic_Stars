
/*
void export_results(double M, double R, double T, double Rho, souble P, double L, double X, double Y, double He3, double C12, double N14, double O16){

}
*/

void export(){
	FILE* file;
	file = fopen("pstar.txt", "w");
	
	fprintf(file, "%s", "Standard Solar Model (BP2004)\n");
	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "astro-ph/0402114\n");
	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "Columns in the Standard Model table (below) represent:\n");
	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "1)  Mass fraction in units of the solar mass\n");
	fprintf(file, "%s", "2)  Radius of the zone in units of one solar radius\n");
	fprintf(file, "%s", "3)  Temperature in units of deg (K)\n");
	fprintf(file, "%s", "4)  Density in units of g/cm^3\n");
	fprintf(file, "%s", "5)  Pressure in units of dyn/cm^2\n");
	fprintf(file, "%s", "6)  Luminosity fraction in units of the solar luminosity\n");
	fprintf(file, "%s", "7)  X(^1H): the hydrogen mass fraction\n");
	fprintf(file, "%s", "8)  X(^4He): the helium 4 mass fraction\n");
	fprintf(file, "%s", "9)  X(^3He): the helium 3 mass fraction\n");
	fprintf(file, "%s", "10) X(^12C): the carbon 12 mass fraction\n");
	fprintf(file, "%s", "11) X(^14N): the nitrogen 14 mass fraction\n");
	fprintf(file, "%s", "12) X(^16O): the oxygen 16 mass fraction\n");
	fprintf(file, "%s", "\n");
	fprintf(file, "%s", "The Table begins here.\n");
	fprintf(file, "%s", "\n");
    
	fclose(file);
}

void export_theta_i(double i, double theta, int write_after){
	FILE* file;
	if(write_after)
		file = fopen("theta_i.dat", "aw");	
	else
		file = fopen("theta_i.dat", "w");
	
	fprintf(file, "%lf %lf\n", i, theta);
	
	fclose(file);
}

void export_beta(double i, double theta, int write_after){
	FILE* file;
	if(write_after)
		file = fopen("beta.dat", "aw");	
	else
		file = fopen("beta.dat", "w");
	
	fprintf(file, "%lf %lf\n", i, theta);
	
	fclose(file);
}
