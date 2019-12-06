int lane_emden(size_t* nx, double** x, double** theta, double n, double dx){
	int ierr;
	size_t size = *nx;
	double dx2 = dx * dx;
	double A = 0, B = 0;
	/* TODO : you might need some definitions */
	
	export_theta_i(0, 1.0);
	export_theta_i(1, 1.0);
	
	for (size_t i = 1;; ++i){
		/* If we reach the size of the arrays , we need to re- allocate memory */
		if(i >= size - 2){
			/* We increase the size by a factor REALLOC_FATOR) */
			size_t new_size = (size_t) (size * REALLOC_FACTOR);
			ierr = resize(new_size, x);
			/* Checking i f memory could be allocated */
			if(!ierr){
				return ierr;
			}
			/* Checking i f memory could be allocated */
			ierr = resize(new_size, theta );
			if(!ierr){
				return ierr;
			}
			size = new_size;
		}
		
		/* 2nd order numerical scheme */
		
		A = 1/((1+dx/(2*(*x)[i]))*(1+dx/(2*(*x)[i])));
        B = (1+dx/(2*(*x)[i]))*(1+dx/(2*(*x)[i]));
        
        (*theta)[i + 1] = (*theta)[i] + A * (B * ((*theta)[i] - (*theta)[i - 1])-dx2 + pow((*theta)[i],(*x)[i])); /* DONE ? : you need to code the scheme here */
        
        export_theta_i(i, (*theta)[i+1]);
	
		/* We reached the outer boundary */
		
		if((*theta) [ i + 1] <= 0.){
			*nx = i + 2;
			/* We realloc again to the actual size , in order to free the useless memory */
			ierr = resize(*nx, x);
			if (!ierr ){
				return ierr;
			}
			ierr = resize(*nx, theta );
			if(!ierr){
				return ierr;
			}
		break;
		}
	}
	return TRUE;
}
