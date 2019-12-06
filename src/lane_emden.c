int lane_emden(size_t* nx, double** x, double** theta, double n, double dx){
	int ierr;
	size_t size = *nx;
	double dx2 = dx * dx;
	/* TODO : you might need some definitions */
	for (size_t k = 1;; ++k){
		/* If we reach the size of the arrays , we need to re- allocate memory */
		if(k >= size - 2){
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
		
		(*theta)[k + 1] = compute_next_theta((*theta)[k], (*theta)[k-1], dx, *x[k], n);/* DONE ? : you need to code the scheme here /*
	
		/* We reached the outer boundary */
		
		if((*theta) [ k + 1] <= 0.){
			*nx = k + 2;
			/* We realloc again to the actual size , in order to free the useless memory */
			ierr = resize(*nx, x);
			if (!ierr ){
				return ierr;
			}
			ierr = resize(*nx, theta );
			if(!ierr){
				return ierr;
			}
		break ;
		}
	}
	return TRUE;
} 
