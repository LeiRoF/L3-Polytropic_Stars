/* Resize */

int resize(size_t size, double** data){
	/* This method is used to resize an array and handling properly whenever we cannot allocate the memory */
	double* new_data = realloc (*data, size * sizeof(double));
	if(new_data){
		*data = new_data;
	} else {
		return FALSE;
	}
	return TRUE;
}
