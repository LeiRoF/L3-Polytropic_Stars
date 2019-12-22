#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#define REALLOC_FACTOR 1.25
#define TRUE 1
#define FALSE 0

const double Msun = 1.9891e30;	// kg
const double Rsun = 6.95508e8;	// m
const double G = 6.673e-11;		// N*m^2*kg^(-2)
const double kb = 1.3806503e-23;	// J*K^(-1)
const double a = 7.565767e-16;	// J*m^{-3}*K^{-4}
const double amu = 1.66053873e-27;	// kg
const double pi = 3.14159265358979323846;

const double h = 6.62607015e-34;
const double c = 299792458;
const double mu_e = 2;
const double Na = 6.02214086e23;
const double m_p = 1e-3/Na;

#include "export_results.c"
#include "compute.c"
#include "utils.c"
#include "rk.c"
#include "lane_emden.c"
