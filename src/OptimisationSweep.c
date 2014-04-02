/*
 *  OptimisationSweep.c
 *  PurkinjePlasticity
 *
 *  Created by David Higgins on 17/03/2014.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

//#include "OptimisationSweep.h"
#include "Synapse.h"

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#ifdef OPTIMISATION_PROGRAM
void print_state(size_t iter, gsl_multifit_fdfsolver * s);
int print_jacobian(gsl_multifit_fdfsolver * s);

/*int main( int argc, char *argv[] ){
	Synapse *syn;
	syn = initialise_parameter_optimisation_sweep(argc, argv);
	printf("Initialistation complete\n");
	
	perform_parameter_optimisation_sim(syn);
	
    //printf("help\n");
	printf("Freeing memory and exiting...\n");
	return finalise(0, syn);
}*/


int main( int argc, char *argv[] ){
	Synapse *syn;
	struct fitting_data data_struct;
	const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmder;
	//const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder; // type of solver
	gsl_multifit_function_fdf f; // function to fit
	gsl_multifit_fdfsolver * s; // solver
	gsl_vector_view x; // initial guess
	
	int status;
	unsigned int iter = 0;
	const size_t n = 17;
	const size_t p = 8; //10;//8;//9;
	
	syn = initialise_parameter_optimisation_sweep(argc, argv);
	print_params();
	printf("Initialistation complete\n");
	
	//float test = 1;
	//test += 0.000000014901161193847656250000;
	//printf("DEBUG: float %.20f, double %.20lf\n", test, test);
	
	// tau, D, C_pf, C_cs, N_pf, theta_d, gamma_d, gamma_p
	//double x_init[8] = {185,(80./dt),0.07,0.6,0.2,0.522,2.3809e-4,7.9365e-5};
    double x_init[8] = {0,0,0,0,0,0,0,0};
	//double x_init[9] = {0,0,0,0,0,0,0,0,0};
	//double x_init[10] = {0,0,0,0,0,0,0,0,0,0};
	//double x_init[1] = {1.0};//{185.};//,(int)(80./dt),0.07,0.6,0.2,0.522,2.3809e-4,7.9365e-5};
	
	data_struct.syn = syn;
	
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	printf ("s is a '%s' solver\n", 
			gsl_multifit_fdfsolver_name (s));
	
	f.f = &cost_function;
	f.df = NULL;
	f.fdf = NULL;
	f.n = n;
	f.p = p;
	f.params = &data_struct;

	x = gsl_vector_view_array(x_init, p); // setup initial guess
	
    times_through_cost_function = 0;
	status = gsl_multifit_fdfsolver_set (s, &f, &x.vector); // setup solver
	printf("status = %s\n", gsl_strerror(status));
	
	print_jacobian(s);
	
	print_state(iter, s);
	
	printf("Beginning iteration loop\n");
	do{
		iter++;
        times_through_cost_function = 0;
		status = gsl_multifit_fdfsolver_iterate(s);
		
		printf("status = %s\n", gsl_strerror(status));
		
        print_jacobian(s);
        
		print_state(iter, s);
		
		if(status){
			printf("Breaking\n");
			break;
		}
		
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-6, 1e-6);		
	} 
	while ( (status == GSL_CONTINUE) && (iter < 10));
		
	
	if (status == GSL_SUCCESS){
		printf("Final Success\n");
	}
	else{
		printf("Final Error\n");
	}

    print_params();
	
	printf("Freeing memory and exiting...\n");
	gsl_multifit_fdfsolver_free (s); // free solver memory
	return finalise(0, syn);
}


int print_jacobian(gsl_multifit_fdfsolver * s){
	printf("Printing Jacobian:\n");
	//printf("%f\n", gsl_vector_get(s->J, 0));
	
	int status, n = 0;
	
	gsl_matrix * m = s->J;
	
	//gsl_matrix_fprintf(stdout, s->J, "%g");
	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = fprintf(stdout, "%g ", gsl_matrix_get(m, i, j))) < 0)
				return -1;
			n += status;
		}
		
		if ((status = fprintf(stdout, "\n")) < 0)
			return -1;
		n += status;
	}
	printf("End of Jacobian\n");
	
	return n;
}


void print_state(size_t iter, gsl_multifit_fdfsolver * s){
		printf("iter: %3u %g %g %g %g %g %g %g %g |f(x)| = %g\n",
			   (unsigned int)iter,
			   gsl_vector_get(s->x, 0),
			   gsl_vector_get(s->x, 1),
			   gsl_vector_get(s->x, 2),
			   gsl_vector_get(s->x, 3),
			   gsl_vector_get(s->x, 4),
			   gsl_vector_get(s->x, 5),
			   gsl_vector_get(s->x, 6),
			   gsl_vector_get(s->x, 7),
			   gsl_blas_dnrm2(s->f));
		/*printf("iter: %3u %f |f(x)| = %g\n",
		   (unsigned int)iter,
		   gsl_vector_get(s->x, 0),
		   gsl_blas_dnrm2(s->f));*/
}
#endif