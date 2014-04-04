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
void print_final_fit(gsl_multifit_fdfsolver * s);
void print_final_params(gsl_multifit_fdfsolver *s);

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
    
    // Choose type of nonlinear solver here
	const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmder;
	//const gsl_multifit_fdfsolver_type * T = gsl_multifit_fdfsolver_lmsder;
    
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
    double x_init[8] = {1,1,1e-2,1e-2,1e-2,1e-3,1e-5,1e-5};
	//double x_init[9] = {0,0,0,0,0,0,0,0,0};
	//double x_init[10] = {0,0,0,0,0,0,0,0,0,0};
	//double x_init[1] = {1.0};//{185.};//,(int)(80./dt),0.07,0.6,0.2,0.522,2.3809e-4,7.9365e-5};
	
	data_struct.syn = syn;
	
	s = gsl_multifit_fdfsolver_alloc (T, n, p);
	printf ("s is a '%s' solver\n", 
			gsl_multifit_fdfsolver_name (s));
	
	f.f = &cost_function;
	f.df = &calculate_jacobian; //NULL;
	f.fdf = NULL;
	f.n = n;
	f.p = p;
	f.params = &data_struct;

	x = gsl_vector_view_array(x_init, p); // setup initial guess
	
    times_through_cost_function = 0;
    printf("Setting up solver \n");
	status = gsl_multifit_fdfsolver_set (s, &f, &x.vector); // setup solver
	printf("Solver setup complete, status = %s\n", gsl_strerror(status));
	
	print_jacobian(s);
	
	print_state(iter, s);
	
	printf("Beginning iteration loop\n");
	do{
		iter++;
        times_through_cost_function = 0;
        times_through_cost_function_jacobian = 0;
        printf("\nNew iteration\n");
		status = gsl_multifit_fdfsolver_iterate(s);
		printf("end of iteration update, ");
		printf("status = %s\n", gsl_strerror(status));
		
        print_jacobian(s);
        
		print_state(iter, s);
		
		if(status){
			printf("Breaking\n");
			break;
		}
		
        printf("Calculating delta test\n");
		status = gsl_multifit_test_delta(s->dx, s->x, 1e-12, 1e-12);
		printf("Delta test completed\n");
	} 
	while ( (status == GSL_CONTINUE) && (iter < 10));
		
	printf("------------------------\n");
	if (status == GSL_SUCCESS){
		printf("Final Success\n");
	}
	else{
		printf("Final Error\n");
	}

    print_params();
    print_state(iter, s);
    
    gsl_vector * temp = gsl_vector_alloc(n);
    cost_function(s->x, &data_struct, temp);
    
    print_final_params(s);
    
    print_final_fit(s);
	
	printf("Freeing memory and exiting...\n");
	gsl_multifit_fdfsolver_free (s); // free solver memory
	return finalise(0, syn);
}


int print_jacobian(gsl_multifit_fdfsolver * s){
	printf("Printing Jacobian:\n");
	//printf("%f\n", gsl_vector_get(s->J, 0));
	
	int status, n = 0;
	
	gsl_matrix * m = s->J;
	
    double column_norm[m->size2];
    for(size_t j = 0; j < m->size2; j++){
        column_norm[j] = 0.;
    }
    
	//gsl_matrix_fprintf(stdout, s->J, "%g");
	for (size_t i = 0; i < m->size1; i++) {
		for (size_t j = 0; j < m->size2; j++) {
			if ((status = fprintf(stdout, "%g ", gsl_matrix_get(m, i, j))) < 0)
				return -1;
			n += status;
            column_norm[j] += gsl_matrix_get(m, i, j);
		}
		
		if ((status = fprintf(stdout, "\n")) < 0)
			return -1;
		n += status;
	}
    
    printf("-------------------\n");
    for (size_t j = 0; j < m->size2; j++) {
        printf("%g ", column_norm[j] );
    }
    
	printf("\nEnd of Jacobian\n");
	
	return n;
}


void print_state(size_t iter, gsl_multifit_fdfsolver * s){
		/*printf("iter: %3u %g %g %g %g %g %g %g %g |f(x)| = %g\n",
			   (unsigned int)iter,
			   gsl_vector_get(s->x, 0),
			   gsl_vector_get(s->x, 1),
			   gsl_vector_get(s->x, 2),
			   gsl_vector_get(s->x, 3),
			   gsl_vector_get(s->x, 4),
			   gsl_vector_get(s->x, 5),
			   gsl_vector_get(s->x, 6),
			   gsl_vector_get(s->x, 7),
			   gsl_blas_dnrm2(s->f));*/
    printf("iter: %3u", (int)iter);
    for(int i = 0; i < s->x->size; i++){
        printf(" %g", gsl_vector_get(s->x, i));
        
    }
    printf(" |f(x)| = %g\n", gsl_blas_dnrm2(s->f));
		/*printf("iter: %3u %f |f(x)| = %g\n",
		   (unsigned int)iter,
		   gsl_vector_get(s->x, 0),
		   gsl_blas_dnrm2(s->f));*/
}


void print_final_params(gsl_multifit_fdfsolver *s){
    gsl_vector * x = s->x;
    printf("in print final params, setting params\n");
    set_optimisation_sim_params(x);
    printf("in print final params, printing params\n");
    print_params();
}

void print_final_fit(gsl_multifit_fdfsolver * s){
    gsl_vector * cost = s->f;
    
    const double objective_dw[17] = { // some of these were departures from 1 and others were absolute values
        1.04, /* Safo */
        1.108,
        1-0.16,
        1-0.28,
        1-0.208,
        1,
        1.148, /*end safo*/
		0.937, /* 2xPF  (7)*/
		1.6868986114, /* 5xPF 200Hz */
		1.2690685961, /* 33Hz */
		1.0449483297, /* 16Hz */
		1-0.325, /*Bidoret pairs  (11)*/
		1-0.35,
		1-0.31,
		1-0.15,
		1-0.08, /* end Bidoret pairs */
		1 /* depol and single pf (16) */}; // these are the target dw values
    
    printf("\tObjective\tCost\n");
    for(int i = 0; i < 17; i++){
        printf("\t %g\t %g\n", objective_dw[i], gsl_vector_get(cost, i));
    }
    printf("END\n");
}
#endif

