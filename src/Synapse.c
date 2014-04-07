#include "GeneralIncludes.h"
#include "Synapse.h"
#include "DataTools.h"


#include "NumericalTools.h"
#include "SpikeTrains.h"


#ifdef PR_OPTIMISATION_PROGRAM
//TODO: These are the Polak-Ribiere functions
void pr_fdf(const gsl_vector *x, void * params, double * f, gsl_vector *df){
	*f = cost_function(x, params);
	calculate_gradient(x, params, df);
}

void calculate_gradient(const gsl_vector * x_orig, void * data, gsl_vector * J){ /* PR_OPTIMISATION_PROGRAM */
	printf("Numerically calculating gradient vector\n");
    //int rows = (int)(J->size);
    int cols = (int)(J->size); // number of parameters
    printf("DEBUG: size %d\n", cols);
    
    double x_local;
    double df_dx;
	double basic_cost;
	double new_cost;
	double diff;
	double diff_norm;
	
    gsl_vector * new_x = gsl_vector_alloc(cols);
    //gsl_vector * f_base = gsl_vector_alloc(rows);
    //gsl_vector * f_delta = gsl_vector_alloc(rows);
    
    const double dx[8] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6};
    //const double dx[9] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6, 0.5};
    //const double dx[10] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6, 0.5, 1};
    
    // make copy of original x_values to which we will add dx
    gsl_vector_memcpy(new_x, x_orig);
    
    // setup baseline f, for calculation of df/dx
    basic_cost = cost_function(new_x, data);
	times_through_cost_function_jacobian++;
    
    for (int j = 0; j < cols; j++) { // loop over parameters (cols)
        // save original value of x, was previously doing this inline but this looks faster
        x_local = gsl_vector_get(x_orig, j);
        
        // increment element j of x by dx
        gsl_vector_set(new_x, j, ( x_local + dx[j] ) );
        
        // calculate cost_function
        new_cost = cost_function(new_x, data);
        
		diff = basic_cost - new_cost;
        diff_norm = diff/dx[j];
        //printf("old %lf, new %lf, diff %lf, diff_norm %lf\n", old, new, diff, diff_norm);
            
        // calculate df/dx
        df_dx = (basic_cost - new_cost) / dx[j];
            
		// save in Jacobian matrix
        gsl_vector_set(J, j, df_dx);
            
		printf("j %d, old %lf, new %lf, diff %lf, diff_norm %lf, df_dx %lf, matrix el %lf\n", j, basic_cost, new_cost, diff, diff_norm, df_dx, gsl_vector_get(J, j));
        
        // reset new_x to original guess values
        //gsl_vector_memcpy(&new_x, x);
        gsl_vector_set(new_x, j, x_local);
        
        times_through_cost_function_jacobian++;
    }
    
    // Zero out the effects of Safo7, just to see if it helps numerical solution
    /*for (int j = 0; j < cols; j++){
	 gsl_matrix_set(J, 6, j, 0);
	 }*/
    
    // Are these calls to free necessary? Surely the function stack will be completely destroyed. (But perhaps the vectors reside in the GSL library space)
    gsl_vector_free(new_x);
    //gsl_vector_free(f_base);
    //gsl_vector_free(f_delta);
    printf("Finished numerically calculating gradient\n");
    //return GSL_SUCCESS;
} /* PR_OPTIMISATION_PROGRAM */

double cost_function(const gsl_vector * x, void * data){ /* PR_OPTIMISATION_PROGRAM */
	int i;
	Synapse * syn;
	syn = ((struct fitting_data *) data)->syn;
	int signal_change = 0;
    double final_cost = 0;
	
    const double cost_coeffs[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    //const double cost_coeffs[17] = {3,1,1,3,1,1,3,1,1,1,1,3,3,1,1,1,1}; // version used for first set of figs: priority max dw on Safo and Bidoret
	//const double cost_coeffs[17] = {1,1,1,1,1,1,1,1,3,3,1,3,3,1,1,1,1}; // prioritise PF and Bidoret max dw
	//const double cost_coeffs[17] = {3,1,1,3,1,1,3,1,3,3,1,1,1,1,1,1,1}; // prioritise PF and Safo max dw
	//const double cost_coeffs[17] = {2000000,3000000,1000000,3000000,1000000,1000000,3000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000};
	
    /*const double objective_error_bars[17] = {0.0600*0.0600,
	 0.0440*0.0440,
	 0.0400*0.0400,
	 0.0400*0.0400,
	 0.0800*0.0800,
	 0.0800*0.0800,
	 0.0700*0.0700,
	 0.0512*0.0512, */ /* 2xPF */
	/* 0.1120*0.1120,
	 0.0733*0.0733,
	 0.1546*0.1546,
	 0.0500*0.0500,
	 0.0500*0.0500,
	 0.1000*0.1000,
	 0.1200*0.1200,
	 0.0300*0.0300,
	 0.0600*0.0600};*/
    
    const double objective_error_bars[17] = {0.0600,
        0.0440,
        0.0400,
        0.0400,
        0.0800,
        0.0800,
        0.0700,
        0.0512, /* 2xPF */
        0.1120,
        0.0733,
        0.1546,
        0.0500,
        0.0500,
        0.1000,
        0.1200,
        0.0300,
        0.0600};
    
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
	
	
	double simulated_dw[17]; // these will be the values we acutally obtain
	double cost[17];
	
    printf("Times through cost function %d, from gradient function %d\n", times_through_cost_function, times_through_cost_function_jacobian);
    times_through_cost_function++;
    
	// set new params based on what gsl sends
	set_optimisation_sim_params(x);
	//print_params();
    
	// update sim
	printf("DEBUG: performing sim...\n ");
	perform_parameter_optimisation_sim(syn);
	printf("done\n");
	
	// calculate cost based on (sim weight change - experimental weight change)^2 / sigma^2
	printf("Cost calculation: \n\t cost \t objective \t simulation \t rho final \t weighted cost \n");
	double norm = 0;
	for(i = 0; i < no_synapses; i++){
		simulated_dw[i] = syn[i].rho[simulation_duration-1] / 0.5; // divide by 0.5 to normalise
		//cost[i] = ( ( (objective_dw[i] - simulated_dw[i]) * (objective_dw[i] - simulated_dw[i]) ) / objective_error_bars[i] );
		cost[i] = ( ( (objective_dw[i] - simulated_dw[i]) ) / objective_error_bars[i] );
        printf("\t %f\t %f\t %f\t %f\t ", cost[i], objective_dw[i], simulated_dw[i], syn[i].rho[simulation_duration-1]);
		//norm += cost[i]; // * cost[i];
        norm += cost[i] * cost[i];
        cost[i] *= cost_coeffs[i];
        printf("%f \n", cost[i]);
        //gsl_vector_set(f, i, cost[i]);
        if (times_through_cost_function > 1){
            if(fabs(old_simulated_dw[i] - simulated_dw[i]) > 0.0000001)
            {
                signal_change++;
                printf("change detected, i %d, old value %0.15f, new value %0.15f\n", i, old_simulated_dw[i], simulated_dw[i]);
                old_simulated_dw[i] = simulated_dw[i];
            }
        }
        else{
            old_simulated_dw[i] = simulated_dw[i];
        }
	}
	final_cost = norm;
	printf("chi-squared %f\n", norm);
	printf("number of elements of f which changed %d\n", signal_change);
	printf("\n");
	
	//calculate_summary_data(syn); // not needed for parameter optimisation, just nice for debugging
	
	
	return final_cost;
}
#endif /* PR_OPTIMISATION_PROGRAM */


#ifdef LM_OPTIMISATION_PROGRAM
//TODO: These are the Levenberg-Marquardt functions
int calculate_jacobian(const gsl_vector * x_orig, void * data, gsl_matrix * J){ /* LM_OPTIMISATION_PROGRAM */
    printf("Numerically calculating Jacobian\n");
    int rows = (int)(J->size1);
    int cols = (int)(J->size2);
    printf("DEBUG: size1 %d, size2 %d\n", rows, cols);
    
    double x_local;
    double df_dx;
    gsl_vector * new_x = gsl_vector_alloc(cols);
    gsl_vector * f_base = gsl_vector_alloc(rows);
    gsl_vector * f_delta = gsl_vector_alloc(rows);
    
    const double dx[8] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6};
    //const double dx[9] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6, 0.5};
    //const double dx[10] = {1, 1, 0.001, 0.001, 0.001, 0.001, 1e-6, 1e-6, 0.5, 1};
    
    // make copy of original x_values to which we will add dx
    gsl_vector_memcpy(new_x, x_orig);
    
    // setup baseline f, for calculation of df/dx
    cost_function(new_x, data, f_base);
    
    for (int j = 0; j < cols; j++) { // loop over parameters (cols)
        // save original value of x, was previously doing this inline but this looks faster
        x_local = gsl_vector_get(x_orig, j);
        
        // increment element j of x by dx
        gsl_vector_set(new_x, j, ( x_local + dx[j] ) );
        
        // calculate cost_function
        cost_function(new_x, data, f_delta);
        
        for (int i = 0; i < rows; i++) {
            double old = gsl_vector_get(f_base, i);
            double new = gsl_vector_get(f_delta, i);
            double diff = old - new;
            double diff_norm = diff/dx[i];
            //printf("old %lf, new %lf, diff %lf, diff_norm %lf\n", old, new, diff, diff_norm);
            
            // calculate df/dx
            df_dx = (gsl_vector_get(f_base, i) - gsl_vector_get(f_delta, i)) / dx[i];
            
            // save in Jacobian matrix
            gsl_matrix_set(J, i, j, df_dx);
            
            printf("i %d, j %d, old %lf, new %lf, diff %lf, diff_norm %lf, df_dx %lf, matrix el %lf\n", i, j, old, new, diff, diff_norm, df_dx, gsl_matrix_get(J, i, j));
        }
        
        // reset new_x to original guess values
        //gsl_vector_memcpy(&new_x, x);
        gsl_vector_set(new_x, j, x_local);
        
        times_through_cost_function_jacobian++;
    }
    
    // Zero out the effects of Safo7, just to see if it helps numerical solution
    /*for (int j = 0; j < cols; j++){
	 gsl_matrix_set(J, 6, j, 0);
	 }*/
    
    // Are these calls to free necessary? Surely the function stack will be completely destroyed. (But perhaps the vectors reside in the GSL library space)
    gsl_vector_free(new_x);
    gsl_vector_free(f_base);
    gsl_vector_free(f_delta);
    printf("Finished numerically calculating Jacobian\n");
    return GSL_SUCCESS;
} /* LM_OPTIMISATION_PROGRAM */

//float* cost_function(float *cost, Synapse *syn){
int cost_function(const gsl_vector * x, void * data, gsl_vector * f){ /* LM_OPTIMISATION_PROGRAM */
	int i;
	Synapse * syn;
	syn = ((struct fitting_data *) data)->syn;
	int signal_change = 0;
    
    const double cost_coeffs[17] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    //const double cost_coeffs[17] = {3,1,1,3,1,1,3,1,1,1,1,3,3,1,1,1,1}; // version used for first set of figs: priority max dw on Safo and Bidoret
	//const double cost_coeffs[17] = {1,1,1,1,1,1,1,1,3,3,1,3,3,1,1,1,1}; // prioritise PF and Bidoret max dw
	//const double cost_coeffs[17] = {3,1,1,3,1,1,3,1,3,3,1,1,1,1,1,1,1}; // prioritise PF and Safo max dw
	//const double cost_coeffs[17] = {2000000,3000000,1000000,3000000,1000000,1000000,3000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000};
	
    /*const double objective_error_bars[17] = {0.0600*0.0600,
	 0.0440*0.0440,
	 0.0400*0.0400,
	 0.0400*0.0400,
	 0.0800*0.0800,
	 0.0800*0.0800,
	 0.0700*0.0700,
	 0.0512*0.0512, */ /* 2xPF */
	/* 0.1120*0.1120,
	 0.0733*0.0733,
	 0.1546*0.1546,
	 0.0500*0.0500,
	 0.0500*0.0500,
	 0.1000*0.1000,
	 0.1200*0.1200,
	 0.0300*0.0300,
	 0.0600*0.0600};*/
    
    const double objective_error_bars[17] = {0.0600,
        0.0440,
        0.0400,
        0.0400,
        0.0800,
        0.0800,
        0.0700,
        0.0512, /* 2xPF */
        0.1120,
        0.0733,
        0.1546,
        0.0500,
        0.0500,
        0.1000,
        0.1200,
        0.0300,
        0.0600};
    
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
	
	
	double simulated_dw[17]; // these will be the values we acutally obtain
	double cost[17];
	
    printf("Times through cost function %d, from jacobian %d\n", times_through_cost_function, times_through_cost_function_jacobian);
    times_through_cost_function++;
    
	// set new params based on what gsl sends
	set_optimisation_sim_params(x);
	//print_params();
    
	// update sim
	printf("DEBUG: performing sim...\n ");
	perform_parameter_optimisation_sim(syn);
	printf("done\n");
	
	// calculate cost based on (sim weight change - experimental weight change)^2 / sigma^2
	printf("Cost calculation: \n\t cost \t objective \t simulation \t rho final \t weighted cost \n");
	double norm = 0;
	for(i = 0; i < no_synapses; i++){
		simulated_dw[i] = syn[i].rho[simulation_duration-1] / 0.5; // divide by 0.5 to normalise
		//cost[i] = ( ( (objective_dw[i] - simulated_dw[i]) * (objective_dw[i] - simulated_dw[i]) ) / objective_error_bars[i] );
		cost[i] = ( ( (objective_dw[i] - simulated_dw[i]) ) / objective_error_bars[i] );
        printf("\t %f\t %f\t %f\t %f\t ", cost[i], objective_dw[i], simulated_dw[i], syn[i].rho[simulation_duration-1]);
		//norm += cost[i]; // * cost[i];
        norm += cost[i] * cost[i];
        cost[i] *= cost_coeffs[i];
        printf("%f \n", cost[i]);
        gsl_vector_set(f, i, cost[i]);
        if (times_through_cost_function > 1){
            if(fabs(old_simulated_dw[i] - simulated_dw[i]) > 0.0000001)
            {
                signal_change++;
                printf("change detected, i %d, old value %0.15f, new value %0.15f\n", i, old_simulated_dw[i], simulated_dw[i]);
                old_simulated_dw[i] = simulated_dw[i];
            }
        }
        else{
            old_simulated_dw[i] = simulated_dw[i];
        }
	}
	printf("chi-squared %f\n", norm);
	printf("number of elements of f which changed %d\n", signal_change);
	printf("\n");
	
	//calculate_summary_data(syn); // not needed for parameter optimisation, just nice for debugging
	
	
	return GSL_SUCCESS;
}
#endif /* LM_OPTIMISATION_PROGRAM */


//TODO: end of optimisation method specific functions
int perform_parameter_optimisation_sim(Synapse *syn){
	int t, i;
	
	// It looks like we can avoid resetting rho(0), c(0), etc because of forward method
	
	siT = 0; // siT is a global which detects the time within each update function
	
	print_params();
    printf("DEBUG: simulation_duration %ld\n", simulation_duration);
    
	// Loop over discrete time steps up to simulation_duration
	//for (t = siT; t < (simulation_duration-1); t++){
    for (i = 0; i < no_synapses; i++){
        // Update each synapse
		
        siT = 0;
    
        //for (i = 0; i < no_synapses; i++)
        for (t = 0; t < (simulation_duration-1); t++)
		{
			//i = 12;
			//printf("syn(%d) t %d\n", i, t);
			//updatePreSynapticVoltageTrace(&syn[i]);
			updatePreSynapticNOConcentration(&syn[i]);
			updateCalciumConcentration(&syn[i]);
			updateSynapticEfficacy(&syn[i]);
			//if (i == 12)
				//printf("t: %d, c: %f, rho: %f, NO: %f\n", siT, syn[i].c[siT-time_of_last_save], syn[i].rho[siT], syn[i].NO_pre[siT]);
            siT++;
		}
        //printf("DEBUG: RHO %f\n", syn[0].rho[siT]);
        /*if(siT == 100000){
            printf("DEBUG: RHO %f\n", syn[0].rho[siT]);
            
        }*/
		//siT++;
	}
	
	/*char outfile[FILE_NAME_LENGTH];
	sprintf(outfile, outfilepattern, syn[i].ID);
	printf("writing...%s\n", outfile);
	saveSynapseOutputFile(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iCaSpikeDelay, iNOSpikeDelay, fTau, fTauC, dRhoFixed, poisson_param, initial_random_seed);
	printf("Final rho values:\n");
	for (i = 0; i < no_synapses; i++){
		printf(" %f %f %f \n", syn[i].rho[siT], syn[i].rho[simulation_duration-1], syn[i].rho[simulation_duration-2]);
	}
    printf("%ld %ld\n", siT, simulation_duration);
	printf("\n");*/
	for (i = 0; i < no_synapses; i++){
		printf("syn(%d) t: %ld, c: %f, rho: %f\n", i, siT, syn[i].c[siT], syn[i].rho[siT]);
	}
	
	return 0;
}


void set_optimisation_sim_params(const gsl_vector * x){
    //double param_multiplier[8] = {1e-8, 1e-8, 1e-5, 1e-4, 1e-5, 1e-4, 1e2, 1e3}; //{1,1e-6,1,1,1,1,1000,1000};
	double param_multiplier[8] = {1,1,1,1,1,1,1,1};// {1e-8, 1e-10, 1e-6, 1e-6, 1e-6, 1e-7, 1e0, 1e1}; //{1,1e-6,1,1,1,1,1000,1000};
    //double param_multiplier[9] = {1e-8, 1e-10, 1e-6, 1e-6, 1e-6, 1e-7, 1e0, 1e1, 1e-8}; //{1,1e-6,1,1,1,1,1000,1000};
	//double param_multiplier[10] = {1e-10, 1e-10, 1e-6, 1e-6, 1e-6, 1e-7, 1e0, 1e1, 1e-8, 1e-10}; //{1,1e-6,1,1,1,1,1000,1000};
	double temp_reader;
    double delay_as_double;
    
    printf("DEBUG, setting new params. ");
    print_params();
    
	fTauC = (param_multiplier[0] * fTauCfixed + gsl_vector_get(x,0)) / param_multiplier[0]; //gsl_vector_get(x, 0);
	
	temp_reader = iCaSpikeDelay;
	delay_as_double = (param_multiplier[1] * iCaSpikeDelayFixed + gsl_vector_get(x,1)) / param_multiplier[1]; //gsl_vector_get(x,1);
	iCaSpikeDelay = (int) delay_as_double; //gsl_vector_get(x, 1);
	iNOSpikeDelay = iCaSpikeDelay; //gsl_vector_get(x, 1);
	printf("DEBUG C delay difference %g %g\n", (temp_reader - iCaSpikeDelay), delay_as_double);
	
	dCpre = (param_multiplier[2] * dCpreFixed + gsl_vector_get(x,2)) / param_multiplier[2]; //gsl_vector_get(x, 2);
	dCpost = (param_multiplier[3] * dCpostFixed + gsl_vector_get(x,3)) / param_multiplier[3]; //gsl_vector_get(x, 3);
	lfNMDARjump = (param_multiplier[4] * lfNMDARjumpFixed + gsl_vector_get(x,4)) / param_multiplier[4]; //gsl_vector_get(x, 4);
	
	dThetaD = (param_multiplier[5] * dThetaDfixed + gsl_vector_get(x,5)) / param_multiplier[5]; //gsl_vector_get(x, 5);

	temp_reader = dGammaD;
	dGammaD = (param_multiplier[6] * dGammaDfixed + gsl_vector_get(x,6)) / param_multiplier[6]; //gsl_vector_get(x, 6);
	dGammaP = (param_multiplier[7] * dGammaPfixed + gsl_vector_get(x,7)) / param_multiplier[7]; //gsl_vector_get(x, 7);
	
	//lfTauV = (param_multiplier[8] * lfTauVfixed + gsl_vector_get(x,8)) / param_multiplier[8]; //gsl_vector_get(x, 7);
	lfTauV = lfTauVfixed;
    lfTauNMDAR = fTauC;
    //lfTauNMDAR = (param_multiplier[9] * lfTauNMDARfixed + gsl_vector_get(x,9)) / param_multiplier[9]; //fTauC; //gsl_vector_get(x, 0);
    //dCdepol = (param_multiplier[8] * dCdepolFixed + gsl_vector_get(x,8)) / param_multiplier[8];
	
	printf("DEBUG gammaD difference %g\n", (temp_reader - dGammaD));
	
    printf("DEBUG, end of setting new params. ");
    print_params();
	/*V_MAX 1
	V_JUMP 1
	TAU_V 70.*/
}

void print_params(){
	printf("Parameters: tauC %0.30f, DC %d, DN, %d, Cpre %0.30f, Cpost %0.30f, Npre %0.30f, thetaD %0.30f, gammaD %0.30f, gammaP %0.30f, tau_v %0.30f, tauNO %0.30f Cdepol %0.30f\n", fTauC, iCaSpikeDelay, iNOSpikeDelay, dCpre, dCpost, lfNMDARjump, dThetaD, dGammaD, dGammaP, lfTauV, lfTauNMDAR, dCdepol);
}

void calculate_summary_data(Synapse *syn){
	int i, j;
	
	summary_outname = "output/summary_safo.dat";
	summary_outfile = fopen(summary_outname, "a");
	fprintf(summary_outfile, "\n\n\n\n\n%% LoopOffset, SynID, alpha_d, alpha_p, GammaD, GammaP, LTP zone, LTD zone, AmountLTP, AmountLTD\n");
	
	printf("Calculating summary data...\n");
	// Calculate alpha_d and alpha_p
	double alpha_d[no_synapses];
	double alpha_p[no_synapses];
	double above_NO_d[no_synapses];
	double above_NO_p[no_synapses];
	float theta_d = dThetaD;
	float theta_p = dThetaP;
	double ltp[no_synapses];
	double ltd[no_synapses];
	for(i = 0; i < no_synapses; i++){
		alpha_d[i] = 0;
		alpha_p[i] = 0;
		above_NO_d[i] = 0;
		above_NO_p[i] = 0;
		ltp[i] = 0;
		ltd[i] = 0;
	}
	int discard = 0;
	for(j = discard; j < (simulation_duration-1); j++){ //discard first 'discard' ms
		for(i = 0; i < no_synapses; i++){
			if ( syn[i].c[j] > theta_d){
				alpha_d[i]++;
			}
			if( syn[i].c[j] > theta_p){
				alpha_p[i]++;
			}
			if ( syn[i].NO_pre[j] > syn[i].no_threshold[j] ){
				if( syn[i].c[j] > theta_d ){
					above_NO_d[i]++;
				}
				else{ 
					above_NO_p[i]++;
				}
			}
			ltp[i] += syn[i].ltp[j];
			ltd[i] += syn[i].ltd[j];
		}
	}
	long t_total = simulation_duration - 3000;
	for(i = 0; i < no_synapses; i++){
		printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %f, LTD zone: %f, LTP: %lf, LTD: %lf\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
		alpha_d[i] /= t_total;
		alpha_p[i] /= t_total;
		/*above_NO_p[i] /= t_total;
		 above_NO_d[i] /= t_total;*/
		/*ltp[i] /= t_total;
		 ltd[i] /= t_total;*/
		printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %lf, LTD zone: %lf, LTP: %lf, LTD: %lf\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
		fprintf(summary_outfile, "%d, %d, %f, %f, %f, %f, %lf, %lf, %f, %f\n", loop_index, i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
	}
	
	printf("done.\n");
	fclose(summary_outfile);
}


Synapse* initialise_parameter_optimisation_sweep(int argc, char *argv[]){
	Synapse *syn;
	int i;
	siT = 0; // carry over global var from checkpointing method
	
	loadSimulationParameters(argc, argv);
	no_synapses = 17; // number of protocols we're trying
    dGammaDfixed = dGammaD;
    dGammaPfixed = dGammaP;
    
    dCpreFixed = dCpre;
    dCpostFixed = dCpost;
    lfNMDARjumpFixed = lfNMDARjump;
    dThetaDfixed = dThetaD;
    
    fTauCfixed = fTauC;
    lfTauNMDARfixed = lfTauNMDAR;
    
    iCaSpikeDelayFixed = iCaSpikeDelay;
    iNOSpikeDelayFixed = iNOSpikeDelay;
	
	lfTauVfixed = lfTauV;
    dCdepolFixed = dCdepol;
	
	for(i = 0; i < 17; i++){
		old_simulated_dw[i] = 0.0;
	}
	i = 0;
	
	// Memory allocation for array of synapses
	syn = (Synapse *) malloc( no_synapses * sizeof(Synapse));
	if (syn == NULL){
		perror("Memory allocation failure (syn array)\n");
		fprintf(logfile, "ERROR: Memory allocation failure (syn array)\n");
	}
	else{
		fprintf(logfile, "syn array successfully assigned %d synapses\n", no_synapses);
	}
	
	synapse_memory_init(syn);
	for (i = 0; i < no_synapses; i++){
		syn[i].rho[siT] = initial_rho;
		syn[i].c[siT] = initial_c; 
		//syn[i].V_pre[siT] = 0; //TODO: use variables to assign initial values to V_pre and NO_pre
		syn[i].NO_pre[siT] = 0;
	}
	
	fprintf(logfile, "Initialising spike times\n");
    // Safo and Regehr: 7 protocols, loop_index/dt gives offset from -300 in timesteps
	loop_index = 0;
	train30(syn[0].preT, syn[0].postT, simulation_duration);
	loop_index = 150 / dt;
	train30(syn[1].preT, syn[1].postT, simulation_duration);
	loop_index = 250 / dt;
	train30(syn[2].preT, syn[2].postT, simulation_duration);
	loop_index = 350 / dt;
	train30(syn[3].preT, syn[3].postT, simulation_duration);
	loop_index = 450 / dt;
	train30(syn[4].preT, syn[4].postT, simulation_duration);
	loop_index = 600 / dt;
	train30(syn[5].preT, syn[5].postT, simulation_duration);
	loop_index = 800 / dt;
	train30(syn[6].preT, syn[6].postT, simulation_duration);
	
	// 3xPF at 200Hz no change
	//trains_no_pf_stims = 3;
    // experimental data was for 2xPF, no change (slight depression)
    trains_no_pf_stims = 2;
	loop_index = 4. / dt;
	train32(syn[7].preT, syn[7].postT, simulation_duration);
	
	// 5xPF at 200Hz LTP
	trains_no_pf_stims = 5;
	loop_index = 4. / dt;
	train32(syn[8].preT, syn[8].postT, simulation_duration);
	// 5xPF at 33.3Hz LTP
	trains_no_pf_stims = 5;
	loop_index = 29. / dt;
	train32(syn[9].preT, syn[9].postT, simulation_duration);
	// 5xPF at 16.6Hz LTP
	trains_no_pf_stims = 5;
	loop_index = 59. / dt;
	train32(syn[10].preT, syn[10].postT, simulation_duration);
	trains_no_pf_stims = -1; // restore to defaults
	
	// Bidoret: 5 data points
	loop_index = 0;
    syn[11].uses_depol = 1;
	train34(syn[11].preT, syn[11].postT, simulation_duration);
	loop_index = 4. / dt;
    syn[12].uses_depol = 1;
	train34(syn[12].preT, syn[12].postT, simulation_duration);
	loop_index = 14. / dt;
    syn[13].uses_depol = 1;
	train34(syn[13].preT, syn[13].postT, simulation_duration);
	loop_index = 29. / dt;
    syn[14].uses_depol = 1;
	train34(syn[14].preT, syn[14].postT, simulation_duration);
	loop_index = 59. / dt;
    syn[15].uses_depol = 1;
	train34(syn[15].preT, syn[15].postT, simulation_duration);
	
	// Bidoret, single PF stim
	trains_no_pf_stims = 1;
    syn[16].uses_depol = 1;
	train34(syn[16].preT, syn[16].postT, simulation_duration);
	
	
    fprintf(logfile, "Spike times initialised\n");
	
	return syn;
}


#ifdef SIM_LOOP_PROGRAM
//TODO: main program loop
int main( int argc, char *argv[] ){
	int index_loop_counter;
	float loop_increment = 0.1; // ms
	
	summary_outname = "output/summary_safo.dat";
	summary_outfile = fopen(summary_outname, "a");
	fprintf(summary_outfile, "\n\n\n\n\n%% LoopOffset, SynID, alpha_d, alpha_p, GammaD, GammaP, LTP zone, LTD zone, AmountLTP, AmountLTD\n");
	
	// Don't forget, index_loop_counter is in units of DT
	//for(index_loop_counter = 0; index_loop_counter< SAFO_STEPS; index_loop_counter+=1000){
	for(index_loop_counter = 40; index_loop_counter< BIDORET_STEPS; index_loop_counter+=10000){
	//for(index_loop_counter = 0; index_loop_counter< PF_LOOP_STEPS; index_loop_counter+=100){
		printf("beginning loop %d\n", index_loop_counter);
		
		int i;
		long j, t;
		char outfile[FILE_NAME_LENGTH];

        // Initialise checkpointing
		/* Checkpoint_init:
		 1. Load parameters
		 2. openLogFile() (now in load params)
		 3. Reserve memory for array of synapses
		 4. Load Reset values
		 5. Reserve memory for members of Synapse(s)
		 */
		Synapse *syn;
		syn = checkpoint_init(argc, argv, syn);
		fflush(logfile);
        
        if (index_loop_counter > 0){
			loop_index = (index_loop_counter * loop_increment) / dt;
			//bidoret_index = (index_loop_counter * loop_increment) / dt;
		}
		else {
			loop_index = 0;
			//bidoret_index = 0;
		}
		//printf("safo index %d\n", safo_index);
		//printf("bidoret index %d\n", bidoret_index);
		printf("loop index %d\n", loop_index);
		
		//    for (k = 0; i < no_synapses; i++){
		//        fprintf(logfile, "DEBUG:: main: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
		//    }
		//    printf("DEBUG:: MEM TEST\n");
		//    printf("syn.ID: %d\n", syn[0].ID);
		fflush(stdout);



		// Load pre- and post- synaptic spike times into arrays for each synapse
		loadInitialSpikeTimes(syn);

		//CONSIDER: whether epsillon bound around thetaP=0 is necessary or not
		// Add Epsillon to thetaP (0) to ensure Ca switches off
		// if epsillon is 0.000001 then expect Ca to reach 0 from 1 in less than 14*tau
		/*if (dThetaP == 0){
		 dThetaP += 0.000001;
		 }
		 if (dThetaD == 0){
		 dThetaD += 0.000001;
		 }*/
	
		// Main simulation loop
		fprintf(logfile, "Entering main simulation loop\n");
		printf("Entering main simulation loop\n");
		//perform_parameter_optimisation_sim(syn);
		// Loop over discrete time steps up to simulation_duration
		for (t = siT; t < (simulation_duration-1); t++){
			//checkpoint_save(syn);
			// Update each synapse
			for (i = 0; i < no_synapses; i++){
				//printf("syn(%d) ", i);
				//updatePreSynapticVoltageTrace(&syn[i]);
				updatePreSynapticNOConcentration(&syn[i]);
				updateCalciumConcentration(&syn[i]);
				updateSynapticEfficacy(&syn[i]);
				//printf("t: %d, c: %f, rho: %f, NO: %f\n", siT, syn[i].c[siT-time_of_last_save], syn[i].rho[siT], syn[i].NO_pre[siT]);
			}
			//checkpoint_save(syn); // disable checkpointing
			siT++;
		}
		printf("DEBUG:: SIM OVER\n");
		checkpoint_save(syn);
		for (i = 0; i < no_synapses; i++){
			printf("syn(%d) t: %ld, c: %f, rho: %f\n", i, siT, syn[i].c[siT], syn[i].rho[siT]);
		}
		fprintf(logfile, "Simulation complete\n");
		printf("Simulation complete\n");

		// Debugging output after simulation has completed
		/*for (j = 0; j < (simulation_duration); j++){
			for (i = 0; i < no_synapses; i++){
				fprintf(logfile, "syn(%d).preT(%d): %u, postT(%d): %u, c: %f, rho: %f\n", i, j, syn[i].preT[j], j, syn[i].postT[j], syn[i].c[j], syn[i].rho[j]);
			}
		}
		fprintf(logfile, "siT: %d\n", siT);*/

		// Output to files loop
		if (!checkpointing){
			for (i = 0; i < no_synapses; i++){
				//sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
				sprintf(outfile, outfilepattern, syn[i].ID);
				printf("writing...%s\n", outfile);
				// not saving output file here
				//saveSynapseOutputFile(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iCaSpikeDelay, iNOSpikeDelay, fTau, fTauC, dRhoFixed, poisson_param, initial_random_seed);
			}
		}

        printf("Calculating summary data...\n");
		// Calculate alpha_d and alpha_p
		double alpha_d[no_synapses];
		double alpha_p[no_synapses];
		double above_NO_d[no_synapses];
		double above_NO_p[no_synapses];
		float theta_d = dThetaD;
		float theta_p = dThetaP;
		double ltp[no_synapses];
		double ltd[no_synapses];
		for(i = 0; i < no_synapses; i++){
			alpha_d[i] = 0;
			alpha_p[i] = 0;
			above_NO_d[i] = 0;
			above_NO_p[i] = 0;
			ltp[i] = 0;
			ltd[i] = 0;
		}
        int discard = 0;
		for(j = discard; j < (simulation_duration-1); j++){ //discard first 'discard' ms
			for(i = 0; i < no_synapses; i++){
				if ( syn[i].c[j] > theta_d){
					alpha_d[i]++;
				}
				if( syn[i].c[j] > theta_p){
					alpha_p[i]++;
				}
				if ( syn[i].NO_pre[j] > syn[i].no_threshold[j] ){
					if( syn[i].c[j] > theta_d ){
						above_NO_d[i]++;
					}
					else{ 
						above_NO_p[i]++;
					}
				}
				ltp[i] += syn[i].ltp[j];
				ltd[i] += syn[i].ltd[j];
			}
		}
		long t_total = simulation_duration - 3000;
		for(i = 0; i < no_synapses; i++){
			printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %f, LTD zone: %f, LTP: %lf, LTD: %lf\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
			alpha_d[i] /= t_total;
			alpha_p[i] /= t_total;
			/*above_NO_p[i] /= t_total;
			above_NO_d[i] /= t_total;*/
			/*ltp[i] /= t_total;
			ltd[i] /= t_total;*/
			printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %lf, LTD zone: %lf, LTP: %lf, LTD: %lf\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
			fprintf(summary_outfile, "%d, %d, %f, %f, %f, %f, %lf, %lf, %f, %f\n", loop_index, i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
		}

		fflush(summary_outfile);
		// Free memory and exit
		//return finalise(0, syn);
		finalise(0, syn);
		
		printf("ending loop %d\n", index_loop_counter);
		//end of safo loop
	}
	
	fclose(summary_outfile);
	return 0;
}
#endif /* SIM_LOOP_PROGRAM */

//TODO: underlying synaptic efficacy evolution code
// Calculate synaptic efficacy for next time step
void updateSynapticEfficacy(Synapse *syn){
    double rho, drho;//, minTheta, rand_no, noise;
	
    rho = (*syn).rho[siT];
	drho = 0;
	
	(*syn).no_threshold[siT] = fmax(( fThetaNO * ( 1 - ((*syn).c[siT] / fThetaNO2) ) ), 0);
	//CONSIDER: do I need an epsillon bound above NO_theshold when threshold=0?
	/*if((*syn).no_threshold[siT] == 0){
		(*syn).no_threshold[siT] += 0.000001;
	}*/
	
	//TODO: which version of the LTP rule do we wish to implement?
	//if ( h((*syn).c[siT], dThetaP) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
	if ( !(h((*syn).c[siT], dThetaD)) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
		// rho update is weight dependent here
		(*syn).ltp[siT] = ((dGammaP * (1 - rho)) / fTau) * dt; // had previously omitted divide by tau here
	}
	else{
		(*syn).ltp[siT] = 0;
	}
	
	if ( h((*syn).c[siT], dThetaD) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
		// rho update is weight dependent here
		(*syn).ltd[siT] = ((dGammaD * rho) / fTau) * dt; // had previously omitted divide by tau here
	}
	else{
		(*syn).ltd[siT] = 0;
	}
	
    //drho = (-rho * (1.0 - rho) * (dRhoFixed - rho)) + (dGammaP * (1 - rho) * h((*syn).c[siT], dThetaP) * h((*syn).NO_pre[siT], (fThetaNO*(1-((*syn).c[siT]/fThetaNO2))) )) - (dGammaD * rho * h((*syn).c[siT], dThetaD) * h((*syn).NO_pre[siT], (fThetaNO*(1-((*syn).c[siT]/fThetaNO2))) ));
	//drho = (-rho * (1.0 - rho) * (dRhoFixed - rho)) + (*syn).ltp[siT] - (*syn).ltd[siT];
	drho = (*syn).ltp[siT] - (*syn).ltd[siT]; // dt included in calculation of ltp and ltd
	
    // Add noise
    /*minTheta = fmin(dThetaP, dThetaD);
    if (h((*syn).c[siT], minTheta) > 0){ // Noise is on
        rand_no = (double) gasdev(&random_seed);
        noise = dSigma * sqrt(iTau) * rand_no;
        printf("\nNoise is active, rand_no: %f, noise: %f\n", rand_no, noise);
        drho += noise;
    }*/
	

	// drho /= fTau; // moved into calculation of LTP/D above (otherwise those vars omitted tau dependence)
    if ( (rho + drho) > 0 ){
        (*syn).rho[siT + 1] = rho + drho; // Euler forward method // dt included in calculation of ltp and ltd
    }
    else{
		// hard lower bound for numerical reasons
        (*syn).rho[siT + 1] = 0.0;
    }
}


// Simple Heaviside implemenation for comparing calcium
// concentration with a threshold value
/*BOOL h(Synapse *syn, double theta){
    if ( (*syn).c[siT] < theta)
        return 0;
    else
        return 1;
}*/
BOOL h(float c, double theta){
	//switched to strictly greater than to see if it helps with C_cs only situation
    if ( c >= theta )
        return 1;
    else
        return 0;
}

// Calculate synaptic calcium concentration for next time step
void updateCalciumConcentration(Synapse *syn){
    double c, dc;
    c = (*syn).c[siT];
	
	//TODO: change whether we fix Ca after a post spike or allow it to evolve dynamically here
	
	// Allow normal evolution of Ca from spikes
	/*dc = (-c / fTauC) + calciumFromPreSynapticSpikes(syn) + calciumFromPostSynapticSpikes(syn);
	 (*syn).c[siT + 1] = c + dc; // Euler forward method
	 */
	
	// Fix Ca at constant level on PC depolarisation
	/*if (calciumFromPostSynapticSpikes(syn) > 0){ 
		// post-synaptic depolarisation: fix Ca concentration at constant upper-level
		(*syn).c[siT + 1] = ((double) (*syn).postT[siT]) * dCpost;
	}
	else{
		dc = (-c / fTauC) + calciumFromPreSynapticSpikes(syn) + calciumFromPostSynapticSpikes(syn);
		(*syn).c[siT + 1] = c + dc; // Euler forward method
	}*/
	
	// PC depolarisation enforces lower bound on Ca, but PF can still modify Ca above this level
	if ((*syn).postT[siT - 1] == 1){
		// We've already applied the depolarisation dependent calcium influx on a previous timestep, now use normal
		// dynamics for PF dependent calcium and try to correct for leakage of depolarisation dependent calcium.
		dc = (-(c-dCpost) / fTauC);
		//dc = (-c / fTauC) + calciumFromPreSynapticSpikes(syn) + (dCpost / fTauC);
		(*syn).c[siT + 1] = c + (dc * dt) + calciumFromPreSynapticSpikes(syn); // Euler forward method
	}
	else{
		 dc = (-c / fTauC);
		 (*syn).c[siT + 1] = c + (dc * dt) + calciumFromPreSynapticSpikes(syn) + calciumFromPostSynapticSpikes(syn); // Euler forward method
	}
}


// Calculate contribution to next synaptic calcium concentration
// from pre-synaptic spikes
// Note: there is a delay iCaSpikeDelay before calcium from a
// pre-synaptic spike enters the synaptic cleft
double calciumFromPreSynapticSpikes(Synapse *syn){
    double d;

    //printf("preT: %u ", (*syn).preT[siT]);

    if (siT < iCaSpikeDelay){
        d = 0.0;
    }
    else if( (siT >= iCaSpikeDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iCaSpikeDelay]) * dCpre;
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in calciumFromPreSynapticSpikes()\n");
		printf("ERROR: unexpected situation in calciumFromPreSynapticSpikes()\n");
    }

    return d;
}


// Calculate contribution to next synaptic calcium concentration
// from post-synaptic spikes
double calciumFromPostSynapticSpikes(Synapse *syn){
    double d;
    //printf("postT: %u ", (*syn).postT[siT]);
    if ((*syn).uses_depol == 1){
        d = ((double) (*syn).postT[siT]) * dCdepol;
    }
    else{
        d = ((double) (*syn).postT[siT]) * dCpost;
    }
    return d;
}


// Update VGCC_avail variable, representing proportion of pre-synaptic
// VGCCs which are open
// Note: there is a delay iVGCCOpeningDelay before Magnesium block
// is released by a pre-synaptic spike
// This variable tracks the proportion of NMDARs which have Glutatmate attached (6/12/13) 
// ie. which are activated [ignore the function and variable names]
void updatePreSynapticVoltageTrace(Synapse *syn){
	//PF NMDARs
	double vgcc, dvgcc;
	
	vgcc = (*syn).V_pre[siT];
	dvgcc = (-vgcc / lfTauV);
	
	(*syn).V_pre[siT + 1] = vgcc + (dvgcc * dt) + voltageTraceFromPreSynapticSpikes(syn);
	
	if ((*syn).V_pre[siT + 1] > fVmax){
		(*syn).V_pre[siT + 1] = fVmax;
	}
}


float voltageTraceFromPreSynapticSpikes(Synapse *syn){
    float d;
	
    if (siT < iVOpeningDelay){
        d = 0.0;
    }
    else if( (siT >= iVOpeningDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iVOpeningDelay]) * lfVjump;
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in voltageTraceFromPreSynapticSpikes()");
    }
	
    return d;
}


// NMDAR state leads to NO concentration
// This variable tracks the concentration of NO released via
// convolution of NMDAR activation level with spikes
void updatePreSynapticNOConcentration(Synapse *syn){
	double no, dno;
	
    //printf("DEBUG NO_pre %lf \n", (*syn).NO_pre[siT]);
    
	no = (*syn).NO_pre[siT];
	dno = (-no / lfTauNMDAR);
	
	//if ((no + dno) < fNMDARmax){
	(*syn).NO_pre[siT + 1] = no + (dno * dt) + nmdarFromPreSynapticSpikes(syn);
	//}
	/*else{
		(*syn).NO_pre[siT + 1] = fNMDARmax;
	}*/
}


// A pre-synaptic spike combined with an already depolarised membrane leads to an influx of calcium, which in turn leads to a release of NO
double nmdarFromPreSynapticSpikes(Synapse *syn){
	double d;
	//CONSIDER: could potentially look-ahead to next VGCC_avail value ([siT+1]), if I process them in the right order
	
	if (siT < iNOSpikeDelay){
        d = 0.0;
    }
    else if( (siT >= iNOSpikeDelay) && ( siT < (simulation_duration - 1) ) ){
        // Simple model, no dependence on presynaptic potential unblocking of NMDAR
		d = ((double) (*syn).preT[siT - iNOSpikeDelay]) * lfNMDARjump;
		
		// Dependence on activation of presynaptic NMDAR
		// we need the equal delay on the spike and on the V state in order to have consistency
		//d = ((double) (*syn).preT[siT - iNOSpikeDelay]) * lfNMDARjump * (*syn).V_pre[siT - iNOSpikeDelay];
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in nmdarFromPreSynapticSpikes()");
    }
	
	// modified here to account for lack of V_pre
	//d = (*syn).V_pre[siT] * ((double)(*syn).preT[siT]) * fNMDARjump;
	//removed following to add delay on NO influx
	//d = ((double)(*syn).preT[siT]) * lfNMDARjump;
	
	return d;
}


// Setup spike times (hard-coded version)
// preT[i] = 1 means a spike occurs at time i
// preT[i] = 0 implies no spike at time i
void loadInitialSpikeTimes(Synapse *syn){
    int i;
    fprintf(logfile, "Initialising spike times\n");
//    fflush(logfile);
//    fprintf(logfile, "DEBUG:: syn(%d).preT[0] is %d\n", 0, (*syn).preT[0]);
//    fflush(logfile);
//    syn[0].preT[0] = 1;
//    syn[0].postT[0] = 0;
//    printf("DEBUG:: first spikes\n");
//    for (i = 1; i < simulation_duration; i++){
//        syn[0].preT[i] = 0;
//        syn[0].postT[i] = 0;
//    }
    for (i = 0; i < no_synapses; i++){
        (*train_fn)(syn[i].preT, syn[i].postT, simulation_duration);
    }
    fprintf(logfile, "Spike times initialised\n");
    //fflush(logfile);
}


void synapse_memory_init(Synapse *syn){
    int i;
    double * local_c;
    double * local_rho;
    unsigned int * local_preT;
    unsigned int * local_postT;
	double * local_v_pre;
	double * local_no_pre;
	double * local_ltp;
	double * local_ltd;
	double * local_no_threshold;
    //Synapse * local_synapse;
    fprintf(logfile, "Synapse simulator initialising.\n");

    for (i = 0; i < no_synapses; i++){
//        // Memory allocation for each synapse
//        local_synapse = (Synapse *) malloc( sizeof(Synapse) );
//        if (local_synapse == NULL){
//            perror("Memory allocation error (Synapse)\n");
//            fprintf(logfile, "ERROR: Memory allocation failure (Synapse)\n");
//        }
//        else{
//            (syn[i]) = local_synapse;
//            fprintf(logfile, "syn(%d) successfully assigned\n", i);
//        }
        // default to not using depolarisation;
        (syn[i]).uses_depol = 0;
        // Set synapse ID
        (syn[i]).ID = siID;
        siID++;
        fprintf(logfile, "Set synaptic id to: %d\n", (syn[i]).ID);

        // Memory allocation for c(t) array
        local_c = (double *) calloc( (simulation_duration), sizeof(double));
		//local_c = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_c == NULL){
            perror("Memory allocation failure (c)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (c)\n");
        }
        else{//removed (*syn) to allow for array based syn[0]
            (syn[i]).c = local_c;
            syn[i].c[0] = initial_c; // I shorten c[] to newly remaining length of sim if loading from a checkpoint!
            fprintf(logfile, "syn(%d).c successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
        }
        // Memory allocation for rho(t) array
        local_rho = (double *) calloc(simulation_duration, sizeof(double) );
        if (local_rho == NULL){
            perror("Memory allocation failure (rho)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (rho)\n");
        }
        else{
            (syn[i]).rho = local_rho;
            //syn[i].rho[0] = initial_rho;
            fprintf(logfile, "syn(%d).rho successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).rho(0): %lf\n", i, syn[i].rho[0]);
        }

        // Memory allocation for preT(t) array
        // CONSIDER: using calloc instead of malloc for spike time arrays (defaults to 0)
        local_preT = (unsigned int *) calloc( (simulation_duration), sizeof(unsigned int) );
        if (local_preT == NULL){
            perror("Memory allocation failure (preT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (preT)\n");
        }
        else{
            (syn[i]).preT = local_preT;
            fprintf(logfile, "syn(%d).preT successfully assigned\n", i);
            //(syn[i]).preT[0] = 99;  //
            //fprintf(logfile, "DEBUG:: syn(%d).preT[0] is %d\n", i, syn[i].preT[0]);
            //fflush(logfile);
        }
        // Memory allocation for postT(t) array
        local_postT = (unsigned int *) calloc( (simulation_duration), sizeof(unsigned int) );
        if (local_postT == NULL){
            perror("Memory allocation failure (postT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (postT)\n");
        }
        else{
            (syn[i]).postT = local_postT;
            fprintf(logfile, "syn(%d).postT successfully assigned\n", i);
        }
		
		// Memory allocation for VGCC and NO variables
		local_v_pre = (double *) calloc( (simulation_duration), sizeof(double) );
        if (local_v_pre == NULL){
            perror("Memory allocation failure (V_pre)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (V_pre)\n");
        }
        else{
            (syn[i]).V_pre = local_v_pre;
            fprintf(logfile, "syn(%d).V_pre successfully assigned\n", i);
        }
		local_no_pre = (double *) calloc( (simulation_duration), sizeof(double) );
        if (local_no_pre == NULL){
            perror("Memory allocation failure (NO_pre)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (NO_pre)\n");
        }
        else{
            (syn[i]).NO_pre = local_no_pre;
            fprintf(logfile, "syn(%d).NO_pre successfully assigned\n", i);
        }
		
		// Memory allocation for debugging LTP and LTD variables
		local_ltp = (double *) calloc( (simulation_duration), sizeof(double) );
        if (local_ltp == NULL){
            perror("Memory allocation failure (ltp)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (ltp)\n");
        }
        else{
            (syn[i]).ltp = local_ltp;
            fprintf(logfile, "syn(%d).ltp successfully assigned\n", i);
        }
		local_ltd = (double *) calloc( (simulation_duration), sizeof(double) );
        if (local_ltd == NULL){
            perror("Memory allocation failure (ltd)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (ltd)\n");
        }
        else{
            (syn[i]).ltd = local_ltd;
            fprintf(logfile, "syn(%d).ltd successfully assigned\n", i);
        }
		local_no_threshold = (double *) calloc( (simulation_duration), sizeof(double) );
        if (local_no_threshold == NULL){
            perror("Memory allocation failure (no_threshold)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (no_threshold)\n");
        }
        else{
            (syn[i]).no_threshold = local_no_threshold;
            fprintf(logfile, "syn(%d).no_threshold successfully assigned\n", i);
        }
    }
    fprintf(logfile, "Initialisation of simulator complete\n");
}


int finalise(int status, Synapse *syn){
    int i;
    if (status == 0){
        fprintf(logfile, "Synapse simulator exiting successfully\n");
        for (i = 0; i < no_synapses; i++){
            free((syn[i]).c);
            free((syn[i]).rho);
            free((syn[i]).preT);
            free((syn[i]).postT);
			free((syn[i]).V_pre);
			free((syn[i]).NO_pre);
			free((syn[i]).ltp);
			free((syn[i]).ltd);
			free((syn[i]).no_threshold);
        }
        free(syn);
        fprintf(logfile, "Memory freed\n");
        fprintf(logfile, "Exiting\n");
        closeLogFile(logfile);
        return 0;
    }
    else{
        fprintf(logfile, "An error occurred: exiting\n");
        closeLogFile(logfile);
		exit(1);
        return 1;
    }
}
