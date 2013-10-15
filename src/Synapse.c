#include "GeneralIncludes.h"
#include "Synapse.h"
#include "DataTools.h"


#include "NumericalTools.h"
#include "SpikeTrains.h"

int main( int argc, char *argv[] ){
    int i, j, t;
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

//    for (k = 0; i < no_synapses; i++){
//        fprintf(logfile, "DEBUG:: main: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
//    }
//    printf("DEBUG:: MEM TEST\n");
//    printf("syn.ID: %d\n", syn[0].ID);
    fflush(stdout);



    // Load pre- and post- synaptic spike times into arrays for each synapse
    loadInitialSpikeTimes(syn);

	//TODO: decide whether epsillon bound around thetaP=0 is necessary or not
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
    // Loop over discrete time steps up to simulation_duration
    for (t = siT; t < (simulation_duration-1); t++){
        //checkpoint_save(syn);
        // Update each synapse
        for (i = 0; i < no_synapses; i++){
            printf("syn(%d) ", i);
			//updatePreSynapticVoltageTrace(&syn[i]);
			updatePreSynapticNOConcentration(&syn[i]);
            updateCalciumConcentration(&syn[i]);
            updateSynapticEfficacy(&syn[i]);
            printf("t: %d, c: %f, rho: %f, NO: %f\n", siT, syn[i].c[siT-time_of_last_save], syn[i].rho[siT], syn[i].NO_pre[siT]);
        }
        checkpoint_save(syn);
        siT++;
    }
    printf("DEBUG:: SIM OVER\n");
    checkpoint_save(syn);
    for (i = 0; i < no_synapses; i++){
        printf("syn(%d) t: %d, c: %f, rho: %f\n", i, siT, syn[i].c[siT], syn[i].rho[siT]);
    }
    fprintf(logfile, "Simulation complete\n");
    printf("Simulation complete\n");

    // Debugging output after simulation has completed
    for (j = 0; j < (simulation_duration); j++){
        for (i = 0; i < no_synapses; i++){
            fprintf(logfile, "syn(%d).preT(%d): %u, postT(%d): %u, c: %f, rho: %f\n", i, j, syn[i].preT[j], j, syn[i].postT[j], syn[i].c[j], syn[i].rho[j]);
        }
    }
    fprintf(logfile, "siT: %d\n", siT);

    // Output to files loop
    if (!checkpointing){
        for (i = 0; i < no_synapses; i++){
            //sprintf(outfile, "output/01_syn_%.3d.dat", syn[i].ID);
            sprintf(outfile, outfilepattern, syn[i].ID);
            printf("writing...%s\n", outfile);
            saveSynapseOutputFile(outfile, &syn[i], siT, dCpre, dCpost, dThetaD, dThetaP, dGammaD, dGammaP, dSigma, iPreSpikeDelay, fTau, fTauC, dRhoFixed, poisson_param, initial_random_seed);
        }
    }

    // Calculate alpha_d and alpha_p
    float alpha_d[no_synapses];
    float alpha_p[no_synapses];
	float above_NO_d[no_synapses];
	float above_NO_p[no_synapses];
    float theta_d = dThetaD;
    float theta_p = dThetaP;
	float ltp[no_synapses];
	float ltd[no_synapses];
    for(i = 0; i < no_synapses; i++){
        alpha_d[i] = 0;
        alpha_p[i] = 0;
		above_NO_d[i] = 0;
		above_NO_p[i] = 0;
		ltp[i] = 0;
		ltd[i] = 0;
    }
    for(j = 3000; j < (simulation_duration-1); j++){ //discard first 3000ms
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
    int t_total = simulation_duration - 3000;
    for(i = 0; i < no_synapses; i++){
		printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %f, LTD zone: %f, LTP: %f, LTD: %f\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
		alpha_d[i] /= t_total;
        alpha_p[i] /= t_total;
		above_NO_p[i] /= t_total;
		above_NO_d[i] /= t_total;
		ltp[i] /= t_total;
		ltd[i] /= t_total;
        printf("Syn(%d), alpha_d: %f, alpha_p: %f, GammaD: %f, GammaP: %f, LTP zone: %f, LTD zone: %f, LTP: %f, LTD: %f\n", i, alpha_d[i], alpha_p[i], (alpha_d[i]*dGammaD), (alpha_p[i]*dGammaP), above_NO_p[i], above_NO_d[i], ltp[i], ltd[i]);
    }

    // Free memory and exit
    return finalise(0, syn);
}


// Calculate synaptic efficacy for next time step
void updateSynapticEfficacy(Synapse *syn){
    double rho, drho;//, minTheta, rand_no, noise;
	
    rho = (*syn).rho[siT];
	drho = 0;
	
	(*syn).no_threshold[siT] = fmax(( fThetaNO * ( 1 - ((*syn).c[siT] / fThetaNO2) ) ), 0);
	//TODO: do I need an epsillon bound above NO_theshold when threshold=0?
	/*if((*syn).no_threshold[siT] == 0){
		(*syn).no_threshold[siT] += 0.000001;
	}*/
	
	//TODO: which version of the LTP rule do we wish to implement?
	//if ( h((*syn).c[siT], dThetaP) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
	if ( !(h((*syn).c[siT], dThetaD)) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
		//TODO: remove bistability from rho update?
		(*syn).ltp[siT] = (dGammaP * (1 - rho)); 
	}
	else{
		(*syn).ltp[siT] = 0;
	}
	
	if ( h((*syn).c[siT], dThetaD) && h((*syn).NO_pre[siT], (*syn).no_threshold[siT] ) ){
		//TODO: remove bistability from rho update?
		(*syn).ltd[siT] = (dGammaD * rho);
	}
	else{
		(*syn).ltd[siT] = 0;
	}
	
    //drho = (-rho * (1.0 - rho) * (dRhoFixed - rho)) + (dGammaP * (1 - rho) * h((*syn).c[siT], dThetaP) * h((*syn).NO_pre[siT], (fThetaNO*(1-((*syn).c[siT]/fThetaNO2))) )) - (dGammaD * rho * h((*syn).c[siT], dThetaD) * h((*syn).NO_pre[siT], (fThetaNO*(1-((*syn).c[siT]/fThetaNO2))) ));
	//drho = (-rho * (1.0 - rho) * (dRhoFixed - rho)) + (*syn).ltp[siT] - (*syn).ltd[siT];
	drho = (*syn).ltp[siT] - (*syn).ltd[siT];
	
    // Add noise
    /*minTheta = fmin(dThetaP, dThetaD);
    if (h((*syn).c[siT], minTheta) > 0){ // Noise is on
        rand_no = (double) gasdev(&random_seed);
        noise = dSigma * sqrt(iTau) * rand_no;
        printf("\nNoise is active, rand_no: %f, noise: %f\n", rand_no, noise);
        drho += noise;
    }*/
	
	//TODO: should I introduce an explicit dt? (system wide change)

    drho /= fTau;
    if ( (rho + drho) > 0 ){
        (*syn).rho[siT + 1] = rho + drho; // Euler forward method
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
		dc = (-(c-dCpost) / fTauC) + calciumFromPreSynapticSpikes(syn);
		//dc = (-c / fTauC) + calciumFromPreSynapticSpikes(syn) + (dCpost / fTauC);
		(*syn).c[siT + 1] = c + dc; // Euler forward method
	}
	else{
		 dc = (-c / fTauC) + calciumFromPreSynapticSpikes(syn) + calciumFromPostSynapticSpikes(syn);
		 (*syn).c[siT + 1] = c + dc; // Euler forward method
	}
}


// Calculate contribution to next synaptic calcium concentration
// from pre-synaptic spikes
// Note: there is a delay iPreSpikeDelay before calcium from a
// pre-synaptic spike enters the synaptic cleft
double calciumFromPreSynapticSpikes(Synapse *syn){
    double d;

    printf("preT: %u ", (*syn).preT[siT]);

    if (siT < iPreSpikeDelay){
        d = 0.0;
    }
    else if( (siT >= iPreSpikeDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iPreSpikeDelay]) * dCpre;
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in calciumFromPreSynapticSpikes()");
    }

    return d;
}


// Calculate contribution to next synaptic calcium concentration
// from post-synaptic spikes
double calciumFromPostSynapticSpikes(Synapse *syn){
    double d;
    printf("postT: %u ", (*syn).postT[siT]);
    d = ((double) (*syn).postT[siT]) * dCpost;
    return d;
}


// Update VGCC_avail variable, representing proportion of pre-synaptic
// VGCCs which are open
// Note: there is a delay iVGCCOpeningDelay before Magnesium block
// is released by a pre-synaptic spike
/*void updatePreSynapticVoltageTrace(Synapse *syn){
	//PF NMDARs
	float vgcc, dvgcc;
	
	vgcc = (*syn).V_pre[siT];
	dvgcc = (-vgcc / (float)iTauV) + voltageTraceFromPreSynapticSpikes(syn);
	if ((vgcc + dvgcc) < fVmax){
		(*syn).V_pre[siT + 1] = vgcc + dvgcc;
	}
	else{
		(*syn).V_pre[siT + 1] = fVmax;
	}
}*/


/*float voltageTraceFromPreSynapticSpikes(Synapse *syn){
    float d;
	
    if (siT < iVOpeningDelay){
        d = 0.0;
    }
    else if( (siT >= iVOpeningDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iVOpeningDelay]) * fVjump;
    }
    else{ // This shouldn't happen!
        fprintf(logfile, "ERROR: unexpected situation in voltageTraceFromPreSynapticSpikes()");
    }
	
    return d;
}*/


// NMDAR state leads to NO concentration
void updatePreSynapticNOConcentration(Synapse *syn){
	double no, dno;
	
	no = (*syn).NO_pre[siT];
	dno = (-no / lfTauNMDAR) + nmdarFromPreSynapticSpikes(syn);
	
	//if ((no + dno) < fNMDARmax){
	(*syn).NO_pre[siT + 1] = no + dno;
	//}
	/*else{
		(*syn).NO_pre[siT + 1] = fNMDARmax;
	}*/
}


// A pre-synaptic spike combined with an already depolarised membrane leads to an influx of calcium, which in turn leads to a release of NO
double nmdarFromPreSynapticSpikes(Synapse *syn){
	double d;
	//CONSIDER: could potentially look-ahead to next VGCC_avail value ([siT+1]), if I process them in the right order
	
	if (siT < iPreSpikeDelay){
        d = 0.0;
    }
    else if( (siT >= iPreSpikeDelay) && ( siT < (simulation_duration - 1) ) ){
        d = ((double) (*syn).preT[siT - iPreSpikeDelay]) * lfNMDARjump;
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
	//float * local_v_pre;
	double * local_no_pre;
	float * local_ltp;
	float * local_ltd;
	float * local_no_threshold;
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
        // Set synapse ID
        (syn[i]).ID = siID;
        siID++;
        fprintf(logfile, "Set synaptic id to: %d\n", (syn[i]).ID);

        // Memory allocation for c(t) array
        local_c = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_c == NULL){
            perror("Memory allocation failure (c)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (c)\n");
        }
        else{//removed (*syn) to allow for array based syn[0]
            (syn[i]).c = local_c;
			// TODO: check if this hardcoded 0 is ok
            syn[i].c[0] = initial_c; // I shorten c[] to newly remaining length of sim if loading from a checkpoint!
            fprintf(logfile, "syn(%d).c successfully assigned\n", i);
            //fprintf(logfile, "DEBUG:: syn(%d).c(0): %lf\n", i, syn[i].c[0]);
        }
        // Memory allocation for rho(t) array
        local_rho = (double *) malloc( (simulation_duration) * sizeof(double) );
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
        local_preT = (unsigned int *) malloc( (simulation_duration) * sizeof(unsigned int) );
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
        local_postT = (unsigned int *) malloc( (simulation_duration) * sizeof(unsigned int) );
        if (local_postT == NULL){
            perror("Memory allocation failure (postT)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (postT)\n");
        }
        else{
            (syn[i]).postT = local_postT;
            fprintf(logfile, "syn(%d).postT successfully assigned\n", i);
        }
		
		// Memory allocation for VGCC and NO variables
		/*local_v_pre = (float *) malloc( (simulation_duration) * sizeof(float) );
        if (local_v_pre == NULL){
            perror("Memory allocation failure (V_pre)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (V_pre)\n");
        }
        else{
            (syn[i]).V_pre = local_v_pre;
            fprintf(logfile, "syn(%d).V_pre successfully assigned\n", i);
        }*/
		local_no_pre = (double *) malloc( (simulation_duration) * sizeof(double) );
        if (local_no_pre == NULL){
            perror("Memory allocation failure (NO_pre)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (NO_pre)\n");
        }
        else{
            (syn[i]).NO_pre = local_no_pre;
            fprintf(logfile, "syn(%d).NO_pre successfully assigned\n", i);
        }
		
		// Memory allocation for debugging LTP and LTD variables
		local_ltp = (float *) malloc( (simulation_duration) * sizeof(float) );
        if (local_ltp == NULL){
            perror("Memory allocation failure (ltp)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (ltp)\n");
        }
        else{
            (syn[i]).ltp = local_ltp;
            fprintf(logfile, "syn(%d).ltp successfully assigned\n", i);
        }
		local_ltd = (float *) malloc( (simulation_duration) * sizeof(float) );
        if (local_ltd == NULL){
            perror("Memory allocation failure (ltd)\n");
            fprintf(logfile, "ERROR: Memory allocation failure (ltd)\n");
        }
        else{
            (syn[i]).ltd = local_ltd;
            fprintf(logfile, "syn(%d).ltd successfully assigned\n", i);
        }
		local_no_threshold = (float *) malloc( (simulation_duration) * sizeof(float) );
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
			//free((syn[i]).V_pre);
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
        return 1;
    }
}
