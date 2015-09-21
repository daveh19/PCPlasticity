#ifndef SYNAPSE_H_
#define SYNAPSE_H_

#include "GeneralIncludes.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit_nlin.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

//#define SIM_LOOP_PROGRAM
//#define BASIC_ALL_EXPERIMENTS_SWEEP
//#define LM_OPTIMISATION_PROGRAM
//#define PR_MINIMISATION_PROGRAM
#define NM_MINIMISATION_PROGRAM

#ifdef PR_MINIMISATION_PROGRAM
	#define MINIMISATION_PROGRAM
#endif
#ifdef NM_MINIMISATION_PROGRAM
	#define MINIMISATION_PROGRAM
#endif
#ifdef BASIC_ALL_EXPERIMENTS_SWEEP
	#define MINIMISATION_PROGRAM
#endif

#define SAFO_STEPS (8001) /*(8001)*/
#define BIDORET_STEPS (800)
#define PF_LOOP_STEPS (601)

/*
    In order to enable save and resume from a checkpoint export the
    following variables from the command shell:
        SYNAPSE_CHECKFILE=filename to save checkpoint information to
        SYNAPSE_RESTART=y if this is a resume, unset it otherwise
*/

/* 
	Reminder: the way the timings are coded, viewed from time t (after an update)
     a delay D which was applied to the variables before the update is actually seen
     as a delay of D+1	 (so minimal delay right now is 1ms)
*/
long int simulation_duration; // measured in timesteps
int no_synapses;

double initial_c;
double initial_rho;

long initial_random_seed;
long random_seed;

//int iTau;
//int iTauC;
float fTau;
double fTauC;
double dt; // using 1 = 1ms, so 0.1=0.1ms

double dRhoFixed;

double dCpre;
double dCpost;
double dCdepol;
double dThetaD;
double dThetaP;
double dGammaD;
double dGammaP;

double dSigma;
int iCaSpikeDelay; // measured in timesteps // delay on Calcium increase
int iNOSpikeDelay; // measured in timesteps // delay on NO increase
double poisson_param;

int iVOpeningDelay;
double lfTauV;

double lfTauNMDAR;
double lfVjump;

double lfNMDARjump;

double fVmax;
double fNMDARmax;

double fThetaNO;
double fThetaNO2;

long int siT; // no longer static
int siID; // no longer static
long int time_of_last_save;
long int resume_offset;

typedef int BOOL;

BOOL checkpointing;
//BOOL a_restart;

FILE* summary_outfile;
char* summary_outname;
char summary_outarray[FILE_NAME_LENGTH];

FILE* logfile;
char* logfilename;
char* outfilepattern;
char logfilearray[FILE_NAME_LENGTH];
char outfilearray[FILE_NAME_LENGTH];

//int (*train_fn)(unsigned int *, unsigned int *, unsigned int);
int (*train_fn)(unsigned int *, unsigned int *, unsigned long);

typedef struct Synapse{
        double * rho;
        double * c;
        unsigned int * preT;
        unsigned int * postT;
		double * V_pre;
		double * NO_pre;
	double * ltp;
	double * ltd;
	double * no_threshold;
        int ID;
    int uses_depol;
} Synapse;

struct fitting_data {
	Synapse * syn;
};

double old_simulated_dw[17];// = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

int times_through_cost_function;
int times_through_cost_function_jacobian;

double dGammaDfixed;
double dGammaPfixed;

double dCpreFixed;
double dCpostFixed;
double dCdepolFixed;
double lfNMDARjumpFixed;
double dThetaDfixed;

double fTauCfixed;
double lfTauNMDARfixed;

double iCaSpikeDelayFixed;
int iNOSpikeDelayFixed; 

double lfTauVfixed;
double lfVjumpFixed;

void synapse_memory_init(Synapse *);
int finalise(int, Synapse *);
void loadInitialSpikeTimes(Synapse *);
double calciumFromPreSynapticSpikes(Synapse *);
double calciumFromPostSynapticSpikes(Synapse *);
void updateCalciumConcentration(Synapse *);
//BOOL h(Synapse *, double);
BOOL h(float, double);
void updateSynapticEfficacy(Synapse *);
float voltageTraceFromPreSynapticSpikes(Synapse *syn);
void updatePreSynapticVoltageTrace(Synapse *syn);
double nmdarFromPreSynapticSpikes(Synapse *syn);
void updatePreSynapticNOConcentration(Synapse *syn);

void print_params();
Synapse* initialise_parameter_optimisation_sweep(int argc, char *argv[]);
void set_optimisation_sim_params(const gsl_vector * x);
void calculate_summary_data(Synapse *syn);
int perform_parameter_optimisation_sim(Synapse *syn);
#ifdef LM_OPTIMISATION_PROGRAM
	int calculate_jacobian(const gsl_vector * x, void * data, gsl_matrix * J);
	int cost_function(const gsl_vector * x, void * data, gsl_vector * f);
#endif /* LM_OPTIMISATION_PROGRAM */
#ifdef PR_MINIMISATION_PROGRAM
	void pr_fdf(const gsl_vector *x, void * params, double * f, gsl_vector *df);
	void calculate_gradient(const gsl_vector * x, void * data, gsl_vector * J);
#endif /* PR_MINIMISATION_PROGRAM */
#ifdef MINIMISATION_PROGRAM
	double cost_function(const gsl_vector * x, void * data);
#endif /* MINIMISATION_PROGRAM */
//float* cost_function(float *cost, Synapse *syn);

#endif /*SYNAPSE_H_*/
