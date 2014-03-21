#ifndef SYNAPSE_H_
#define SYNAPSE_H_

#include "GeneralIncludes.h"

#include <gsl/gsl_multifit_nlin.h>

#define OPTIMISATION_PROGRAM

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
float fTauC;
double dt; // using 1 = 1ms, so 0.1=0.1ms

double dRhoFixed;

double dCpre;
double dCpost;
double dThetaD;
double dThetaP;
double dGammaD; //CONSIDER: does this really need to be double?
double dGammaP; //CONSIDER: does this really need to be double?

double dSigma;
int iCaSpikeDelay; // measured in timesteps // delay on Calcium increase
int iNOSpikeDelay; // measured in timesteps // delay on NO increase
double poisson_param;

int iVOpeningDelay;
double lfTauV;

double lfTauNMDAR;
double lfVjump;

double lfNMDARjump;

float fVmax;
float fNMDARmax;

float fThetaNO;
float fThetaNO2;

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
	float * no_threshold;
        int ID;
} Synapse;

struct fitting_data {
	Synapse * syn;
};

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

void set_optimisation_sim_params(const gsl_vector * x);
void calculate_summary_data(Synapse *syn);
int perform_parameter_optimisation_sim(Synapse *syn);
//float* cost_function(float *cost, Synapse *syn);
int cost_function(const gsl_vector * x, void * data, gsl_vector * f);
Synapse* initialise_parameter_optimisation_sweep(int argc, char *argv[]);

#endif /*SYNAPSE_H_*/
