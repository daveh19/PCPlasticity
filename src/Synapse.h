#ifndef SYNAPSE_H_
#define SYNAPSE_H_

#include "GeneralIncludes.h"

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
int simulation_duration;
int no_synapses;

double initial_c;
double initial_rho;

long initial_random_seed;
long random_seed;

int iTau;
int iTauC;

double dRhoFixed;

double dCpre;
double dCpost;
double dThetaD;
double dThetaP;
double dGammaD; //CONSIDER: does this really need to be double?
double dGammaP; //CONSIDER: does this really need to be double?

double dSigma;
int iPreSpikeDelay;
double poisson_param;

int iVOpeningDelay;
int iTauV;
int iTauNMDAR;
float fVjump;
float fNMDARjump;
float fVmax;
float fNMDARmax;
float fThetaNO;
float fThetaNO2;

int siT; // no longer static
int siID; // no longer static
int time_of_last_save;
int resume_offset;

typedef int BOOL;

BOOL checkpointing;
//BOOL a_restart;

FILE* logfile;
char* logfilename;
char* outfilepattern;
char logfilearray[FILE_NAME_LENGTH];
char outfilearray[FILE_NAME_LENGTH];

int (*train_fn)(unsigned int *, unsigned int *, unsigned int);

typedef struct Synapse{
        double * rho;
        double * c;
        unsigned int * preT;
        unsigned int * postT;
		float * V_pre;
		float * NO_pre;
        int ID;
} Synapse;


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
float nmdarFromPreSynapticSpikes(Synapse *syn);
void updatePreSynapticNOConcentration(Synapse *syn);

#endif /*SYNAPSE_H_*/
