/*
 *  scratchpad.c
 *  PurkinjePlasticity
 *
 *  Created by David Higgins on 01/06/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "scratchpad.h"

//----------Purkinje trains------------

// train28 Lev-Ram'03 LTD protocol: PF+CF at 1Hz for 30 secs, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
//TODO: make PF into doublet
int train28(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength, no_stims;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 0; // measured in ms
    wavelength = 1000; //(int)(1.0 / rho);
	no_stims = 30;
	
    for (i = 0; i < no_stims; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
		//        for(j = 0; j < (wavelength-1); j++){
		//            preT[(i*wavelength) + j] = 0;
		//        }
        if((i*wavelength)+dt < simulation_duration){
            postT[(i*wavelength) + dt] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	//    for (;i < simulation_duration; i++){
	//        if( i > time_of_last_save){
	//            preT[i] = 0;
	//            postT[i] = 0;
	//        }
	//    }
	
    return 0;
}

// train29 Lev-Ram'03 LTP protocol: PF only at 1Hz for 300 stims, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, everything measured in ms
//TODO: make PF into doublet
int train29(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, wavelength, no_stims;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    //dt = 0; // measured in ms
    wavelength = 1000; //(int)(1.0 / rho);
	no_stims = 300;
	
    for (i = 0; i < no_stims; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
		//        for(j = 0; j < (wavelength-1); j++){
		//            preT[(i*wavelength) + j] = 0;
		//        }
        /*if((i*wavelength)+dt < simulation_duration){
		 postT[(i*wavelength) + dt] = 1;
		 }*/
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	//    for (;i < simulation_duration; i++){
	//        if( i > time_of_last_save){
	//            preT[i] = 0;
	//            postT[i] = 0;
	//        }
	//    }
	
    return 0;
}


// train30 Wang'00 protocol: X PF stims with fixed gap between them, followed by CF stim at offset,
// repeated at intervals Y times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
int train30(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, wavelength, recurrent_no_stims;
	//int dt;
	int no_pf_stims;
	int cf_offset; // ms
	int inter_pf_gap; // ms
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    //dt = 0; // measured in ms
    wavelength = 1000; //(int)(1.0 / rho);
	recurrent_no_stims = 30; // repeat protocol x times
	
	inter_pf_gap = 10; // ms
	no_pf_stims = 3;
	cf_offset = 62;
	
    for (i = 0; i < recurrent_no_stims; i++){
		for (int j = 0; j < no_pf_stims; j++){
			if (((i*wavelength) + (j*inter_pf_gap)) < simulation_duration){
				preT[(i*wavelength) + (j*inter_pf_gap)] = 1;
			}
		}
		
        if((i*wavelength)+cf_offset < simulation_duration){
            postT[(i*wavelength) + cf_offset] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	
    return 0;
}


// train31 Wang'00 protocol: X PF stims with fixed gap between them, no CF stim,
// repeated at intervals Y times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
int train31(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, wavelength, recurrent_no_stims;
	//int dt;
	int no_pf_stims;
	//int cf_offset; // ms
	int inter_pf_gap; // ms
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    //dt = 0; // measured in ms
    wavelength = 1000; //(int)(1.0 / rho);
	recurrent_no_stims = 30; // repeat protocol x times
	
	inter_pf_gap = 10; // ms
	no_pf_stims = 3;
	//cf_offset = 62;
	
    for (i = 0; i < recurrent_no_stims; i++){
		for (int j = 0; j < no_pf_stims; j++){
			if (((i*wavelength) + (j*inter_pf_gap)) < simulation_duration){
				preT[(i*wavelength) + (j*inter_pf_gap)] = 1;
			}
		}
		
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	
    return 0;
}