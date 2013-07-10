#include "SpikeTrains.h"
#include "Synapse.h"
#include "NumericalTools.h" //for poidev()

// Setup spike times (hard-coded version)
// preT[i] = 1 means a spike occurs at time i
// preT[i] = 0 implies no spike at time i

// train1 is one pre- and no post- synaptic spikes
int train1(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    if (time_of_last_save < 0){
        preT[0] = 1;
        postT[0] = 0;
    }
    for (i = 1; i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = 0;
            postT[i] = 0;
        }
    }

    return 0;
}

// train2 is 20 pre-synaptic spikes and no post-synaptic spikes
int train2(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    for (i = 0; i < 20; i++){
        if( (i > time_of_last_save) && (i < simulation_duration) ){
            preT[i] = 1;
            postT[i] = 0;
        }
    }
    for (;i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = 0;
            postT[i] = 0;
        }
    }

    return 0;
}


// train3 is Poisson(n) distributed pre-synaptic spiking
int train3(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    int no_spikes = 0;

    for (i = 0;i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = (int) poidev(poisson_param, &random_seed);
            no_spikes += preT[i];
            postT[i] = 0;
        }
    }

    fprintf(logfile, "Number of poisson distributed pre-spikes: %d\n", no_spikes);
    return 0;
}


// train4 is Poisson(n) distributed post-synaptic spiking
int train4(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    int no_spikes = 0;

    for (i = 0;i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = 0;
            postT[i] = (int) poidev(poisson_param, &random_seed);
            no_spikes += postT[i];
        }
    }

    fprintf(logfile, "Number of poisson distributed post-spikes: %d\n", no_spikes);
    return 0;
}


// train5 is Poisson(n) distributed pre- and post-synaptic spiking
int train5(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    int no_pre_spikes = 0;
    int no_post_spikes = 0;

    for (i = 0;i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = (int) poidev(poisson_param, &random_seed);
            postT[i] = (int) poidev(poisson_param, &random_seed);
            no_pre_spikes += preT[i];
            no_post_spikes += postT[i];
        }
    }

    fprintf(logfile, "Number of poisson distributed pre-spikes: %d\n", no_pre_spikes);
    fprintf(logfile, "Number of poisson distributed post-spikes: %d\n", no_post_spikes);
    return 0;
}


// - Pasted from Pfister code

//TODO: Sjostrom style trains
// train6 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=0.1Hz, dt=+10ms
int train6(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 10000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
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


// train7 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=0.1Hz, dt=-10ms
int train7(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 10000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
//        for(j = 0; j < (wavelength-1); j++){
//            preT[(i*wavelength) + j] = 0;
//        }
        if((i*wavelength)+dt < simulation_duration){
            preT[(i*wavelength) + dt] = 1;
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


// train8 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=10Hz, dt=+10ms
int train8(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 100; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
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


// train9 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=10Hz, dt=-10ms
int train9(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 100; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
//        for(j = 0; j < (wavelength-1); j++){
//            preT[(i*wavelength) + j] = 0;
//        }
        if((i*wavelength)+dt < simulation_duration){
            preT[(i*wavelength) + dt] = 1;
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


// train10 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=20Hz, dt=+10ms
int train10(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 50; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
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


// train11 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=20Hz, dt=-10ms
int train11(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 50; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
//        for(j = 0; j < (wavelength-1); j++){
//            preT[(i*wavelength) + j] = 0;
//        }
        if((i*wavelength)+dt < simulation_duration){
            preT[(i*wavelength) + dt] = 1;
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


// train12 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=40Hz, dt=+10ms
int train12(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 25; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
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


// train13 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=40Hz, dt=-10ms
int train13(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 25; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
//        for(j = 0; j < (wavelength-1); j++){
//            preT[(i*wavelength) + j] = 0;
//        }
        if((i*wavelength)+dt < simulation_duration){
            preT[(i*wavelength) + dt] = 1;
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


// train14 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=50Hz, dt=+10ms
int train14(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 20; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
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


// train15 is Sjoestrom style frequency based stimulation for 60 repetitions
// Assumes preT and postT initialised to 0 using calloc()
// f=50Hz, dt=-10ms
int train15(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, dt, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    dt = 10; // measures in ms
    wavelength = 20; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
//        for(j = 0; j < (wavelength-1); j++){
//            preT[(i*wavelength) + j] = 0;
//        }
        if((i*wavelength)+dt < simulation_duration){
            preT[(i*wavelength) + dt] = 1;
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


// Wang 2005 triplet protocols
// train16 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=5, t2=5
int train16(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 5; // measures in ms
    t2 = 5;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            postT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            preT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train17 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=10, t2=10
int train17(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 10; // measures in ms
    t2 = 10;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            postT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            preT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train18 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=5, t2=15
int train18(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 5; // measures in ms
    t2 = 15;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            postT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            preT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train19 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=15, t2=5
int train19(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 15; // measures in ms
    t2 = 5;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            postT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            preT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// Wang 2005 triplet protocols Post-Pre-Post
// train20 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=5, t2=5
int train20(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 5; // measures in ms
    t2 = 5;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            preT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            postT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train21 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=10, t2=10
int train21(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 10; // measures in ms
    t2 = 10;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            preT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            postT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train22 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=5, t2=15
int train22(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 5; // measures in ms
    t2 = 15;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            preT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            postT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// train23 is Wang style triplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// t1=15, t2=5
int train23(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 15; // measures in ms
    t2 = 5;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
        if((i*wavelength)+t1 < simulation_duration){
            preT[(i*wavelength)+t1] = 1;
        }
        if((i*wavelength)+(t1+t2) < simulation_duration){
            postT[(i*wavelength) + (t1+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    return 0;
}


// Wang 2005 - Quadruplets
// train24 is Wang style quadruplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// T = 25
// Post-Pre-Pre-Post
int train24(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 6; // measures in ms
    t2 = 6;
    static int T = 7;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            postT[(i*wavelength)] = 1;
        }
        if(((i*wavelength)+t1) < simulation_duration){
            preT[(i*wavelength)+t1] = 1;
        }
        if(((i*wavelength)+(T)) < simulation_duration){
            preT[(i*wavelength) + (T)] = 1;
        }
        if(((i*wavelength)+(T+t2)) < simulation_duration){
            postT[(i*wavelength) + (T+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    T += 5;

    return 0;
}


// train25 is Wang style quadruplet based stimulation for 60 repetitions @ 1Hz
// Assumes preT and postT initialised to 0 using calloc()
// T = 25
// Pre-Post-Post-Pre
int train25(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i, t1, t2, wavelength;
    //float freq, rho;
    //freq = 1;
    //rho = 1000; // freq / 1000; // Convert Hz frequency to per millisecond
    t1 = 6; // measures in ms
    t2 = 6;
    static int T = 7;
    wavelength = 1000; //(int)(1.0 / rho);

    for (i = 0; i < 60; i++){
        //if( (i > time_of_last_save) && (i < simulation_duration) ){
        if ((i*wavelength) < simulation_duration){
            preT[(i*wavelength)] = 1;
        }
        if(((i*wavelength)+t1) < simulation_duration){
            postT[(i*wavelength)+t1] = 1;
        }
        if(((i*wavelength)+(T)) < simulation_duration){
            postT[(i*wavelength) + (T)] = 1;
        }
        if(  ((i*wavelength)+(T+t2)) < simulation_duration){
            preT[(i*wavelength) + (T+t2)] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }

    T += 5;

    return 0;
}


//TODO: reminder train26 is different in Pfister code compared to all other codebases
// Dual spike shot noise simulation
// Pre-spike occurs as poisson process, post-spike occurs T timesteps after pre-spike
int train26(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    int no_spikes = 0;
    int T = 23; // 693; // Delay between pre-spike and subsequent post-spike
    printf("DEBUG: here\n");

    for (i = 0;i < simulation_duration; i++){
        if( i > time_of_last_save){
            preT[i] = (int) poidev(poisson_param, &random_seed);
            no_spikes += preT[i];
            if(i < T){
                postT[i] = 0;
            }
            else{
                postT[i] = preT[i-T];
            }
        }
    }

    fprintf(logfile, "Number of poisson distributed pre-spikes: %d\n", no_spikes);
    return 0;
}


// Two shot noise processes, one with a single spike and one with a split 'dual' shape
// Single: pre-spike occurs as a poisson process, rate n1
// Dual: Pre-spike occurs as poisson process, post-spike occurs T timesteps after pre-spike
int train27(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
    int i;
    int no_simple_spikes = 0;
    int no_dual_spikes = 0;
    int T = 693; // Delay between pre-spike and subsequent post-spike

    for (i = 0;i < simulation_duration; i++){
        if( i > time_of_last_save){
            // Simple pre-spike
            preT[i] = (int) poidev(poisson_param, &random_seed);
            no_simple_spikes += preT[i];

            // Dual pre- post- spike
            postT[i] = (int) poidev(poisson_param, &random_seed);
            no_dual_spikes += postT[i];
            if (i > T){
                preT[i - T] += postT[i];
            }
        }
    }

    fprintf(logfile, "Number of poisson distributed simple pre-spikes: %d\n", no_simple_spikes);
    fprintf(logfile, "Number of poisson distributed dual pre- post-spikes: %d\n", no_dual_spikes);
    return 0;
}


//----------Purkinje trains------------

// train28 Lev-Ram'03 LTD protocol: PF+CF at 1Hz for 30 secs, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
//Modified: to assume that the PF stim is actually a paired-pulse
int train28(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	no_pf_stims = 4;
	cf_offset = 0;
	
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

// train29 Lev-Ram'03 LTP protocol: PF only at 1Hz for 300 stims, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, everything measured in ms
//Modified: to assume that the PF stim is actually a paired-pulse
int train29(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	recurrent_no_stims = 300; // repeat protocol x times
	
	inter_pf_gap = 10; // ms
	no_pf_stims = 4;
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


// train30 Wang'00 LTD protocol: X PF stims with fixed gap between them, followed by CF stim at offset,
// repeated at intervals 30 times, then no further inputs.
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
	
	// Wang LTD values
	inter_pf_gap = 10; //5; //10; // ms W:10
	no_pf_stims = 3; //2; //3; W:3
	cf_offset = 62; //normal wang ltd protocol
	// Generic CF before PF
	//inter_pf_gap = 5; //5; //10; // ms 
	//no_pf_stims = 2; //2; //3;
	//cf_offset = -30;
	// Safo LTD values
	inter_pf_gap = 10;
	no_pf_stims = 7;
	cf_offset = 500;
	cf_offset += 30; // 30 accounts for middle (4th) spike
	// Mathy values
	inter_pf_gap = 10;
	no_pf_stims = 5;
	cf_offset = 40;
	recurrent_no_stims = 25;
	wavelength = 2000;
	
    for (i = 0; i < recurrent_no_stims; i++){
		for (int j = 0; j < no_pf_stims; j++){
			if (((i*wavelength) + (j*inter_pf_gap)) < simulation_duration){
				preT[(i*wavelength) + (j*inter_pf_gap)] = 1;
			}
		}

        if(((i*wavelength)+cf_offset < simulation_duration) && ((i*wavelength)+cf_offset > -1)){
            postT[(i*wavelength) + cf_offset] = 1;
        }
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	
    return 0;
}


// train31 Wang'00 (pseudo-)LTP protocol: X PF stims with fixed gap between them, no CF stim,
// repeated at intervals 30 times, then no further inputs.
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


// train32 Mariano's LTP protocol: X PF stims with fixed gap between them, no CF stim,
// repeated at intervals 300 times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
int train32(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	recurrent_no_stims = 300; // repeat protocol x times
	
	inter_pf_gap = 5; // ms
	no_pf_stims = 5;
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


// train33 (LTD) Mariano's LTP protocol but with added CF stim: X PF stims with fixed gap between them, CF stim coincident with 5th PF stim,
// repeated at intervals 300 times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
int train33(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	recurrent_no_stims = 300; // repeat protocol x times
	
	inter_pf_gap = 5; // ms
	no_pf_stims = 5;
	cf_offset = 40;
	
	// try new parameters
	no_pf_stims = 2;
	cf_offset = 0;
	
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


// train34 Bidoret'09 LTD protocol: PC depolarisation of 120ms, 2PF stims one falling 60ms after start of PC depolarisation,
// the other Xms before, protocol repeated at intervals Y times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, everything measured in ms
int train34(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	recurrent_no_stims = 120; // repeat protocol x times
	
	inter_pf_gap = 30; //42;//5; // ms
	no_pf_stims = 2;
	//cf_offset = 40;
	
    for (i = 0; i < recurrent_no_stims; i++){
		for ( int j = 0; j < 120; j++){ // PC depolarisation
			if ((i*wavelength + j) < simulation_duration){
				postT[(i*wavelength) + j] = 1;
			}
		}
		for (int j = 0; j < no_pf_stims; j++){
			if (((i*wavelength) + 60 - (j*inter_pf_gap)) < simulation_duration){
				preT[(i*wavelength) + 60 - (j*inter_pf_gap)] = 1;
			}
		}
		
        //printf("DEBUG: i: %d, j: %d\n", (i*wavelength), ((i*wavelength)+dt));
        //}
    }
	
    return 0;
}


// train35 Predicted LTD protocol protocol: X PF stims with fixed gap between them, no CF stim,
// repeated at intervals 300 times, then no further inputs.
// Assumes preT and postT initialised to 0 using calloc()
// wavelength is gap between stims, dt is offset between pre and post stims, everything measured in ms
int train35(unsigned int * preT, unsigned int * postT, unsigned int simulation_duration){
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
	recurrent_no_stims = 300; // repeat protocol x times
	
	inter_pf_gap = 10; // ms
	no_pf_stims = 2;
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