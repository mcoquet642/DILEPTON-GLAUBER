#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace RandomNumberGenerator{
    
    //GSL RANDOM NUMBER GENERATORS AND SEED
    const gsl_rng_type *Type;
    gsl_rng **Generator;
    long int *MySEED;
    
    //INITIALIZATION OF RANDOM NUMBER GENERATOR
    void Init(long int BASE_SEED){
        
        // SETUP //
        gsl_rng_env_setup();
        
        Type=gsl_rng_default;
        
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        MySEED=new long int[NumberOfOpenMPThreads];
        Generator=new gsl_rng*[NumberOfOpenMPThreads];
        
        for(int tID=0;tID<NumberOfOpenMPThreads;tID++){
            
            Generator[tID]=gsl_rng_alloc(Type);
            
            long int SEED=BASE_SEED+tID;
            
            MySEED[tID]=SEED;
            
            gsl_rng_set(Generator[tID],SEED);
            
        }
        
    }
    
    // UNIFORM DISTRIBUTED RANDOM NUMBER in [0,1) //
    double rng(){
        
        int tID=omp_get_thread_num();
        return gsl_rng_uniform(Generator[tID]);
    }
    
}

