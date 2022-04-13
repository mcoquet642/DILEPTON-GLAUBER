// DEFINE RANDOM NUMBER GENERATOR //
namespace RandomNumberGenerator {
    
    //INITIALIZATION OF RANDOM NUMBER GENERATOR
    void Init(long int BASE_SEED){
        srand48(BASE_SEED);
    }
    
    double rng(){
        return drand48();
    }
    
}


