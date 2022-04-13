#include "NUCLEARDATA/NuclearData.cpp"
#include "NUCLEARDATA/H2.cpp"
#include "NUCLEARDATA/He3.cpp"

class Nucleus{
        
    
    public:
    
    // ATOMIC NUMBER //
    int A;
    
    // WOOD SAXON PARAMETERS //
    double NucleusRadius;
    double SurfaceDiffusiveness;
    
    // PARTICIPANT STATUS AND COLLISION NUMBER //
    int *ParticipantStatus;
    int *CollisionNumber;
    
    // NUMBER OF COLLISIONS //
    int NumberOfCollisions;
    int NumberOfParticipants;
    
    // SAMPLE A NUCLEON POSITION //
    void GlauberSamplePosition(double *x){
        
        int Accept=0;
        
        while(Accept==0){
            
            // SET RANDOM POSITIONS //
            x[0]=8.0*NucleusRadius*(RandomNumberGenerator::rng()-0.5);
            x[1]=8.0*NucleusRadius*(RandomNumberGenerator::rng()-0.5);
            x[2]=8.0*NucleusRadius*(RandomNumberGenerator::rng()-0.5);
            
            // GET ACCEPTANCE PROBABILITY //
            double xR=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
            double Probability=(1.0+exp(-NucleusRadius/SurfaceDiffusiveness))/(1.0+exp((xR-NucleusRadius)/SurfaceDiffusiveness));
            
            // CHECK ACCEPTANCE CRITERION //
            if(RandomNumberGenerator::rng()<Probability){
                Accept=1;
            }
            
        }
        
    }
    
    public:
    // NUCLEON POSITIONS //
    double **x; double xBar[3];
    
    // SAMPLE NUCLEON POSITIONS FROM WOOD SAXON DISTRIBUTION //
    void SetNucleonPositions(){
        
        // ALLOCATE POSITION ARRAY //
        x=new double*[A];
        
        for(int n=0;n<A;n++){
            x[n]=new double[3];
        }
        
        //////////////////////////////////////
        // SAMPLE POSITIONS OF EACH NUCLEON //
        //////////////////////////////////////
        
        // p //
        if(A==1){
            
            x[0][0]=0.0; x[0][1]=0.0; x[0][2]=0.0;
        }
        
        // d //
        else if(A==2){
            
            // SET DEUTERON DATASET //
            H2::Init();
            
            // GET NUCLEON POSITIONS //
            H2::GetNucleonPositions(x[0],x[1]);
            
        }
        
        // He3 //
        else if(A==3){
            
            // SET HELIUM THREE DATASET //
            He3::Init();
            
            // GET NUCLEON POSITIONS //
            He3::GetNucleonPositions(x[0],x[1],x[2]);
            
        }
        
        // GLAUBER MODEL FOR A>3 //
        else{
            
            for(int n=0;n<A;n++){
                
                // GLAUBER SAMPLING //
                GlauberSamplePosition(x[n]);
                
            }
        }
        
        // GET CENTER OF MASS //
        xBar[0]=0.0; xBar[1]=0.0; xBar[2]=0.0;
        
        for(int n=0;n<A;n++){
            xBar[0]+=x[n][0]/double(A);   xBar[1]+=x[n][1]/double(A);   xBar[2]+=x[n][2]/double(A);
        }
        
        // SHIFT CENTER OF MASS AND IMPACT PARAMETER //
        for(int n=0;n<A;n++){
            x[n][0]-=xBar[0];   x[n][1]-=xBar[1];   x[n][2]-=xBar[2];
            
        }
        
        xBar[0]=0.0; xBar[1]=0.0; xBar[2]=0.0;
        
    }

    
    // SINGLE NUCLEON THICKNESS //
    double ComputeNucleonThickness(double x,double y,double xCenter,double yCenter){
        
        double SqrDistance=(x-xCenter)*(x-xCenter)+(y-yCenter)*(y-yCenter);
        
        return exp(-0.5*SqrDistance/BG)/(2.0*M_PI*BG);
    }
    
    
    public:
    
    // CREATE OUTPUT OF NUCLEON POSITIONS //
    void OutputNucleonPositions(double b0X, double b0Y, std::string fname){
                
        // OPEN OUTPUT STREAM //
        std::ofstream OutStream;
        
        OutStream.open(fname);
        
        // CREATE OUTPUT //
        OutStream << "#NUCLEON POSITIONS" << std::endl;
        
        for(int n=0;n<A;n++){
            OutStream << x[n][0]+b0X << " " << x[n][1]+b0Y << " " << x[n][2] << " " << ParticipantStatus[n] << " " << CollisionNumber[n] << std::endl;
        }
        
        // CLOSE OUTPUT STREAM //
        OutStream.close();
    }
    

    // GET NUCLEON THICKNESS //
    double GetThickness(double x,double y){
        
        double TValue=0.0;
        
        // ADD CONTRIBUTION FROM ALL NUCLEONS //
        for(int n=0;n<A;n++){
            TValue+=ComputeNucleonThickness(x,y,this->x[n][0],this->x[n][1]);
        }
        
        return TValue;
    }
    
    // CONSTRUCTOR //
    Nucleus(int AtomicNumberArg){
        
        // SET NUCLEUS SIZE //
        A=AtomicNumberArg; 
        NucleusRadius=NuclearData::GetRadius(A);
        SurfaceDiffusiveness=NuclearData::GetSurface(A);
        
        // SET NUCLEON POSITIONS //
        SetNucleonPositions();
        
        ParticipantStatus=new int[A];
        CollisionNumber=new int[A];
        
        
    }
    
    
    // DESTRUCTOR //
    ~Nucleus(){
        
        for(int n=0;n<A;n++){
            delete[] x[n];
        }
        
        delete[] x;
        
        delete ParticipantStatus;
        delete CollisionNumber;
        
    }
    
};
