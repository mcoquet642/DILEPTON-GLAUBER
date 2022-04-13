namespace  SaturationModel {
    
    // DEFAULT PARAMETERS //
    const double BG = 4.0*M_HBARC*M_HBARC;
    const double Sigma0 = 2.0*M_PI*BG;
    const double alphas = 0.3;
    
    const double Nc=3.0;
    double CA=Nc; double CF=(Nc*Nc-1.0)/(2.0*Nc);
    
    
    // VARIABLE PARAMETERS FOR HEAVY-ION FIT //
    double Qs0=0.63;
    double lambda=0.36;
    double NormalizationFactor;
    
    void SetNormalizationFactor(double dNchBydS,double CInfty,double eta_over_s,double nuEff){
        NormalizationFactor=3.16/((7.5*dNchBydS)*std::pow(CInfty/0.87,3.0/4.0)*std::pow(eta_over_s/0.16,1.0/3.0)*std::pow(nuEff/40.0,1.0/3.0));
    }
    
    // SET COLLIDING NUCLEI, CENTER OF MASS ENERGY AND RAPIDITY //
    int A1; int A2;
    double SqrtS_NN;
    double yRap;
    
    // SPACING //
    double aS;
    
    
#include "../MCGLAUBER/Nucleus.cpp"
#include "../MCGLAUBER/Participants.cpp"
    
    namespace NuclearParameters {
        
        double RA;
        double RB;
        
        Nucleus *NA;
        Nucleus *NB;
        
        void SampleNuclei() {
            
            // SAMPLE NUCLEI //
            NA = new Nucleus(A1);
            NB = new Nucleus(A2);
            
            // SET RADII //
            RA = NA->NucleusRadius;
            RB = NB->NucleusRadius;
            
        }
        
        
    }
    
    
    // DEFINE STRUCTURE FOR PASSING IMPACT PARAMETER //
    struct ImpactParameters{
        
        double bXA; double bYA;
        double bXB; double bYB;
        
    };
    
    
    struct DistributionValues{
        
        double TA; double TB;
        
    };
    
    
    
    namespace IPModel {
        
        double x0 = 0.000304;   double beta = 1.0;
        
        // Saturation parameters from https://arxiv.org/pdf/1711.11360.pdf
        double QsSqrAdjProton(double x){
            return CA/CF*std::pow(x0/x,lambda)*std::pow(1.0-x,beta)*Qs0*Qs0;
        }
        
        // GLUON DISTRIBUTION OF NUCLEUS A //
        double PhiA(double x, double k, double TA){
            return M_PI*(Nc*Nc-1.0)/(alphas*Nc)*(k*k)/(QsSqrAdjProton(x)*Sigma0*TA)*exp(-k*k/(QsSqrAdjProton(x)*Sigma0*TA));
        }
        
        // GLUON DISTRIBUTION OF NUCLEUS B //
        double PhiB(double x, double k, double TB){
            return M_PI*(Nc*Nc-1.0)/(alphas*Nc)*(k*k)/(QsSqrAdjProton(x)*Sigma0*TB)*exp(-k*k/(QsSqrAdjProton(x)*Sigma0*TB));
        }
        
        
        double QsAvg(double T) {
            return SqrtS_NN * x0 * std::pow( (CA/CF) * (Qs0/(SqrtS_NN*x0)) * (Qs0/(SqrtS_NN*x0)) * (Sigma0 * T) , 1.0/(2.0+lambda) );
        }
        
        
    }
    
    
    // CALCULATE ENERGY DENSITY IN UNITS GEV^3 -- BASED ON ANALYTICAL RESULT //
    double dEdyd2bAnalyticalResult (double TA, double TB) {
        
        using namespace IPModel;
        
        double QA = QsAvg(TA);
        double QB = QsAvg(TB);
        
        if (QA < 1e-30 || QB < 1e-30) { return 0.0; }
        else { return (Nc*Nc-1.0)/(16.0*M_PI*alphas*Nc*std::pow(M_PI, 1.0/2.0)) * (QA*QA)*(QB*QB) * (2.0*QA*QA*QA*QA+7.0*QA*QA*QB*QB+2.0*QB*QB*QB*QB) / (std::pow(QA*QA+QB*QB,5.0/2.0)); }
        
    }
    
    double Calculate_dEdyd2b_Analytical(double bXA,double bYA,double bXB,double bYB){
        
        using namespace NuclearParameters;
        
        double TA = NuclearParameters::NA->GetThickness(bXA,bYA);
        double TB = NuclearParameters::NB->GetThickness(bXB,bYB);
        
        return dEdyd2bAnalyticalResult(TA, TB);
        
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////

    
    void GenerateEvent(int EventID,int Ns,double *dEdyd2b){
        
        using namespace NuclearParameters;
        
        // Set RNG //
        RandomNumberGenerator::Init(EventID);
        
        //std::cerr << "#GENERATING EVENT " << EventID << std::endl;
        
        // CREATE NUCLEI //
        SampleNuclei();
        
        // SET GRID DIMENSIONS //
        double bXmin = -2.0*std::max(RA,RB);
        double bXmax = +2.0*std::max(RA,RB);

        double bYmin = -2.0*std::max(RA,RB);
        double bYmax = +2.0*std::max(RA,RB);
        
        // SET LATTICE SPACING //
        aS=4.0*std::max(RA,RB)/double(Ns-1);

        // SAMPLE IMPACT PARAMETER CIRCLE //
        double b0X=0.0; double b0Y=0.0;
//------------------------------------------------------------------------------------------
        int accepted = 0;
        while (accepted == 0) {

            b0X = 8.0*std::max(RA,RB)*(RandomNumberGenerator::rng()-0.5);
            b0Y = 8.0*std::max(RA,RB)*(RandomNumberGenerator::rng()-0.5);

            // CHECK THAT EVENT HAS ACCEPTABLE IMPACT PARAMETER AND NPart>0 //
            if (b0X*b0X+b0Y*b0Y < 16.0*std::max(RA,RB)*std::max(RA,RB)) {

                CheckParticipants(NA,NB,b0X,b0Y);
                
                int NPart=NA->NumberOfParticipants+NB->NumberOfParticipants;
                
                if(NPart>0){
                    accepted = 1;
                }
                
            }

        }
        
        // CALCULATE IMPACT PARAMETER //
        double b=sqrt(b0X*b0X+b0Y*b0Y);

        // SET IMPACT PARAMETER IN X DIRECTION //
        b0X=b; b0Y=0.0;
        
//------------------------------------------------------------------------------------------
        // SETUP GRID //
        double ETot=0.0;
        
//        #pragma omp parallel for collapse(2) reduction(+: ETot)
        for (int i = 0; i < Ns; i++) {
            for (int j = 0; j < Ns; j++) {
                
                double Current_dEdyd2b_fm3 = std::pow(NormalizationFactor,3.0/2.0)*Calculate_dEdyd2b_Analytical(bXmin + (bXmax-bXmin)/(Ns-1)*j+b0X/2, bYmin + (bYmax-bYmin)/(Ns-1)*i+b0Y/2, bXmin + (bXmax-bXmin)/(Ns-1)*j-b0X/2, bYmin + (bYmax-bYmin)/(Ns-1)*i-b0Y/2)/(M_HBARC*M_HBARC*M_HBARC);
                
                dEdyd2b[i+Ns*j]=Current_dEdyd2b_fm3;
                ETot+=Current_dEdyd2b_fm3;
                
            }
        }
        

//        std::cerr << "#dE/dy=" << (bXmax-bXmin)/(Ns-1)*(bYmax-bYmin)/(Ns-1)*ETot*M_HBARC << " [GeV]" << std::endl;
        
    }
    
    void SetupEventGenerator(int A1Val,int A2Val,double SqrtS_NNVal,double yRapVal){
        
        // SET COLLISION PARAMETERS //
        A1=A1Val; A2=A2Val; SqrtS_NN=SqrtS_NNVal; yRap=yRapVal;
        
    }

    
}
