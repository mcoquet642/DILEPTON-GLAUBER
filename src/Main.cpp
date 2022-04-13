// STANDARD INCLUSIONS //
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sstream>
#include <iterator>
#include <cmath>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <cstring>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


// DEFINE TO ENABLE OF DISABLE MPI //
#ifndef ENABLE_MPI
#define ENABLE_MPI 0
#endif

// MPI //
#if(ENABLE_MPI==1)
#include "mpi.h"

namespace MPIBasic {
    int ID; int NumberOfNodes;
}
#endif

// EIGEN LIBRARY //
#include "Eigen/Dense"

// COMMANDLINE ARGUMENTS //
#include "IO/cfile.c"

// RANDOM NUMBER GENERATOR //
#include "MISC/GSLRandomGenerator.cpp"

// NON-EQUILIBRIUM CONSTANT //
const double CInfty = 0.87;

// EOS -- e(T)=Pi^2/30 nuEff(T) T^4  //
const double nuEff = 32.0;
//const double nuEff = 40.0;

// NUMBER OF CHARGED HADRONS PER UNIT ENTROPY  -- see https://arxiv.org/pdf/1908.02792.pdf //
const double dNchBydS = 1.0/6.7;
//const double dNchBydS = 1.0/7.5;

const double M_HBARC = 0.197;

const double alphaEM = 1.0/137.0; const double mllSqr = 0.0; const double qFSqrSum = 1.0/9.0+4.0/9.0+1.0/9.0;


// SATURATION MODEL BASED ON SMALL X TMDs //
#include "TMD/EventGenerator.cpp"

// MULTIPLICITY ESTIMATOR //
#include "MultiplicityEstimator.cpp"

// CENTRALITY SELECTION //
#include "CentralitySelection.cpp"

// BACKGROUND EVOLUTION FOR DILEPTON CALCULATION //
#include "HydroAttractor.cpp"
#include "PhaseSpaceDistribution.cpp"

// DILEPTON CALCULATION //
#include "DileptonRates.cpp"
#include "DileptonEstimator.cpp"

int main(int argc, char **argv) {
    
    ///////////////////////////////
    //INITIALIZE MPI ENVIRONMENT //
    ///////////////////////////////
    
    #if(ENABLE_MPI==1)
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //SET MPI ID's AND NUMBER OF NODES
    MPI_Comm_rank(MPI_COMM_WORLD,&MPIBasic::ID);
    MPI_Comm_size(MPI_COMM_WORLD,&MPIBasic::NumberOfNodes);
    
    if(MPIBasic::ID==0){
        std::cerr << "#SIMULATING WITH " << MPIBasic::NumberOfNodes << " MPI NODES" << std::endl;
    }
    #endif
    
    // SET COMMANDLINE ARGUMENTS //
    Konfig CommandlineArguments(argc,argv);
    
    // SET COLLISION SYSTEM AND KINETMATICS //
    int A1=208; int A2=208; double SqrtSNN=2760; double yRap=0.0; double EtaOverS=0.16;
    
    CommandlineArguments.Getval("A1",A1);
    CommandlineArguments.Getval("A2",A2);
    CommandlineArguments.Getval("s",SqrtSNN);
    CommandlineArguments.Getval("y",yRap);
    CommandlineArguments.Getval("etas",EtaOverS);
    
    // SETUP SATURATION MODEL //
    SaturationModel::SetupEventGenerator(A1,A2,SqrtSNN,yRap);
    SaturationModel::SetNormalizationFactor(dNchBydS,CInfty,EtaOverS,nuEff);
    
    // SET EVENT SAMPLE SIZE //
    int Ns=64;
    CommandlineArguments.Getval("Ns",Ns);
    
    // SET EVENT SAMPLE SIZE //
    int NumberOfEvents=262144;
    CommandlineArguments.Getval("NEv",NumberOfEvents);
    
    int doCentralityClass=0;
    CommandlineArguments.Getval("doCent",doCentralityClass);
    
    //////////////////////////
    // CENTRALITY SELECTION //
    //////////////////////////
    
    // SET CENTRALITY CLASSES //
    char CentralityInput[256]="0 5 10 20 30 40 50 60 70 80 90 100";
    CommandlineArguments.Getval("c",CentralityInput);
    
    std::string CentralityString(CentralityInput);
    std::stringstream ss(CentralityInput);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> CentralityStrings(begin, end);
    
    if (doCentralityClass){
        for(size_t i=0;i<CentralityStrings.size();i++){
            CentralitySelection::CentralityClasses.push_back(atof(CentralityStrings.at(i).c_str()));
        }
        
        char CentralityOutput[256]="CentralityClasses.txt";
        CommandlineArguments.Getval("co",CentralityOutput);
        
        
        // CREATE CENTRALITY CLASSES //
        CentralitySelection::CreateCentralityClasses(NumberOfEvents,Ns,EtaOverS);
        
        // CREATE OUTPUT //
        CentralitySelection::Output(CentralityOutput);
    }
    
    
    //////////////////////////
    // DILEPTON PRODUCTION  //
    //////////////////////////
    
    // SET ENERGY DENSITY ARRAY //
    double dEdyd2b[Ns*Ns];
    
    // SET EVENT ID //
    int EventID=1234;
    CommandlineArguments.Getval("ID",EventID);
    
    // GENERATE EVENT //
    SaturationModel::GenerateEvent(EventID,Ns,dEdyd2b);
    
    double dNchdEta= MultiplicityEstimator::GetdNchdEta(dEdyd2b,Ns,SaturationModel::aS,EtaOverS);
    
    if (doCentralityClass){
    	int iCentralityClass=CentralitySelection::GetCentralityClass(dNchdEta);
    	std::cout << "#EVENTID="<< EventID << " MULTIPLICITY dNch/dEta=" << dNchdEta << " CENTRALITY CLASS=" << iCentralityClass << " (" << CentralitySelection::CentralityClasses.at(iCentralityClass) << "-" << CentralitySelection::CentralityClasses.at(iCentralityClass+1) << "%)"  << endl;
    }else{
    	std::cout << "#EVENTID="<< EventID << " MULTIPLICITY dNch/dEta=" << dNchdEta << endl;
    }
    
    // SETUP HYDRO ATTRACTORS FOR DILEPTON CALCULATION //
    HydroAttractor::Setup();
    
    // DILEPTON KINEMATICS //
    double qTMin=0.001; double qTMax=10.0; double TauMin=0.0; double TauMax=40.0;
    
    // DILEPTON RATES //
    double dNlldQdY=0.0; double dNlldQdYPreEq=0.0; double dNlldQdYHydro=0.0;
    
    int NSamples=50000;
    CommandlineArguments.Getval("NSamples",NSamples);
    
    // SET Q RANGE AND Q BINNING //
    double QMin=1.0; double QMax=5.0;
    
    CommandlineArguments.Getval("QMin",QMin);
    CommandlineArguments.Getval("QMax",QMax);
    
    int NQ=15;
    CommandlineArguments.Getval("NQ",NQ);

    std::string fn = "ID" + std::to_string(EventID) + ".txt";
    std::ofstream OutStream; OutStream.open(fn.c_str());

    
    std::cerr << "#CALCULATING DILEPTON PRODUCTION FOR A1=" << A1 << ", A2=" << A2 << ", AT SQRT(SNN)=" << SqrtSNN << ", y=" << yRap << " AND Eta/s=" << EtaOverS << std::endl;
    std::cerr << "#KINEMATIC CUTS ARE qT=" << qTMin << " - " << qTMax << " AND  tau=" << TauMin << "-" << TauMax << " fm" << std::endl;
    
    for(int iQ=0;iQ<NQ;iQ++){
        
        double Q=QMin+iQ*(QMax-QMin)/double(NQ-1);
        double dN,dNPreEq,dNHydro;
        
        #pragma omp parallel for reduction (+ : dNlldQdY,dNlldQdYPreEq,dNlldQdYHydro)
        for(int i=0;i<NSamples;i++){
            
            DileptonEstimator::SampledNdQdy_w_intial_e(Q*Q,qTMin,qTMax,TauMin,TauMax,yRap,dEdyd2b,dN,dNPreEq,dNHydro,Ns,SaturationModel::aS,EtaOverS);
            dNlldQdY+=dN; dNlldQdYPreEq+=dNPreEq; dNlldQdYHydro+=dNHydro;
        }
        
        dNlldQdY/=double(NSamples);
        dNlldQdYPreEq/=double(NSamples);
        dNlldQdYHydro/=double(NSamples);
        
        OutStream << Q << " " << dNlldQdY << " " << dNlldQdYPreEq << " " << dNlldQdYHydro << std::endl;
    }
    OutStream.close();
    
    //FINALIZE MPI
    #if(ENABLE_MPI==1)
    MPI_Finalize();
    #endif
    
    //EXIT
    exit(0);
    
}
