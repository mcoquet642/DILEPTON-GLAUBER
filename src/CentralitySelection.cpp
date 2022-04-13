namespace CentralitySelection{
    
    std::vector<double> CentralityClasses;
    
    std::vector<double> dNchdEtaMin;
    std::vector<double> dNchdEtaMax;
    std::vector<double> dNchdEtaAvg;
    std::vector<double> dNchdEtaAvgErr;
    
    // CREATE CENTRALITY SELECTION //
    void CreateCentralityClasses(int NumberOfEvents,int Ns,double EtaOverS){
        
        // GENERATE EVENT RECORD OF ALL EVENTS //
        std::vector<double> EventRecord;
        
        // SET GLOBAL SEED //
        int GlobalSeed=time(0);
        int NumberOfLocalSamples=0;
        
        #if(ENABLE_MPI==0)
        int NumberOfLocalEvents=NumberOfEvents; int EventID=GlobalSeed;
        #endif
        
        #if(ENABLE_MPI==1)
        MPI_Bcast(&GlobalSeed,1,MPI_INT,0,MPI_COMM_WORLD);
        int NumberOfLocalEvents=NumberOfEvents/MPIBasic::NumberOfNodes; int EventID=GlobalSeed+MPIBasic::ID*NumberOfLocalEvents;
        #endif
        
        // SET ENERGY DENSITY ARRAY //
        double dEdyd2b[Ns*Ns];
        
        while(NumberOfLocalSamples<NumberOfLocalEvents){

            // GENERATE EVENT //
            SaturationModel::GenerateEvent(EventID,Ns,dEdyd2b);
            
            // CALCULATE MULTIPLICITY //
            double dNchdEta=MultiplicityEstimator::GetdNchdEta(dEdyd2b,Ns,SaturationModel::aS,EtaOverS);

            EventRecord.push_back(dNchdEta); NumberOfLocalSamples++; EventID++;
            

        }
        
        // GET GLOBAL EVENT RECORD //
        std::vector<double> GlobalEventRecord;
        
        #if(ENABLE_MPI==0)
        GlobalEventRecord=EventRecord;
        #endif
        
        #if(ENABLE_MPI==1)
        
        // SYNCHRONIZE //
        MPI_Barrier(MPI_COMM_WORLD);
        
        GlobalEventRecord.resize(NumberOfEvents);
        MPI_Allgather(EventRecord.data(),NumberOfLocalEvents, MPI_DOUBLE, GlobalEventRecord.data(),NumberOfLocalEvents, MPI_DOUBLE,MPI_COMM_WORLD);
        #endif
        
        
        // SPLIT GLOBAL EVENT RECORD INTO FOUR SAMPLES  //
        std::vector<double> EventSampleA(GlobalEventRecord.begin() + 0*GlobalEventRecord.size()/4, GlobalEventRecord.begin() + 1*GlobalEventRecord.size()/4);
        std::vector<double> EventSampleB(GlobalEventRecord.begin() + 1*GlobalEventRecord.size()/4, GlobalEventRecord.begin() + 2*GlobalEventRecord.size()/4);
        std::vector<double> EventSampleC(GlobalEventRecord.begin() + 2*GlobalEventRecord.size()/4, GlobalEventRecord.begin() + 3*GlobalEventRecord.size()/4);
        std::vector<double> EventSampleD(GlobalEventRecord.begin() + 3*GlobalEventRecord.size()/4, GlobalEventRecord.end());
        
        // SORT EVENT SAMPLES //
        std::sort(EventSampleA.begin(),EventSampleA.end(),std::greater<double>());
        std::sort(EventSampleB.begin(),EventSampleB.end(),std::greater<double>());
        std::sort(EventSampleC.begin(),EventSampleC.end(),std::greater<double>());
        std::sort(EventSampleD.begin(),EventSampleD.end(),std::greater<double>());
        
        // CREATE MULTIPLICITY HISTOGRAMS //
        dNchdEtaAvg.clear(); dNchdEtaAvgErr.clear();
        
        int MinRecord=int(NumberOfEvents/4*CentralityClasses.at(0)/100.0); int MaxRecord;
        
        for(size_t i=1;i<CentralityClasses.size();i++){
            
            MaxRecord=int(NumberOfEvents/4*CentralityClasses.at(i)/100.0);
            
            double AvgA=0.0; double AvgB=0.0; double AvgC=0.0; double AvgD=0.0;
            
            for(int s=MinRecord;s<MaxRecord;s++){
                AvgA+=EventSampleA.at(s);
            }
            for(int s=MinRecord;s<MaxRecord;s++){
                AvgB+=EventSampleB.at(s);
            }
            
            for(int s=MinRecord;s<MaxRecord;s++){
                AvgC+=EventSampleC.at(s);
            }
            for(int s=MinRecord;s<MaxRecord;s++){
                AvgD+=EventSampleD.at(s);
            }
            
            AvgA/=(MaxRecord-MinRecord);
            AvgB/=(MaxRecord-MinRecord);
            AvgC/=(MaxRecord-MinRecord);
            AvgD/=(MaxRecord-MinRecord);
            
            double Avg=(AvgA+AvgB+AvgC+AvgD)/4.0;
            double Err=std::sqrt(((AvgA-Avg)*(AvgA-Avg) + (AvgB-Avg)*(AvgB-Avg) + (AvgC-Avg)*(AvgC-Avg) + (AvgD-Avg)*(AvgD-Avg))/3.0);
            
            dNchdEtaAvg.push_back(Avg); dNchdEtaAvgErr.push_back(Err);
            
            MinRecord=MaxRecord;
            
        }
        
        
        // SORT EVENTS BY MULTIPLICITY AND DETERMINE BOUNDARIES OF CENTRALITY CLASSES //
        std::sort(GlobalEventRecord.begin(),GlobalEventRecord.end(),std::greater<double>());
        
        MinRecord=int(NumberOfEvents*CentralityClasses.at(0)/100.0);
        
        dNchdEtaMin.clear(); dNchdEtaMax.clear();
        
        for(size_t i=1;i<CentralityClasses.size();i++){
            
            MaxRecord=int(NumberOfEvents*CentralityClasses.at(i)/100.0);
            
            dNchdEtaMax.push_back(GlobalEventRecord.at(MinRecord));
            dNchdEtaMin.push_back(GlobalEventRecord.at(MaxRecord-1));
            
            MinRecord=MaxRecord;
            
        }
        
        #if(ENABLE_MPI==1)
        if(MPIBasic::ID==0){
        #endif
        
            for(size_t i=1;i<CentralityClasses.size();i++){

                std::cerr << CentralityClasses.at(i-1) << " - " << CentralityClasses.at(i) << " " << dNchdEtaMin.at(i-1) << " " << dNchdEtaMax.at(i-1) << " " << dNchdEtaAvg.at(i-1) << " " << dNchdEtaAvgErr.at(i-1) << std::endl;
                
            }
            
        #if(ENABLE_MPI==1)
        }
        #endif
        
    }
    
    // GET CENTRALITY CLASS FOR AN EVENT WITH GIVEN MULTIPLICITY //
    int GetCentralityClass(double dNchdEta){
        
        // CHECK IF MULTIPLICITY IS IN THE MOST CENTRAL BIN //
        if(dNchdEta>=dNchdEtaMin.at(0)){
            return 0;
        }
        
        // CHECK IF MULTIPLICITY IS IN THE MOST PERIPHERAL BIN //
        else if(dNchdEta<dNchdEtaMin.at(dNchdEtaMin.size()-2)){
            return dNchdEtaMin.size()-1;
        }
        
        // CHECK IF MULTIPLICITY IS IN ONE OF THE INTERMEDIATE BINS //
        else{
            for(size_t i=1;i<dNchdEtaMin.size()-1;i++){
                if(dNchdEta>=dNchdEtaMin.at(i) && dNchdEta<dNchdEtaMin.at(i-1)){
                    return i;
                }
            }
        }
        
        // GIVE INVALID RETURN IF THIS IS NOT THE CASE //
        return -1;
        
    }
    
    
    // CREATE OUTPUT OF CENTRALITY CLASSES //
    void Output(std::string fname){
        
        std::ofstream OutStream; OutStream.open(fname.c_str());
        
        for(size_t i=1;i<CentralityClasses.size();i++){
            
            OutStream << CentralityClasses.at(i-1) << " - " << CentralityClasses.at(i) << " " << dNchdEtaMin.at(i-1) << " " << dNchdEtaMax.at(i-1) << " " << dNchdEtaAvg.at(i-1) << " " << dNchdEtaAvgErr.at(i-1) << std::endl;
            
        }
        
        OutStream.close();
    }
    
    
    
}
