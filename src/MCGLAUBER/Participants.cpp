//////////////////////////////////////////////////
//                                              //
//////////////////////////////////////////////////

// COMPUTE COLLISION OVERLAP BETWEEN INDIVIDUAL NUCLEONS //
double ComputeIndividualOverlap(Nucleus *N1,int in1,Nucleus *N2,int in2,double bX,double bY){
    
    double x1=N1->x[in1][0];    double y1=N1->x[in1][1];
    double x2=N2->x[in2][0];    double y2=N2->x[in2][1];
    
    double dsqr=(x1-x2+bX)*(x1-x2+bX)+(y1-y2+bY)*(y1-y2+bY);
    
    return exp(-dsqr/(4.0*BG))/(4.0*M_PI*BG);
    
}

// GET EVENT PROBABILITY //
double ComputeInteractionProbability(double CollisionOverlap,int VERBOSE_MODE){
    
    // NUCLEON-NUCLEON CROSS SECTION //
    double SigmaNN;
    
    // DETERMINE NUCLEON-NUCLEON INELASTIC CROSS-SECTION //
    if(SqrtS_NN==13000){
        SigmaNN=53.5;
    }
    
    else if(SqrtS_NN==7000){
        SigmaNN=39.4;
    }
    
    else if(SqrtS_NN==5020){
        SigmaNN=27.7462;
    }
    
    else if(SqrtS_NN==2760){
        SigmaNN=23.0;
    }
    
    else if(SqrtS_NN==900){
        SigmaNN=15.6844;
    }
    
    else if(SqrtS_NN==200){
        SigmaNN=9.3916;
    }
    
    // CHECK THAT CROSS SECTION IS SPECIFIED //
    else{
        
        // WARNING MESSAGE //
        if(VERBOSE_MODE==1){
            std::cerr << "#WARNING -- CENTER OF MASS ENERGY " << SqrtS_NN << " NOT SPECIFIED -- USING intERPOLATION FORMULA" << std::endl;
        }
        
        double logS=log(SqrtS_NN);
        
        SigmaNN=-44.8024 + 0.000129448*logS + 0.705328*logS*logS + 181.784/logS;
        
    }
    
    // COMMANDLINE OUTPUT //
    if(VERBOSE_MODE==1){
        std::cerr << "#SIGMA NN=" << SigmaNN << " FOR CENTER OF MASS ENERGY " << SqrtS_NN << " GEV" << std::endl;
    }
    
    // COMPUTE INTERACTION PROBABILITY //
    return 1.0-exp(-SigmaNN*CollisionOverlap);
    
}

void CheckParticipants(Nucleus* N1,Nucleus *N2,double bX,double bY){
    
    // NUCLEON-NUCLEON INTERACTION STATUS //
    int **InteractionStatus= new int*[N1->A];
    
    for(int s1=0;s1<N1->A;s1++){
        InteractionStatus[s1]=new int[N2->A];
    }
    
    // COMPUTE NUCLEON-NUCLEON intERACTION PROBABILITIES AND SAMPLE intERACTION //
    for(int s1=0;s1<N1->A;s1++){
        
        for(int s2=0;s2<N2->A;s2++){
            
            // GET NUCLEON-NUCLEON OVERLAP //
            double NNOverlap=ComputeIndividualOverlap(N1,s1,N2,s2,bX,bY);
            
            // COMPUTE INTERACTION PROBABILITY //
            double NNInteractionProbability=ComputeInteractionProbability(NNOverlap,0);
            
            // CHECK THAT intERACTION PROBABILITY IS SUFFICIENTLY LARGE TO BE RELEVANT //
            if(RandomNumberGenerator::rng()<NNInteractionProbability){
                InteractionStatus[s1][s2]=1;
            }
            else{
                InteractionStatus[s1][s2]=0;
            }
            
        }
    }
    
    // CHECK FOR EACH NUCLEON IN NUCLEUS ONE THE NUMBER OF COLLISIONS AND PARTICIPANT STATUS //
    N1->NumberOfParticipants=0; N1->NumberOfCollisions=0;
    
    for(int s1=0;s1<N1->A;s1++){
        
        // NUMBER OF INDIVIDUAL COLLISIONS //
        int NIndColl=0;
        
        for(int s2=0;s2<N2->A;s2++){
            
            if(InteractionStatus[s1][s2]==1){
                NIndColl++;
            }
        }
        
        // UPDATE NUMBER OF COLLISIONS //
        N1->CollisionNumber[s1]=NIndColl; N1->NumberOfCollisions+=NIndColl;
        
        // UPDATE PARTICIPANT STATUS //
        if(NIndColl>0){
            N1->ParticipantStatus[s1]=1; N1->NumberOfParticipants++;
        }
        else{
            N1->ParticipantStatus[s1]=0;
        }
        
    }
    
    // CHECK FOR EACH NUCLEON IN NUCLEUS TWO THE NUMBER OF COLLISIONS AND PARTICIPANT STATUS //
    N2->NumberOfParticipants=0; N2->NumberOfCollisions=0;

    for(int s2=0;s2<N2->A;s2++){
        
        // NUMBER OF INDIVIDUAL COLLISIONS //
        int NIndColl=0;
        
        for(int s1=0;s1<N1->A;s1++){
            
            if(InteractionStatus[s1][s2]==1){
                NIndColl++;
            }
        }
        
        // UPDATE NUMBER OF COLLISIONS //
        N2->CollisionNumber[s2]=NIndColl; N2->NumberOfCollisions+=NIndColl;
        
        // UPDATE PARTICIPANT STATUS //
        if(NIndColl>0){
            N2->ParticipantStatus[s2]=1; N2->NumberOfParticipants++;
        }
        else{
            N2->ParticipantStatus[s2]=0;
        }
        
    }
    
    // CLEAN-UP //
    for(int s1=0;s1<N1->A;s1++){
        delete InteractionStatus[s1];
    }
    
    delete InteractionStatus;
    
    
}
