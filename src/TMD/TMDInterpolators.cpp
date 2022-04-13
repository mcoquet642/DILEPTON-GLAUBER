// GSL INTERPOLATION //



namespace Interpolator{
    
    // MIN/MAX VALUES //
    double TBMin; double TBMax;
    double TAMin; double TAMax;
    
    // VALUES //
    double *TBVals;
    double *TAVals;
    
    double *dEdyd2bVals;
    
    // GSL INTERPOLATION OBJECTS //
    gsl_interp_accel **dEdyd2bXAcc;
    gsl_interp_accel **dEdyd2bKAcc;
    
    gsl_spline2d *dEdyd2bInt;
    
    void Setup(int NTB,int NTA,double TBMinVal,double TBMaxVal,double TAMinVal,double TAMaxVal,std::string fname,int EventID){
        
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#SETTING UP dEdyd2b TABLES" << std::endl;
        
        // CREATE TABLES //
        TAVals=new double[NTA];
        TBVals=new double[NTB];
        
        dEdyd2bVals=new double[NTB*NTA];
        
        // SET TB AND TA VALUES //
        for(int ib=0;ib<NTB;ib++){
            
            // LINEAR //
            //TBVals[ib]=TBMinVal+ib*(TBMaxVal-TBMinVal)/double(NTB-1);
            
            // LOG //
            if (TBMinVal == 0) { TBMinVal = (TBMaxVal-TBMinVal)/double(NTB-1)*1e-2; }
            TBVals[ib]=TBMinVal*std::exp(ib*log(TBMaxVal/TBMinVal)/double(NTB-1))-TBMinVal;
        
        }
        
        for(int ia=0;ia<NTA;ia++){
            
            // LINEAR //
            //TAVals[ia]=TAMinVal+ia*(TAMaxVal-TAMinVal)/double(NTA-1);
            
            //LOG //
            if (TAMinVal == 0) { TAMinVal = (TAMaxVal-TAMinVal)/double(NTA-1)*1e-2; }
            TAVals[ia]=TAMinVal*std::exp(ia*log(TAMaxVal/TAMinVal)/double(NTA-1))-TAMinVal;
        
        }
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        dEdyd2bXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        dEdyd2bKAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            dEdyd2bXAcc[i]= gsl_interp_accel_alloc ();
            dEdyd2bKAcc[i]= gsl_interp_accel_alloc ();
        }
        
        // SETUP GSL SPLINES //
        dEdyd2bInt = gsl_spline2d_alloc(gsl_interp2d_bilinear,NTA,NTB);
        
        
        // CALCULATE dEdyd2b FOR TA AND TB VALUES AND TABULATE RESULTS //
        
        #pragma omp parallel for ordered schedule(dynamic, 8)
        for(int ia=0;ia<NTA;ia++){
            
            // GET TA VALUE //
            double TA=TAVals[ia];
            
            for(int ib=0;ib<NTB;ib++){
                
                // GET TB VALUE //
                double TB=TBVals[ib];
                
                // CALCULATE dEdyd2b FOR GIVEN VALUES OF TA,TB <--- CHANGE HERE //
                //std::cerr << "ia " << ia << " ib " << ib << std::endl;
                double CurrentdEdyd2bValue=Calculate_dEdyd2b(TA,TB);
                
                // WRITE TO TABLE //
                gsl_spline2d_set(dEdyd2bInt,dEdyd2bVals,ia,ib,CurrentdEdyd2bValue);
                
            }
            
        }
        
        // CREATE OUTPUT //
        std::ofstream OutStream;
        OutStream.open(fname.c_str());
        OutStream.precision(15);
        
        for(int ia=0;ia<NTA;ia++){
            for(int ib=0;ib<NTB;ib++){
                
                // GET TA,TB VALUE //
                double TA=TAVals[ia]; double TB=TBVals[ib];
                
                // GET dEdyd2b VALUE //
                double dEdyd2b=gsl_spline2d_get(dEdyd2bInt,dEdyd2bVals,ia,ib);
                
                OutStream << ia << " " << ib << " " << TA << " " << TB << " " << dEdyd2b  << std::endl;
            }
            
            OutStream << std::endl;
        }
        
        OutStream.close();
        
        // SET MIN/MAX //
        TAMin=TAVals[0]; TAMax=TAVals[NTA-1];
        TBMin=TBVals[0]; TBMax=TBVals[NTB-1];
        
        // SETUP GSL INTERPOLATORS //
        gsl_spline2d_init(dEdyd2bInt,TAVals,TBVals,dEdyd2bVals,NTA,NTB);
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#SETUP OF dEdyd2b FROM TA=" << TAMin << " TO " << TAMax << " AND  TB=" << TBMin << " TO " << TBMax << " COMPLETED" << std::endl;
        
        // CLEAN-UP //
        delete[] TAVals;
        delete[] TBVals;
        
        delete[] dEdyd2bVals;
        
    }
    
    
    // LOAD INTERPOLATOR TABLE FROM FILE //
    void LoadTable(int NTB, int NTA, std::string fname, int EventID) {
    
    	
        // COMMANDLINE OUTPUT //
        std::cerr << "#LOADING dEdyd2b TABLES" << std::endl;
        
        // CREATE TABLES //
        TAVals=new double[NTA];
        TBVals=new double[NTB];
        
        dEdyd2bVals=new double[NTB*NTA];
        
        // SETUP GSL INTERPOLATION //
        int NumberOfOpenMPThreads=omp_get_max_threads();
        
        dEdyd2bXAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        dEdyd2bKAcc=new gsl_interp_accel*[NumberOfOpenMPThreads];
        
        #pragma omp parallel for
        for(int i=0;i<NumberOfOpenMPThreads;i++){
            
            dEdyd2bXAcc[i]= gsl_interp_accel_alloc ();
            dEdyd2bKAcc[i]= gsl_interp_accel_alloc ();
        }
        
        // SETUP GSL SPLINES //
        dEdyd2bInt = gsl_spline2d_alloc(gsl_interp2d_bilinear,NTA,NTB);
        
        
        // READING IN TABLUATED RESULTS FROM FILE //
        std::cerr << "#Reading File " << fname << std::endl;
        std::ifstream Read;
        Read.open(fname.c_str());
        
        int iaRead; int ibRead; double TA; double TB; double CurrentdEdyd2bValue;
        
        for(int ia=0;ia<NTA;ia++){
            
            for(int ib=0;ib<NTB;ib++){
                
                // READING THE VALUES FROM FILE //
                
                Read >> iaRead; Read >> ibRead; Read >> TA; Read >> TB; Read >> CurrentdEdyd2bValue;
                
                TAVals[ia] = TA;
                TBVals[ib] = TB;
                
                if (iaRead != ia || ibRead != ib) { std::cerr << "Error: ia or ib from file does not match read ia or ib" << std::endl; }
                
                // WRITE TO TABLE //
                gsl_spline2d_set(dEdyd2bInt,dEdyd2bVals,ia,ib,CurrentdEdyd2bValue);
                
            }
            
        }
        
        Read.close();
        
        // SET MIN/MAX //
        TAMin=TAVals[0]; TAMax=TAVals[NTA-1];
        TBMin=TBVals[0]; TBMax=TBVals[NTB-1];
        
        // SETUP GSL INTERPOLATORS //
        gsl_spline2d_init(dEdyd2bInt,TAVals,TBVals,dEdyd2bVals,NTA,NTB);
        
        // COMMANDLINE OUTPUT //
        std::cerr << "#SETUP OF dEdyd2b FROM TA=" << TAMin << " TO " << TAMax << " AND  TB=" << TBMin << " TO " << TBMax << " COMPLETED" << std::endl;
        
        // CLEAN-UP //
        delete[] TAVals;
        delete[] TBVals;
        
        delete[] dEdyd2bVals;
        
    
    }
    
    double GetdEdyd2b(double TA,double TB){
        
        
        double INTERPOLATED_VALUE;
        
        int tID=omp_get_thread_num();
		        
        // EVALUATION OF GSL INTERPOLATOR //
        if((TA)<(TAMin) || (TA)>(TAMax) || (TB)<(TBMin) || (TB)>(TBMax)){
            
            //std::cerr << "#WARNING -- DATA BASE ACCESSED OUT OF INTERPOLATOR RANGE" << std::endl;
            //std::cerr << "#TAMin=" << TAMin << " TA=" << TA << " TAMax=" << TAMax << std::endl;
            //std::cerr << "#TBMin=" << TBMin << " TB=" << TB << " TBMax=" << TBMax << std::endl;
            
            if (TA*TB > 1e-10) {
            	
            	std::cerr << "#WARNING -- LARGE VALUE THROWN AWAY" << std::endl;
            	std::cerr << "#TAMin=" << TAMin << " TA=" << TA << " TAMax=" << TAMax << std::endl;
            	std::cerr << "#TBMin=" << TBMin << " TB=" << TB << " TBMax=" << TBMax << std::endl;
            
            }            
            
            INTERPOLATED_VALUE=0.0;
        }
        else{
            INTERPOLATED_VALUE=gsl_spline2d_eval(dEdyd2bInt,TA,TB,dEdyd2bXAcc[tID],dEdyd2bKAcc[tID]);
        }
        
        return INTERPOLATED_VALUE;

        
        
    }
}
