namespace MultiplicityEstimator{
    
    
    // ESTIMATOR OF dNch/dEta BASED ON dE/dEtad2b [fm^-3] //
    double GetdNchdEta(double *dEdyd2b,int Ns,double aS, double eta_over_s){
        
        // CALCULATE TOTAL ENTROPY dS/deta //
        double EntropyPreFactor=4.0/3.0*std::pow(CInfty,3.0/4.0)*std::pow(4.0*M_PI*eta_over_s,1.0/3.0)*std::pow(M_PI*M_PI/30.0*nuEff,1.0/3.0);
        
        double dSdEta = 0.0;
        
        #pragma omp parallel for reduction(+ : dSdEta)
        for (int s=0; s<Ns*Ns; s++) {
            dSdEta += EntropyPreFactor*std::pow(dEdyd2b[s], 2.0/3.0);
        }
        
        dSdEta*=aS*aS;
        
        // CALCULATE CHARGED PARTICLE MULTIPLICITY //
        return dSdEta*dNchBydS;
        
    }
    
    
}
