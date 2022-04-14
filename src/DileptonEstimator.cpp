namespace DileptonEstimator{
    
    // SAMPLE DILPETON PRODUCTION -- QSqr INVARIANT MASS SQUARE, qT TRANSVERSE MOMENTUM , EtaQ RAPIDITY OF DILEPTON PAIR //
    void SampledNdQdy_w_intial_e(double QSqr,double qTMin,double qTMax,double TauMin,double TauMax,double EtaQ,double *dEdyd2b,double &dN,double &dNPreEq,double &dNHydro,int Ns,double aS,double eta_over_s){
        
        // SAMPLE INTEGRATION POINT //
        double Jacobian=1.0;
        
        double EtaMin=-8+EtaQ;
        double EtaMax=8+EtaQ;
        double EtaX=EtaMin+(EtaMax-EtaMin)*RandomNumberGenerator::rng();
        double qT=qTMin+(qTMax-qTMin)*RandomNumberGenerator::rng();
        double PhiQ=2.0*M_PI*RandomNumberGenerator::rng();
        
        double Q=sqrt(QSqr);
        
        // ENERGY AND MOMENTUM OF DILEPTON PAIR //
        double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
        double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
        double q0=std::sqrt(QSqr+qAbs*qAbs);
        
        // PSEUDO-RAPIDITY OF DILEPTON PAIR //
        double yQ=atanh(qZ/qAbs);
        
        // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
        Jacobian*=2.0*M_PI*(qTMax-qTMin)*(EtaMax-EtaMin)*qT*Q;
        
        // EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*RandomNumberGenerator::rng();
        Jacobian*=Tau*(TauMax-TauMin)/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);
        
        // GET EVOLUTION OF MACROSCOPIC FIELDS //
        double T,wTilde,e,pL,eQOvereG;
        
        dN=0.0; dNHydro=0.0;dNPreEq=0.0;
        
        for (int i=0; i<Ns*Ns; i++) {
            
            // GET BACKGROUND EVOLUTION AT EACH POINT //
            HydroAttractor::GetLocalValues(dEdyd2b[i],eta_over_s,nuEff,Tau,T,wTilde,e,pL,eQOvereG);
            if (T<1e-12){
                continue;
            }
            
            // CALCULATE DILEPTON PRODUCTION //
            double Xi,Teff,qSupp;
            PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);
            double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
            double dNlld4xd4Q=PreFactor*DileptonRates::SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);
            
            double dNSamp=Jacobian*dNlld4xd4Q*aS*aS;
            
            // DISTINGUSIH BETWEEN PRE-EQUILIBRIUM AND HYDRO //
            dN+=dNSamp;
            
            if(wTilde<1.0){
                dNPreEq+=dNSamp;
            }
            else{
                dNHydro+=dNSamp;
            }
            
        }
        
    }
    
    // SAMPLE DILPETON PRODUCTION -- QSqr INVARIANT MASS SQUARE, qT TRANSVERSE MOMENTUM , EtaQ RAPIDITY OF DILEPTON PAIR //
    void SampledNd2qTdQdy_w_intial_e(double Q,double qT,double TauMin,double TauMax,double EtaQ,double *dEdyd2b,double &dN,double &dNPreEq,double &dNHydro,int Ns,double aS,double eta_over_s){
        
        // GET EVOLUTION OF MACROSCOPIC FIELDS //
        double T,wTilde,e,pL,eQOvereG;

        dN=0.0; dNHydro=0.0;dNPreEq=0.0;
        
        // SET JACOBIAN FOR MONTE-CARLO INTEGRATION //
        double Jacobian=1.0;
        
        double EtaMin=-8+EtaQ;
        double EtaMax=8+EtaQ;

        // SAMPLE INTEGRATION POINT //
        double EtaX=EtaMin+(EtaMax-EtaMin)*RandomNumberGenerator::rng();
        double PhiQ=2.0*M_PI*RandomNumberGenerator::rng();
        double QSqr=Q*Q;
                    
        // ENERGY AND MOMENTUM OF DILEPTON PAIR //
        double qZ=std::sqrt(QSqr+qT*qT)*sinh(EtaQ);
        double qAbs=std::sqrt(qT*qT+(QSqr+qT*qT)*sinh(EtaQ)*sinh(EtaQ));
        double q0=std::sqrt(QSqr+qAbs*qAbs);
                    
        // PSEUDO-RAPIDITY OF DILEPTON PAIR //
        double yQ=atanh(qZ/qAbs);
                    
        // SAMPLE EVOLUTION TIME //
        double Tau=TauMin+(TauMax-TauMin)*RandomNumberGenerator::rng();
        Jacobian*=Tau*(TauMax-TauMin)/(M_HBARC*M_HBARC*M_HBARC*M_HBARC);
        // JACOBIAN  -- d^4Q=QdQ dy d^2qT //
        Jacobian*=(EtaMax-EtaMin)*Q*qT*2*M_PI;
        
        for (int i=0; i<Ns*Ns; i++) {
            
            // GET BACKGROUND EVOLUTION AT EACH POINT //
            HydroAttractor::GetLocalValues(dEdyd2b[i],eta_over_s,nuEff,Tau,T,wTilde,e,pL,eQOvereG);
            if (T<1e-12){
                continue;
            }
            
            // CALCULATE DILEPTON PRODUCTION dN/dydQd2qT FOR ALL BINS //
                    // CALCULATE DILEPTON PRODUCTION //
                    double Xi,Teff,qSupp;
                    PhaseSpaceDistribution::GetPhaseSpaceDistributionParameters(e,pL,eQOvereG,Xi,Teff,qSupp);
                    double PreFactor=alphaEM*alphaEM/(6.0*M_PI*M_PI*M_PI*QSqr)*(1.0+2.0*mllSqr/QSqr)*sqrt(1.0-4.0*mllSqr/QSqr)*qFSqrSum;
                    double dNlld4xd4Q=PreFactor*DileptonRates::SampleTracePi(q0,qT,PhiQ,yQ,EtaX,Xi,Teff,qSupp);
                    
                    double dNSamp=Jacobian*dNlld4xd4Q*aS*aS;
                    
                    // DISTINGUSIH BETWEEN PRE-EQUILIBRIUM AND HYDRO //

                    dN+=dNSamp;
                    
                    if(wTilde<1.0){
                        dNPreEq+=dNSamp;
                        dNHydro+=0;
                    }
                    else{
                        dNHydro+=dNSamp;
                        dNPreEq+=0;
                    }

            
        }
        
    }
    
                    
    
    
}
