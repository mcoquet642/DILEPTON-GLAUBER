namespace  NuclearData {
    
    double ORadiusInFermi=2.608;
	double CuRadiusInFermi=4.2;
    double XeRadiusInFermi=5.420;
    double AuRadiusInFermi=6.38;
    double PbRadiusInFermi=6.62;
    
    double OSurfaceInFermi=0.5130;
	double CuSurfaceInFermi=0.569;
    double XeSurfaceInFermi=0.57;
    double AuSurfaceInFermi=0.535;
    double PbSurfaceInFermi=0.546;
    
    // GET NUCLEAR RADIUS FOR WOOD SAXON POTENTIAL //
    double GetRadius(int A){
        
        // NOT NEEDED //
        if(A==1){
            return -1;
        }
        
        if(A==2){
            return -1;
        }
        
        if(A==3){
            return -1;
        }
        
        // NEEDED //
        if(A==16){
            return ORadiusInFermi;
        }
		
		if(A==64){
            return CuRadiusInFermi;
        }
        
        if(A==129){
            return XeRadiusInFermi;
        }
        
        if(A==197){
            return AuRadiusInFermi;
            
        }
        
        if(A==208){
            return PbRadiusInFermi;
        }
        
        exit(0);
        
    }
    
    // GET NUCLEAR SURFACE FOR WOOD SAXON POTENTIAL //
    double GetSurface(int A){
        
        // NOT NEEDED //
        if(A==1){
            return -1;
        }
        
        if(A==2){
            return -1;
        }
        
        if(A==3){
            return -1;
        }
        
        // NEEDED //
        if(A==16){
            return OSurfaceInFermi;
        }
		
		if(A==64){
            return CuSurfaceInFermi;
        }
        
        if(A==129){
            return XeSurfaceInFermi;
        }
        
        if(A==197){
            return AuSurfaceInFermi;
        }
        
        if(A==208){
            return PbSurfaceInFermi;
        }
        
        exit(0);
        
    }
}
