function xzPSFOverall = getOverallPSF(xzPSFexc)
    xzPSFdet = ScaleDetectionPSF;
    xzPSFOverall = xzPSFdet .* xzPSFexc;
    
    
    

