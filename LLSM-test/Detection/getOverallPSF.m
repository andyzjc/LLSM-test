function xzPSFOverall = getOverallPSF(xzPSFexc)
    getParameters; %modify image parameter here
    CalculatePhysics;

    xzPSFdet = load("Det_PSF_OTF_510_NA1p0_RichardsWolf.mat");
    xzPSFdet = imresize(xzPSFdet.xz_PSF_RW_510nm_NA1p0,[N,N]);
    
    xzPSFOverall = xzPSFdet .* xzPSFexc;