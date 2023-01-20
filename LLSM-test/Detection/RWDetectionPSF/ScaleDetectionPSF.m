function ScaledDetectionPSF = ScaleDetectionPSF
    getParameters; %modify image parameter here
    CalculatePhysics;

    detection = load("Det_PSF_OTF_510_NA1p0_RichardsWolf.mat");
    xzPSFdet = detection.xz_PSF_RW_510nm_NA1p0;
    downscale_factor = 2.35; %approximation from paper 
    downscale_Image_size = round(size(xzPSFdet,1)/downscale_factor);
    
    xzPSFdet = imresize(xzPSFdet,[downscale_Image_size,downscale_Image_size]);
    xzPSFdet = xzPSFdet/max(max(xzPSFdet));
    xzPSFdet = padarray(xzPSFdet,[1,1],0,'pre');
    desireRow = N;
    desireCol = N;
    [rows,cols] = size(xzPSFdet);
    rowsPre = floor((desireRow - rows)/2);
    collsPre = floor((desireCol - cols)/2);
    xzPSFdet = padarray(xzPSFdet, [rowsPre, collsPre], 0, 'both');
    xzPSFdet(desireRow, desireCol) = 0;
    ScaledDetectionPSF = xzPSFdet/max(max(xzPSFdet));

