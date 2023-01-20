function ScaledDetectionOTF = ScaleDetectionOTF
    getParameters; %modify image parameter here
    CalculatePhysics;

    detection = load("Det_PSF_OTF_510_NA1p0_RichardsWolf.mat");
    xzOTFdet = detection.xz_OTF_RW_510nm_NA1p0;
    up_factor = 1.35; %approximation from paper 
    downscale_Image_size = round(size(xzOTFdet,1) * up_factor);
    
    xzOTFdet = imresize(xzOTFdet,[downscale_Image_size,downscale_Image_size]);
    xzOTFdet = xzOTFdet/max(max(xzOTFdet));
    Image_center = (downscale_Image_size+1)/2;
    ScaledDetectionOTF = xzOTFdet(Image_center-(N+1)/2+1:Image_center+(N+1)/2-1,...
                                  Image_center-(N+1)/2+1:Image_center+(N+1)/2-1);
    Temp = abs(fftshift(ifft2(ScaledDetectionOTF)));
    ScaledDetectionOTF = abs(fftshift(fft2(Temp)));

