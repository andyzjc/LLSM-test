function ScaledPSFdet = getDetectionPSF
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    lambda_det = 0.510; % um 
    wavelength_det =  lambda_det / n;
    NAdet = 1;
    
    k_wave_det = 1/wavelength_det;
    k_bound_det = k_xz_scale * k_wave_det; 
    deltak_det = 2 * k_bound_det / N;
    deltax_det = 1/(2 * k_bound_det);
    k_det = k_wave_det * NAdet / n;
    
    %xy
    [ax, ~] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
    kx_det = deltak_det * ax; 
    ky_det = kx_det'; 
    kz_det = sqrt(k_wave_det^2 - kx_det.^2 - ky_det.^2);
    kz_det(kx_det.^2 + ky_det.^2 > k_wave_det^2) = 0;
    x_det = deltax_det * ax;
    y_det = x_det';
    
    %z-propagation
    z_det = (  -(N-1)/2 : (N-1)/2  )  * deltax_det ; % no scaling factor
    KZ_det = (-(N-1)/2 : (N-1)/2 ) * 1/(2*max(z_det)) / (2*k_wave_det);
    
    % for displaying
    KX_det = kx_det(1,:) / (2 * k_wave_det);
    KY_det = KX_det';
    X_det = x_det(1,:)  / wavelength_det;
    Y_det = X_det';
    Z_det = z_det / (2*wavelength_det);
    
    %%
    Pupil = (k_det).^2 >= (kx_det.^2 + ky_det.^2);
    Pupil = Pupil.* k_wave_det./kz_det;
    Pupil(Pupil == Inf) = 0;
    Pupil = fillmissing(Pupil,'constant',0);
    PSFdet = zeros(N,N,N);
    
    for i = 1:length(z_det)
        propagator_det = exp( 2*pi * 1i * kz_det * z_det(i));
        PSFdet(:,:,i) = abs( fftshift( ifft2(ifftshift(Pupil.* propagator_det)) ) ).^2;
    end 
    % PSFdet = PSFdet/max(max(max(PSFdet)));
    
    %% resize  
    up_factorZ = 0.510/0.488; %approximation from paper 
    up_Image_sizeZ =round(N * up_factorZ)+1;
    up_factorXY = 0.510/0.488; %approximation from paper 
    up_Image_sizeXY = round(N * up_factorXY)+1;
    ZScaledPSFdet = zeros(up_Image_sizeZ,N,N);
    
    %rescale z
    for i = 1:N
        xzPSFdet = squeeze(PSFdet(:,i,:))'; 
        ZScaledPSFdet(:,i,:) = imresize(xzPSFdet,[up_Image_sizeZ,N]); % z,x,y
    end
    Image_centerZ = (up_Image_sizeZ+1)/2;
    ZScaledPSFdet = ZScaledPSFdet(Image_centerZ-(N+1)/2+1:Image_centerZ+(N+1)/2-1,:,:);
    
    %rescale xy
    XYZScaledPSFdet = zeros(N,up_Image_sizeXY,up_Image_sizeXY);
    for i = 1:N
        xyPSFdet = squeeze(ZScaledPSFdet(i,:,:)); 
        XYZScaledPSFdet(i,:,:) = imresize(xyPSFdet,[up_Image_sizeXY,up_Image_sizeXY]); % z,x,y
    end
    Image_centerXY = (up_Image_sizeXY+1)/2;
    XYZScaledPSFdet = XYZScaledPSFdet(:,Image_centerXY-(N+1)/2+1:Image_centerXY+(N+1)/2-1,Image_centerXY-(N+1)/2+1:Image_centerXY+(N+1)/2-1);
    % XYZScaledPSFdet = XYZScaledPSFdet/max(max(max(XYZScaledPSFdet)));
    
    ScaledPSFdet = XYZScaledPSFdet;
    
%     xyPSFdet = squeeze(ScaledPSFdet((N+1)/2,:,:));
%     xzPSFdet = squeeze(ScaledPSFdet(:,:,(N+1)/2)); xzPSFdet = xzPSFdet/max(max(xzPSFdet));
%     yzPSFdet = squeeze(ScaledPSFdet(:,(N+1)/2,:)); yzPSFdet = yzPSFdet/max(max(yzPSFdet));
%     zPSFdet = xzPSFdet(:,(N+1)/2); 
%     xzOTFdet = abs(fftshift(fft2(xzPSFdet))); xzOTFdet = xzOTFdet / max(max(xzOTFdet));
%     zOTFdet = xzOTFdet(:,(N+1)/2);
