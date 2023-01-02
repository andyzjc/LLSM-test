function PSFdet = getDetectionPSF
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
    z_det = (  -(N-1)/2 : (N-1)/2  )  * deltax_det * y_scale; % no scaling factor
    KZ_det = (-(N-1)/2 : (N-1)/2 ) * 1/(2*max(z_det)) / (2*k_wave_det);
    
    % for displaying
    KX_det = kx_det(1,:) / (2 * k_wave_det);
    KY_det = KX_det';
    X_det = x_det(1,:)  / wavelength_det;
    Y_det = X_det';
    Z_det = z_det / wavelength_det;
  
    Pupil = (k_det).^2 > (kx_det.^2 + ky_det.^2);
    Pupil = Pupil.* k_wave_det./kz_det;
    Pupil(Pupil == Inf) = 0;
    Pupil = fillmissing(Pupil,'constant',0);

    for i = 1:length(z_det)
        propagator_det = exp( 2*pi * 1i * kz_det * z_det(i));
        PSFdet(:,:,i) = abs( fftshift( ifft2(Pupil.* propagator_det) ) ).^2;
    end  

    PSFdet = PSFdet / max(max(max(PSFdet)));

