function [FRC_lateral,kr_FRC_arr,kr_cutoff] = FourierRingCorrelation(img_1,img_2,pixel_dim)

    %Function to input two images and return the Fourier Ring Correlation and
    %1/7 cutoff frequency
    
    %Inputs
    %img_1 - the first image
    %img_2 - the second image
    %pixel_dim - the lateral pixel dimensions (in microns)
    
    %Outputs
    %FRC_lateral - the magnutide of the FRC values computed along progressively
    %increasing spatial frequencies as defined in: https://en.wikipedia.org/wiki/Fourier_shell_correlation
    
    %kr_FRC_arr - a vector containing the corresponding spatial frequency
    %values for FRC_lateral. The units will be in 1/microns
    
    %kr_cutoff - the spatial frequency at which FRC_lateral drops below a 1/7
    %cutoff threshold. The units will be in 1/microns
    
    % peform 2D FFT on each image
    img_1_FFT = fftshift(fft2(ifftshift(img_1)));
    img_2_FFT = fftshift(fft2(ifftshift(img_2)));
    
    % define number of pixels in each FFT image
    ky_n = size(img_1_FFT,1);
    kx_n = size(img_1_FFT,2);
    
    % calculate voxel size in each FFT images
    dk_voxel = 1./(pixel_dim.*[ky_n, kx_n]);
    
    % transform kx and ky unit in FFT images from pixels to um^-1, and center
    % at 0
    ky_vec = linspace(-floor(ky_n/2)*dk_voxel(1), floor(ky_n/2)*dk_voxel(1), ky_n);
    kx_vec = linspace(-floor(kx_n/2)*dk_voxel(2), floor(kx_n/2)*dk_voxel(2), kx_n);
    
    
    %kx_arr, ky_arr are the kx and ky in um^-1 for each pixel in FFT images
    %(hint, use the meshgrid function)
    [ky_arr,kx_arr] = meshgrid(kx_vec, ky_vec);
    
    %Compute kr_arr - an array of the radial coordinates for pixels in the FFT
    kr_arr = sqrt(kx_arr.^2 + ky_arr.^2);
    
    % calculate the maximum number of kr steps (rings) allowed for the FFT image (set to the minimum of
    % kxmax and kymax in units of pixels)
    kr_max = min([floor(ky_n/2), floor(kx_n/2)]);
    
    % create array of kr values that will be used to perform FRC calculation
    % (these are the ring radii that have been converted to units of spatial
    % frequencies (e.g. scaled by dk_voxel).
    kr_FRC_arr = (1:kr_max)*dk_voxel(1);
    
    %Initialize a vector to hold the FRC values for each ring
    FRC_lateral = zeros([kr_max, 1]);
    
    %Compute the numerator of the FRC for all pixels in the FFT's
    numerator = img_1_FFT.*conj(img_2_FFT);
    
    %Compute the two terms in the denominator of the FRC for all pixels in the FFT's
    denominator_F1 = abs(img_1_FFT).^2;
    denominator_F2 = abs(img_2_FFT).^2;
    
    %Set a tolerance to define the ring thickness (predefined for you at 0.02);
    tolerance = 0.02;
    
    %For loop to compute the FRC values at progressivly increasing ring
    %diameters
    
    for j =1:round(size(kr_FRC_arr,2))
        [num2str(j/round(size(kr_FRC_arr,2))) ' fraction complete']
        
       %Set the ring radius from the kr_FRC_array
        kr = kr_FRC_arr(j);
        
        % create the Fourier ring of the corresponding kr +- tolerance
        FRC_mask =  abs(kr_arr/kr-1)<tolerance;
        
        % calculate the FRC all pixels within the fourier
        % ring using the FRC_mask and the values for numerator and demoninator
        % computed above
        FRC_lateral(j) = sum(FRC_mask.*numerator, 'all') / (sqrt(sum(FRC_mask.*denominator_F1, 'all') .* sum(FRC_mask.*denominator_F2, 'all')));
    
    end
    
    % plot the FRC calculation results
    figure
    plot(kr_FRC_arr,abs(FRC_lateral))
    hold on
    plot([0 kr_FRC_arr(end)],[1/7 1/7]);
    kr_cutoff = min(kr_FRC_arr(abs(FRC_lateral)<1/7));
    title(['kr = ' num2str(kr_cutoff), ' (1/um)']);
    xlabel(['Spatial Frequency (1/microns)']);
    ylabel(['FRC amplitude']);