function fc2 = FC2(Image1,Image2)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % fourier transform 
    FTimage1 = fftshift(fft2(ifftshift(Image1)));
    FTimage2 = fftshift(fft2(ifftshift(Image2)));

    % Fourier ring correlation formula 
    frc(:,:,1) = FTimage1 .* conj(FTimage2);
    frc(:,:,2) = abs(FTimage1).^2;
    frc(:,:,3) = abs(FTimage2).^2;

    % split into 4 quadrant 
    for i = 1:size(frc,3)
        one(:,:,i) = frc(1:(N+1)/2,(N+1)/2:N,i);
        two(:,:,i) = flip(frc(1:(N+1)/2,1:(N+1)/2,i),2);
        three(:,:,i) = flip(flip(frc((N+1)/2:N,1:(N+1)/2,i),1),2);
        four(:,:,i) = flip(frc((N+1)/2:N,(N+1)/2:N,i),1);
    end

    % sum and average 
    numerator = (one(:,:,1) + two(:,:,1) + three(:,:,1)+ four(:,:,1))./4;
    denominator1 = (one(:,:,2) + two(:,:,2) + three(:,:,2)+ four(:,:,2))./4;
    denominator2 = (one(:,:,3) + two(:,:,3) + three(:,:,3)+ four(:,:,3))./4;

    numerator(abs(numerator)<10^-30) = 0;
    denominator1(abs(denominator1)<10^-30) = 0;
    denominator2(abs(denominator2)<10^-30) = 0;

    % calculate fourier correlation
    fc2 = real(numerator) ./ (sqrt(denominator1.*denominator2));
    fc2 = fc2/max(max(fc2));

%     plot
%     figure
%     subplot(2,2,1)
%     imagesc(Image1)
%     title("SNR=1000")
%     axis image
%     colorbar
% 
%     subplot(2,2,2)
%     imagesc(Image2)
%     title("SNR=10")
%     axis image
%     colorbar
% 
%     h1 = subplot(2,2,3);
%     imagesc(KX_exc((N+1)/2:N),flip(KZ_exc((N+1)/2:N)),fc2)
%     xlabel("k_x/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     colormap(hot)
%     colorbar
%     set(h1, 'YDir','normal')
%     axis image
% 
%     h2 = subplot(2,2,4);
%     fc2contour = zeros(size(fc2));
%     fc2contour(fc2>=1/7) = 1;
%     imagesc(KX_exc((N+1)/2:N),flip(KZ_exc((N+1)/2:N)),fc2contour)
%     xlabel("k_x/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     set(h2, 'YDir','normal')
%     colormap(hot)
%     colorbar
%     axis image
