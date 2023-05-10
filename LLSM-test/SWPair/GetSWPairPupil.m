function [SWPupil,SWMask,SWPupilMeta] = GetSWPairPupil(ProfileType,...
                                  NA1Ideal,NA2Ideal,...
                                  deltaNA1,deltaNA2,...
                                  WeightRatio)
    % generates a pair SW pupil functions 
    % ProfileType = 'gaussian','tophat'
    % WeightRatio = I(NA1)/I(NA2)
    % SWPupil and SWMask = NxNx2
    % SWPupil(:,:,1) = NA1Pupil, SWMask(:,:,1) = SWMask
theta = [90,270];
NAideal = [NA1Ideal,NA2Ideal];
deltaNA = [deltaNA1,deltaNA2];
NAmax = NAideal + deltaNA/2;
NAmin = NAideal - deltaNA/2;

getParameters;
CalculatePhysics;

% Define pupil and mask
SWPupil = zeros(N,N,2);
SWMask = zeros(N,N,2);
k_ideal = k_wave * NAideal / n;
k_deltaNA = deltaNA ./ n * k_wave;
k_NAmax = NAmax /n * k_wave; % k
k_NAmin = NAmin /n * k_wave; 

kxposition1 = k_ideal(1) .* cosd(theta) ./deltak; % pixel
kzposition1 = k_ideal(1) .* sind(theta) ./deltak; % pixel
deltaNApixels1 = k_deltaNA(1) ./ deltak;

kxposition2 = k_ideal(2) .* cosd(theta) ./deltak; % pixel
kzposition2 = k_ideal(2) .* sind(theta) ./deltak; % pixel
deltaNApixels2 = k_deltaNA(2) ./ deltak;

SWMask(:,:,1) = ((k_NAmax(1) > sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_NAmin(1) < sqrt(kx_exc.^2 + kz_exc.^2)));
SWMask(:,:,2) = ((k_NAmax(2) > sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_NAmin(2) <= sqrt(kx_exc.^2 + kz_exc.^2)));

if contains(ProfileType,'gaussian')
    % NA1 Pupil
    gaussian1 = exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA(1)/2).^2) );
    for j = 1:length(kxposition1)
        SWPupil( ...
            (N+1)/2 + round(kzposition1(j)),...
            (N+1)/2 + round(kxposition1(j)),1) = 1;
    end
    SWPupil(:,:,1) = conv2(SWPupil(:,:,1),gaussian1,'same');
%     SWPupil(:,:,1) = SWPupil(:,:,1).* SWMask(:,:,1);

    % NA2 Pupil
    if NAmin(2) > 0
        gaussian2 = exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA(2)/2).^2) );
        for j = 1:length(kxposition2)
        SWPupil( ...
            (N+1)/2 + round(kzposition2(j)),...
            (N+1)/2 + round(kxposition2(j)),2) = 1;
        end
        SWPupil(:,:,2) = (1 ./ WeightRatio) .* conv2(SWPupil(:,:,2),gaussian2,'same');
%         SWPupil(:,:,2) = SWPupil(:,:,2).* SWMask(:,:,2);
    else
        % an airy beam
        gaussian2 = exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA(2)).^2) );
        SWPupil( (N+1)/2,(N+1)/2, 2) = 1;
        SWPupil(:,:,2) = (1 ./ WeightRatio) .* conv2(SWPupil(:,:,2),gaussian2,'same');
    end
elseif contains(ProfileType,'tophat') %tophat doesnt need mask -> beam width define by deltaNA, following the non-diffractive condition
    % NA1 Pupil
    for j = 1:length(kxposition1)
    SWPupil( ...
        (N+1)/2 + round(kzposition1(j) - deltaNApixels1/2) : (N+1)/2 + round(kzposition1(j)+ deltaNApixels1/2),...
        (N+1)/2 + round(kxposition1(j)),1 ) = 1;
    end

    % NA2 Pupil
    for j = 1:length(kxposition2)
    SWPupil( ...
        (N+1)/2 + round(kzposition2(j) - deltaNApixels2/2) : (N+1)/2 + round(kzposition2(j)+ deltaNApixels2/2),...
        (N+1)/2 + round(kxposition2(j)),2 ) = 1 / WeightRatio;
    end
else
    error("Incorrect Intensity Profile")
end
SWPupil(SWPupil == inf) = 0;
SWPupil = fillmissing(SWPupil,'constant',0);

SWPupilMeta.NA1 = NA1Ideal;
SWPupilMeta.NA2 = NA2Ideal;
SWPupilMeta.deltaNA1 = deltaNA1;
SWPupilMeta.deltaNA2 = deltaNA2;
SWPupilMeta.NA1max = NAmax(1);
SWPupilMeta.NA2max = NAmax(2);
SWPupilMeta.NA1min = NAmin(1);
SWPupilMeta.NA2min = NAmin(2);
SWPupilMeta.WeightingRatio = WeightRatio;

% apply sigmoid bounding to x and z sample based on pupil
% temp = zeros(1,N);
% temp(1:(N+1)/2-1) = 1./(1+exp(-(-(N+1)/4+1:(N+1)/4-1)));
% temp((N+1)/2:end) = flip(1./(1+exp(-(-(N+1)/4:(N+1)/4-1))));
% [boundx, ~] = meshgrid(temp);

% Sample = fftshift(fft2(ifftshift(SWPupil(:,:,1))));
% Sample = Sample.*boundx;
% SWPupil(:,:,1) = fftshift(ifft2(ifftshift((Sample))));
% SWPupil(:,:,1) = SWPupil(:,:,1)/(max(SWPupil(:,:,1),[],'all'));
% SWPupil(:,:,1) = SWPupil(:,:,1).* WeightRatio;
% 
% Sample = fftshift(fft2(ifftshift(SWPupil(:,:,2))));
% Sample = Sample.*boundx;
% SWPupil(:,:,2) = fftshift(ifft2(ifftshift((Sample))));
% SWPupil(:,:,2) = SWPupil(:,:,2)/(max(SWPupil(:,:,2),[],'all'));

