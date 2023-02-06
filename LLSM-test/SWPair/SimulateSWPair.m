function [PSFCoherent,PSFIncoherent,center] = SimulateSWPair(SWPupil)
    % Stimulate 3D PSF of pair of SW coherently and incoherently
    % SWPupil = NxNx2, SWPupil(:,:,1) = Pupil of SW NA1 (outer)

getParameters;
CalculatePhysics;

Pupil1 = squeeze(SWPupil(:,:,1));
Pupil2 = squeeze(SWPupil(:,:,2));

PSFCoherent= zeros(N,N,N);
PSFIncoherent = zeros(N,N,N);
Pupil_sum = Pupil1 + Pupil2;
    % propagation
    for i = 1:length(y_exc)
        propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));

        Profile_pupil1 = abs( fftshift( ifft2(ifftshift(Pupil1 .* propagator_exc)) ) ).^2;
        Profile_pupil2 = abs( fftshift( ifft2(ifftshift(Pupil2 .* propagator_exc)) ) ).^2;
        PSFIncoherent(:,:,i) = Profile_pupil1 + Profile_pupil2;

        PSFCoherent(:,:,i) = abs( fftshift( ifft2(ifftshift(Pupil_sum .* propagator_exc)) ) ).^2;
    end

[center(1,1),center(1,2)] = max(PSFCoherent,[],'all'); % value, index
[center(2,1),center(2,2)] = max(PSFIncoherent,[],'all'); % value, index
