function [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil)
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

        Profile_pupil1 = abs( fftshift( ifft2(Pupil1 .* propagator_exc) ) ).^2;
        Profile_pupil2 = abs( fftshift( ifft2(Pupil2 .* propagator_exc) ) ).^2;
        PSFIncoherent(:,:,i) = Profile_pupil1 + Profile_pupil2;

        PSFCoherent(:,:,i) = abs( fftshift( ifft2(Pupil_sum .* propagator_exc) ) ).^2;
    end

PSFIncoherent = PSFIncoherent/max(max(max(PSFIncoherent)));
PSFCoherent = PSFCoherent/max(max(max(PSFCoherent)));
