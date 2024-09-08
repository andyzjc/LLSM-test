function [LatticePSF,LatticePSFDithered,center] = SimulateLattice(LatticePupil)

getParameters;
CalculatePhysics;

% temp = zeros(1,N);
% temp(1:(N+1)/2-1) = 1./(1+exp(-(-(N+1)/4+1:(N+1)/4-1))); 
% temp((N+1)/2:end) = flip(1./(1+exp(-(-(N+1)/4:(N+1)/4-1))));
% [boundx, ~] = meshgrid(temp);
% Sample = fftshift(fft2(ifftshift(LatticePupil)));
% Sample = Sample.*boundx;
% LatticePupil = fftshift(ifft2(ifftshift((Sample))));
% LatticePupil = LatticePupil/max(LatticePupil,[],'all');

LatticePSF = zeros(N,N, N);
LatticePSFDithered = zeros(N,N, N);

% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
    LatticePSF(:,:,i) = abs( fftshift( ifft2(ifftshift(LatticePupil .* propagator_exc)) ) ).^2;
    LatticePSFDithered(:,:,i) = meshgrid(mean(LatticePSF(:,:,i),2))';
end  

% find maximum 
[center(1,1),center(1,2)] = max(LatticePSF,[],'all'); % value,  index
[center(2,1),center(2,2)] = max(LatticePSFDithered,[],'all'); % value, index

