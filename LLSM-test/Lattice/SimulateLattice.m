function [LatticePSF,LatticePSFDithered] = SimulateLattice(LatticePupil)

getParameters;
CalculatePhysics;

LatticePSF = zeros(N,N, N);
LatticePSFDithered = zeros(N,N, N);

% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
    LatticePSF(:,:,i) = abs( fftshift( ifft2(LatticePupil .* propagator_exc) ) ).^2;

    LatticePSFDithered(:,:,i) = meshgrid(mean(LatticePSF(:,:,i),2))';
end  

% Normalize
LatticePSF = LatticePSF/max(max(max(LatticePSF)));
LatticePSFDithered = LatticePSFDithered/max(max(max(LatticePSFDithered)));