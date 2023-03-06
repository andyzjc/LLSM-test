function [LatticePSF,LatticePSFDithered,center] = SimulateLattice(LatticePupil)

getParameters;
CalculatePhysics;

LatticePSF = zeros(N,N, N);
LatticePSFDithered = zeros(N,N, N);

% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i .* ky_exc .* y_exc(i));
    LatticePSF(:,:,i) = abs( fftshift( ifft2(ifftshift(LatticePupil .* propagator_exc)) ) ).^2;
    LatticePSFDithered(:,:,i) = meshgrid(mean(LatticePSF(:,:,i),2))';
end  

% find maximum 
[center(1,1),center(1,2)] = max(LatticePSF,[],'all'); % value,  index
[center(2,1),center(2,2)] = max(LatticePSFDithered,[],'all'); % value, index

