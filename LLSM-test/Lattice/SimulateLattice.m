function [LatticePSF,LatticePSFDithered,center] = SimulateLattice(LatticePupil)

getParameters;
CalculatePhysics;

LatticePSF = zeros(N,N, N);
LatticePSFDithered = zeros(N,N, N);

% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
    LatticePSF(:,:,i) = abs( fftshift( ifft2(ifftshift(LatticePupil .* propagator_exc)) ) ).^2;
    LatticePSFDithered(:,:,i) = meshgrid(mean(LatticePSF(:,:,i),2))';
end  

% box = zeros(N,N);
% box(:,(N+1)/8:7*(N+1)/8) = 1;
% for i = 1:length(y_exc)
%     propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
%     temp = fftshift( ifft2(ifftshift(LatticePupil .* propagator_exc)) );
%     temp = box .* temp;
%     newPupil = fftshift(fft2(ifftshift(temp)));
%     LatticePSF(:,:,i) = abs( fftshift( ifft2(ifftshift(newPupil .* propagator_exc)) ) ).^2;
%     LatticePSFDithered(:,:,i) = meshgrid(mean(LatticePSF(:,:,i),2))';
% end  

% % real dither 
% LatticePSFDithered = LatticePSF;
% dither_step = 64;
% dither_period = 30; %um
% for j = 1:dither_step
%     LatticePSFDithered = LatticePSFDithered + ...
%         circshift(LatticePSF,round(dither_period / deltax ),2);
% end

% find maximum 
[center(1,1),center(1,2)] = max(LatticePSF,[],'all'); % value,  index
[center(2,1),center(2,2)] = max(LatticePSFDithered,[],'all'); % value, index

