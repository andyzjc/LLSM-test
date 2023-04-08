function [FSPSF,center] = SimulateFS(FSPupil)
    % Stimulate 3D PSF of pair of SW coherently and incoherently
    % SWPupil = NxNx2, SWPupil(:,:,1) = Pupil of SW NA1 (outer)

getParameters;
CalculatePhysics;

FSPSF = zeros(N,N,N);

% % propagation
% for j = 1:size(FSPupil,3) 
%     Pupil = FSPupil(:,:,j);
%     for i = 1:length(y_exc)
%         propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
%         temp(:,:,i) = abs( fftshift( ifft2(ifftshift(Pupil .* propagator_exc)) ) ).^2;
%     end
%     FSPSF = FSPSF + temp;
% end

% propagation
box = zeros(N,N);
box(:,(N+1)/8:7*(N+1)/8) = 1;
for j = 1:size(FSPupil,3) 
    Pupil = FSPupil(:,:,j);
    for i = 1:length(y_exc)
        propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
        temp = fftshift( ifft2(ifftshift(Pupil .* propagator_exc)) ) ;
        temp = box .* temp;
        newPupil = fftshift(fft2(ifftshift(temp)));
        temp2(:,:,i) = abs( fftshift( ifft2(ifftshift(newPupil .* propagator_exc)) ) ).^2;
    end
    FSPSF = FSPSF + temp2;
end

[center(1,1),center(1,2)] = max(FSPSF,[],'all'); % value, index