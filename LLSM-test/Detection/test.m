clear all
close all

%% 
N = 513; % image pixels
n = 1.33; % medium index of refraction
lambda_det = 0.510; % um 
wavelength_det =  lambda_det / n;
NAdet = 1;
k_xz_scale = 2; % Pupil space sampling factor
y_scale = 1; %Propagation step

k_wave = 1/wavelength_det;
k_bound = k_xz_scale * k_wave; 
deltak = 2 * k_bound / N;
deltax = 1/(2 * k_bound);
k_det = k_wave * NAdet / n;

%xy
[ax, ~] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
kx_det = deltak * ax;
ky_det = kx_det';
kz_det = sqrt(k_wave^2 - kx_det.^2 - ky_det.^2);
kz_det(kx_det.^2 + ky_det.^2 > k_wave^2) = 0;
x_det = deltax * ax;
y_det = x_det';

%z-propagation
z_det = (  -(N-1)/2 : (N-1)/2  )  * deltax * y_scale; % no scaling factor
KZ_det = (-(N-1)/2 : (N-1)/2 ) * 1/(2*max(z_det)) / (2*k_wave);

% for displaying
KX_det = kx_det(1,:) / (2 * k_wave);
KY_det = KX_det';
X_det = x_det(1,:)  / wavelength_det;
Y_det = X_det';
Z_det = z_det / wavelength_det;

Pupil = (k_det).^2 > (kx_det.^2 + ky_det.^2);
% Pupil_fun_det = Pupil_fun_det.* k_wave./ky_det;
% Pupil_fun_det(Pupil_fun_det == Inf) = 0;
% Pupil_fun_det = fillmissing(Pupil_fun_det,'constant',0);

for i = 1:length(z_det)
    propagator_det = exp( 2*pi * 1i * kz_det * z_det(i));
    PSF_det_3d(:,:,i) = abs( fftshift( ifft2(Pupil.* propagator_det) ) ).^2;
end  

xzPSFdetection = squeeze(PSF_det_3d(:,(N+1)/2,:))'; xzPSFdetection = xzPSFdetection/max(max(xzPSFdetection));
zPSFdetection = xzPSFdetection(:,(N+1)/2); zPSFdetection = zPSFdetection/max(zPSFdetection);
xzOTFdetection = abs(fftshift(fft2(xzPSFdetection))); xzOTFdetection = xzOTFdetection / max(max(xzOTFdetection));
zOTFdetection = xzOTFdetection(:,(N+1)/2); zOTFdetection = zOTFdetection/max(zOTFdetection);

xzPSFdet = xzPSFdet/max(max(xzPSFdet));
zPSFdet = xzPSFdet(:,(N+1)/2); zPSFdet = zPSFdet/max(zPSFdet);
xzOTFdet = abs(fftshift(fft2(xzPSFdet))); xzOTFdet = xzOTFdet / max(max(xzOTFdet));
zOTFdet = xzOTFdet(:,(N+1)/2); zOTFdet = zOTFdet/max(zOTFdet);

figure(1)
subplot(1,2,1)
imagesc(X_det,Z_det,xzPSFdetection);
axis image
subplot(1,2,2)
imagesc(X_det,Z_det,xzPSFdet);
axis image

figure(2)
subplot(1,2,1)
imagesc(KX_det,KZ_det,xzOTFdetection);
axis image
subplot(1,2,2)
imagesc(KX_det,KZ_det,xzOTFdet);
axis image

figure(3)
hold on
plot(Z_det,zPSFdetection)
plot(Z_det,zPSFdet)
legend("Propagation","RW")
grid on

% %% try something new
% kx = deltak * (-N-1:N-1);
% ky = deltak * (-N-1:N-1)';
% kz = sqrt(k_wave^2 - kx_det.^2 - ky_det.^2);

