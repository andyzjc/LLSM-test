% normal pupil
k_wave = 1/wavelength_exc;
k_bound = k_xz_scale * k_wave; 
deltak = 2 * k_bound / N; 
deltax = 1/(2 * k_bound);

% excitation
[ax, az] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
kx_exc = deltak * ax;  
kz_exc = kx_exc';

ky_exc = sqrt(k_wave^2 - kx_exc.^2 - kz_exc.^2);
ky_exc(kx_exc.^2 + kz_exc.^2 > k_wave.^2 ) = 0;

x_exc = deltax * ax; 
z_exc = x_exc'; 
y_exc = (  -(N-1)/2 : (N-1)/2 ) * deltax * y_scale;
KY_exc = (-(N-1)/2 : (N-1)/2) * 1/(2*max(y_exc)) / (2*k_wave);

% for displaying
KX_exc = kx_exc(1,:) / (2*k_wave);
KZ_exc = KX_exc';
X_exc = x_exc(1,:)  / wavelength_exc; % value * wavelength = physical value (um)
Z_exc = X_exc'; 
Y_exc = y_exc / wavelength_exc;  

%% high resolution pupil
% k_wave = 1/wavelength_exc;
% NA_bound = 0.60;
% k_bound = (NA_bound / n * k_wave)+0.1; % set maximum pupil size at 0.61 NA
% deltak = 2 * k_bound / N; 
% deltax = 1/(2 * k_bound);
% 
% excitation
% [ax, az] = meshgrid(  -(N-1)/2 : (N-1)/2 ) ; 
% kx_exc = deltak * ax;  
% kz_exc = kx_exc';
% 
% ky_exc = sqrt(k_wave^2 - kx_exc.^2 - kz_exc.^2);
% ky_exc(kx_exc.^2 + kz_exc.^2 > k_wave.^2 ) = 0;
% 
% x_exc = deltax * ax; 
% z_exc = x_exc'; 
% y_exc = (  -(N-1)/2 : (N-1)/2 ) * deltax * y_scale;
% KY_exc = (-(N-1)/2 : (N-1)/2) * deltak / (2*k_wave);
% 
% for displaying
% KX_exc = kx_exc(1,:) / (2*k_wave);
% KZ_exc = KX_exc';
% X_exc = x_exc(1,:)  / wavelength_exc; % value * wavelength = physical value (um)
% Z_exc = X_exc'; 
% Y_exc = y_exc / wavelength_exc;  