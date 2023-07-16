clear all
close all

n1 = 1.33;
n2 = 1.52;
N = 513;
NA = 1.0;
alpha = asind(NA/n1);
wavelength_exc = 0.515; % um

d = 170; % cover glass thickness, um
gamma = 0;

cartesin_array = -(N-1)/2 : (N-1)/2 ;
[ax,ay] = meshgrid(cartesin_array,cartesin_array);

% define normal 
n_x = sind(gamma) .* ones(N,N);
n_y = zeros(N,N);
n_z = cosd(gamma) .* ones(N,N);

% angles in spherical coordinate
% [phi,theta,~] = cart2sph(ax,ay,az);
% theta = (pi/2-theta) .* 180./(pi); % to degree
% phi = phi .* 180 ./ (pi); 
% phi(theta > alpha) = 0; 
% theta(theta >alpha) = 0; 

phi = atan2(ay,ax) .* 180./pi;
theta = (pi/2-atan2(129,sqrt(ax.^2 + ay.^2))) .* 180./pi;
mask = theta <=alpha;
theta = mask.*theta;
theta = theta + 180;
phi = mask.*phi;

%% wave vector 
k_x = sind(theta).*cosd(phi) ; k_x = k_x.*mask;
k_y = sind(theta).*sind(phi) ; k_y = k_y.*mask;
k_z = cosd(theta); k_z = k_z.*mask;
k = sqrt(k_x.^2 + k_y.^2 + k_z.^2); % unit circle
% k_x = k_x./k;
% k_y = k_y./k;
% k_z = k_z./k;

%% angle of incidence 
sin_psi_x = (  - sind(phi).*sind(theta) .* cosd(gamma) );
sin_psi_x = sin_psi_x.* mask;
sin_psi_y = -(sind(gamma).*cosd(theta) - cosd(gamma).*sind(theta).*cosd(phi) );
sin_psi_y = sin_psi_y.* mask;
sin_psi_z = (sind(gamma).*sind(theta) .* sind(phi));
sin_psi_z =sin_psi_z .* mask;
sin_psi = sqrt( sin_psi_x.^2 + sin_psi_y.^2 + sin_psi_z.^2);

cos_psi_x = sind(gamma) .* sind(theta) .* cosd(phi);
cos_psi_x = cos_psi_x.* mask;
cos_psi_y = 0;
cos_psi_z = cosd(gamma) .* cosd(theta);
cos_psi_z = cos_psi_z.* mask;
cos_psi = abs(sqrt( cos_psi_x.^2 + cos_psi_y.^2 + cos_psi_z.^2));

%% s-p-k coordinate
s_x = 1./ sin_psi .* (sin_psi_x);
s_x = s_x.* mask;
s_y = 1./ sin_psi .* (sin_psi_y);
s_y = s_y.* mask;
s_z = 1./ sin_psi .* (sin_psi_z);
s_z = s_z.* mask;
% s = sqrt( s_x.^2 + s_y.^2 + s_z.^2);
% s_x = s_x./k;
% s_y = s_y./k;
% s_z = s_z./k;

p_x = s_y .* k_z - s_z .* k_y;
p_x = p_x.* mask;
p_y = -(s_x .* k_z - s_z .* k_x);
p_y = p_y.* mask;
p_z = s_x .* k_y - s_y .* k_x;
p_z = p_z.* mask;

%% snell's law, refracted ray 
psi_prime = asind(n1/n2 .* sin_psi);
beta = asind(sin_psi) - psi_prime;

%% refracted ray, k'
k_prime_x = cosd(beta) .* k_x - sind(beta) .* p_x; 
k_prime_x = k_prime_x.* mask;
k_prime_y = cosd(beta) .* k_y - sind(beta) .* p_y;
k_prime_y = k_prime_y.* mask;
k_prime_z = cosd(beta) .* k_z - sind(beta) .* p_z;
k_prime_z = k_prime_z.* mask;

%% upper water-cover slip interface
% incident plane
upper_r_x = -n_x;
upper_r_y = -n_y;
upper_r_z = -n_z;

%%
b = [0; 0; -d];
for i = 1:N
    for j = 1:N
        M1 = [k_prime_z(i,j) 0 -k_prime_x(i,j)];
        M2 = [0 k_prime_z(i,j) -k_prime_y(i,j)];
        M3 = [sind(gamma) 0 cosd(gamma)];
        M = [M1;M2;M3];
        if sum(sum(isnan(M))) > 0
            r_e_x(i,j) = 0;
            r_e_y(i,j) = 0;
            r_e_z(i,j) = 0;
        else
            r_e = linsolve(M,b);
            r_e_x(i,j) = r_e(1);
            r_e_y(i,j) = r_e(2);
            r_e_z(i,j) = r_e(3);
        end
    end
end
r_e_x = r_e_x .* mask;
r_e_y = r_e_y .* mask;
r_e_z = r_e_z .* mask;

%% OLC
OPL_c = n2 .* sqrt(r_e_x.^2 + r_e_y.^2 + r_e_z.^2);

delta_l_x = n1 .* sind(theta) .* (r_e_x - (r_e_z)./(k_z) .* k_x);
delta_l_x = delta_l_x .* mask;
delta_l_y = n1 .* sind(theta) .* (r_e_y - (r_e_z)./(k_z) .* k_y);
delta_l_y = delta_l_y .* mask;
delta_l_z = n1 .* sind(theta) .* (r_e_z - (r_e_z)./(k_z) .* k_z);
delta_l_z = delta_l_z .* mask;
delta_l = sqrt( delta_l_x.^2 + delta_l_y.^2 + delta_l_z.^2);

OPLr = OPL_c + delta_l;
OPLi_x = n1 .* (r_e_z ./ (k_z) .* k_x);
OPLi_x = OPLi_x .* mask;
OPLi_y = n1 .* (r_e_z ./ (k_z) .* k_y);
OPLi_y = OPLi_y .* mask;
OPLi_z = n1 .* (r_e_z ./ (k_z) .* k_z);
OPLi_z = OPLi_z .* mask;
OPLi = sqrt( OPLi_x.^2 + OPLi_y.^2 + OPLi_z.^2);

phase_shift = (OPLr - OPLi)/wavelength_exc;

% subtract piston
phase_shift = phase_shift - (n2-n1).*d;

% subtract defocus
phase_shift = phase_shift - 1/2.*d.*n1*(1-n1/n2).*(k_x.^2+k_y.^2);
phase_shift = phase_shift ./ wavelength_exc;

figure
subplot(1,3,1)
title("Phase after subtract piston, unit: 1/\lambda, tilt=" + num2str(gamma) + "Degree")
imagesc(phase_shift); colorbar; 
axis image

subplot(1,3,2)
title("horizontal lincut, y")
plot(phase_shift(:,(N+1)/2));
axis square

subplot(1,3,3)
title("horizontal lincut, x")
plot(phase_shift((N+1)/2,:));
axis square




