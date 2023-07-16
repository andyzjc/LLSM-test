%% 
clear all
close all

n1 = 1.33;
n2 = 1.52;
N = 513;
NA = 1.0;
wavelength_exc = 0.515; % um
d = 170; % cover glass thickness, um
gamma = 60;
scale = 4;

cartesin_array = -(N-1)/2 : (N-1)/2 ;
[ax,ay] = meshgrid(cartesin_array,cartesin_array);

kx = ax;
ky = ay;

mask = sqrt(kx.^2+ky.^2) <= (NA/n1* ((N+1)/scale)); % set NA
kx = kx .* mask;
ky = ky .* mask;

deltaNA_pixel = scale*n1/N;
NAx = cartesin_array.*deltaNA_pixel;
NAy = NAx';

kx = kx./ ((N+1)/2);
ky = ky./ ((N+1)/2);
kz = sqrt(1-kx.^2-ky.^2);
kz = kz.*mask;

k_bound = max(max(ax))./ ((N+1)/2);
deltax = 1/(2*k_bound);
x = cartesin_array.*deltax;
y = x';
z = x;

phi = acos( sind(gamma) .* kx + cosd(gamma) .* sqrt(1-kx.^2-ky.^2) );
phi_prime = asin(n1./n2 .* sin(phi));

dOPL = n2 .* d .* 1./cos(phi_prime) - n1 .* d .* 1./(cos(phi)) + n1 .* d .* sin(phi) .* ( tan(phi) - tan(phi_prime) ) ;
dOPL_substract_piston = dOPL - (n2-n1).*d;

coeff = 1.41;
dOPL_substract_defocus = dOPL_substract_piston - coeff.*1/2.*d.*n1*(1-n1/n2).*(kx.^2+ky.^2);
dOPL_substract_defocus(mask==0) = 0;

temp = 0.1 .* dOPL_substract_defocus;
phase = temp ./ wavelength_exc * 2*pi;
pupil = exp(1i.*phase);
pupil(mask==0) = 0;

% pupil = double(mask);

PSFdet = zeros(N,N,N);
for i = -(N-1)/2 : (N-1)/2
    propagator_det = exp( 2*pi * 1i * kz * i);
    PSFdet(:,:,i+(N+1)/2) = abs( fftshift( ifft2(ifftshift(pupil.* propagator_det)) ) ).^2;
end 
PSFdet = PSFdet/max(max(max(PSFdet))); 

%%
subplot(2,3,1)
imagesc(NAx,NAy,dOPL)
xlabel("NA")
ylabel("NA")
title("dOPL, unit um")
colorbar
axis image

subplot(2,3,2)
imagesc(NAx,NAy,phase)
xlabel("NA")
ylabel("NA")
title("phase,tilt=" + num2str(gamma) + "Degree")
colorbar
axis image

subplot(2,3,3)
plot(NAx,phase((N+1)/2,:))
hold on
plot(NAy,phase(:,(N+1)/2))
hold off
xlabel("NA")
legend("x","y")
axis square

subplot(2,3,4)
imagesc(x,y,squeeze(max(PSFdet,[],3)))
title("MIP x,y plane")
xlabel("x")
ylabel("y")
colorbar
clim([0,1])
axis image

subplot(2,3,5)
imagesc(y,z,squeeze(max(PSFdet,[],2)))
title("MIP y,z plane")
xlabel("y")
ylabel("z")
colorbar
clim([0,1])
axis image

subplot(2,3,6)
imagesc(x,z,squeeze(max(PSFdet,[],1)))
title("MIP x,z plane")
xlabel("x")
ylabel("z")
colorbar
clim([0,1])
axis image

