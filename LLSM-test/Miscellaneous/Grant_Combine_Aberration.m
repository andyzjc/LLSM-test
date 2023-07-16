%% Simulation without aberrration
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% detection 
% detection
PSFdet = getDetectionPSF;
PSFdetmax = (max(max(max(PSFdet))));
PSFdet = PSFdet./PSFdetmax;

%% Aberrated detection
NAdet = 1.0;
k_apertureNA = NAdet * k_wave / n;
Pupil = (k_apertureNA).^2 >= abs((kx_exc.^2 + kz_exc.^2));

% get single zernike function
PupilNA = 1.0;

% PhaseAmplitude = 1;
Phase1 = GetSingleZmodePupil(2,2,PupilNA);
Phase2 = GetSingleZmodePupil(3,1,PupilNA);

Phase = (6*Phase1 + 4*Phase2) .* wavelength_exc / (2*pi);
AberratedPupil = Pupil.*exp(1i.*Phase);

[AberratedPSFdet,~,~] =  SimulateLattice(AberratedPupil);
AberratedPSFdet = AberratedPSFdet/PSFdetmax;

% AberratedPSFdet = PSFdet;
xz = squeeze(AberratedPSFdet(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(AberratedPSFdet(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(AberratedPSFdet((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,flip(xy',1))
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,flip(yz',1))
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

%% LLS
PupilNA = 0.65;

NA1 = 0.58;
deltaNA=0.04;
LatticeType = 'hex';
Latticeweighting = [1 1 1 1 1 1];
Phase1 = GetSingleZmodePupil(2,2,PupilNA);
Phase2 = GetSingleZmodePupil(3,1,PupilNA);

Phase = (6*Phase1 + 4*Phase2) .* wavelength_exc / (2*pi);

[LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,'tophat', ...
                                                             NA1,deltaNA, ...
                                                             0.8,0.0,...
                                                             Latticeweighting);


%% Unaberrated, excitation
[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
LatticePSF = LatticePSF/max(LatticePSF,[],'all');
LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
% grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,LatticePSF,LatticePSFDithered,LatticePupil,AberratedPSFdet,pwd)

xz = squeeze(LatticePSF(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(LatticePSF(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(LatticePSF((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

xz = squeeze(LatticePSFDithered(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(LatticePSFDithered(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(LatticePSFDithered((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

%% Aberrated excitation
Aberrated_LatticePupil = LatticePupil.*exp(1i.*Phase);
[AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticecenter] = SimulateLattice(Aberrated_LatticePupil);
AberratedLatticePSF = AberratedLatticePSF/Latticecenter(1,1);
AberratedLatticePSFDithered = AberratedLatticePSFDithered/Latticecenter(2,1);
%grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,AberratedLatticePSF,AberratedLatticePSFDithered,Aberrated_LatticePupil,AberratedPSFdet,pwd)

xz = squeeze(AberratedLatticePSF(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(AberratedLatticePSF(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(AberratedLatticePSF((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

xz = squeeze(AberratedLatticePSFDithered(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(AberratedLatticePSFDithered(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(AberratedLatticePSFDithered((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

%% Unaberrated Overall
xz = squeeze(AberratedLatticePSFDithered(:,:,(N+1)/2)) .* squeeze(PSFdet(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(AberratedLatticePSFDithered(:,(N+1)/2,:)) .* squeeze(PSFdet(:,(N+1)/2,:)); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(AberratedLatticePSFDithered((N+1)/2,:,:)) .* squeeze(PSFdet((N+1)/2,:,:)); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

%% Aberrated Overall
xz = squeeze(AberratedLatticePSFDithered(:,:,(N+1)/2)) .* squeeze(AberratedPSFdet(:,:,(N+1)/2)); 
xzSthrel = max(xz,[],'all');
xz = xz/xzSthrel;

yz = squeeze(AberratedLatticePSFDithered(:,(N+1)/2,:)) .* flip(squeeze(AberratedPSFdet(:,(N+1)/2,:))',1); 
yzSthrel = max(yz,[],'all');
yz = yz/yzSthrel;

xy = squeeze(AberratedLatticePSFDithered((N+1)/2,:,:)) .* flip(squeeze(AberratedPSFdet((N+1)/2,:,:))',1); 
xySthrel = max(xy,[],'all');
xy = xy/xySthrel;

figure
subplot(1,3,1)
imagesc(X_exc,Y_exc,xy)
title("x,y plane, Strehl=" + num2str(xySthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("y(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)

subplot(1,3,2)
imagesc(Y_exc,Z_exc,yz)
title("y,z plane, Strehl=" + num2str(yzSthrel))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image

colormap(hot)

subplot(1,3,3)
imagesc(X_exc,Z_exc,xz)
title("x,z plane, Strehl=" + num2str(xzSthrel))
xlabel("x(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colorbar
clim([0,1])
axis image
colormap(hot)
