%% Simulation without aberrration
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Incoherent sum, sum intensity 
% generate pupil 
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('tophat',0.58,0.58/2,...
                                                           0.04,0.04*2,...
                                                           1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[PSFCoherent_I,~,~] = SimulateSWPair(SWPupil);
PSFCoherent_I = PSFCoherent_I/max(PSFCoherent_I,[],'all');

PSFCoherent = PSFCoherent_I;

yzPSFCoherent = squeeze(PSFCoherent(:,(N+1)/2,:));
yPSFCoherent = yzPSFCoherent((N+1)/2,:);
plot(Y_exc,yPSFCoherent)
hold on
% findpeaks(yPSFCoherent,'MinPeakHeight',0.3);
[~,locs] = findpeaks(yPSFCoherent,'MinPeakHeight',0.3);
shift_period_pixel = round(mean(diff(locs))/2);
ShiftPSF_Incoherent_sum = PSFCoherent + circshift(PSFCoherent,shift_period_pixel,3);
plot(Y_exc,yPSFCoherent+circshift(yPSFCoherent,shift_period_pixel))

FocalPlanexzPSF = PSFCoherent(:,:,(N+1)/2);
FocalPlanexzOTF = fftshift(ifft2(ifftshift(FocalPlanexzPSF)));
FocalPlanexzOTF = FocalPlanexzOTF/max(FocalPlanexzOTF,[],'all');
FocalPlanezOTF = FocalPlanexzOTF(:,(N+1)/2);

% Shifted_FocalPlanexzPSF1 = PSFCoherent(:,:,(N+1)/2+shift_period_pixel);
% Shifted_FocalPlanexzOTF = fftshift(ifft2(ifftshift(Shifted_FocalPlanexzPSF1)));
% Shifted_FocalPlanexzOTF = Shifted_FocalPlanexzOTF/max(Shifted_FocalPlanexzOTF,[],'all');
% Shifted_FocalPlanezOTF = Shifted_FocalPlanexzOTF(:,(N+1)/2);

Sum_focalPlane_xzPSF = squeeze(ShiftPSF_Incoherent_sum(:,:,(N+1)/2+round(shift_period_pixel/2)));
ShiftxzOTF_Incoherent_sum = fftshift(ifft2(ifftshift(Sum_focalPlane_xzPSF)));
ShiftxzOTF_Incoherent_sum = ShiftxzOTF_Incoherent_sum/max(ShiftxzOTF_Incoherent_sum,[],'all');
ShiftzOTF_Incoherent_sum = ShiftxzOTF_Incoherent_sum(:,(N+1)/2);

figure
imagesc(Y_exc,Z_exc,squeeze(ShiftPSF_Incoherent_sum(:,(N+1)/2,:)))
ylim([-40,40]);
xlim([-40,40]);

figure
imagesc(X_exc,Z_exc,Sum_focalPlane_xzPSF) %xzPSF

fig = figure;
plot(KZ_exc,abs(FocalPlanezOTF))
hold on
% plot(KZ_exc,abs(Shifted_FocalPlanezOTF))
plot(KZ_exc,abs(ShiftzOTF_Incoherent_sum))
title("abs(OTF)")
legend("Unshifted beam","new focal")
grid on

fig = figure;
plot(KZ_exc,real(FocalPlanezOTF))
hold on
% plot(KZ_exc,real(Shifted_FocalPlanezOTF))
plot(KZ_exc,real(ShiftzOTF_Incoherent_sum))
title("real(OTF)")
legend("Unshifted beam","new focal")
grid on

fig = figure;
plot(KZ_exc,angle(FocalPlanezOTF))
hold on
% plot(KZ_exc,angle(Shifted_FocalPlanezOTF))
plot(KZ_exc,angle(ShiftzOTF_Incoherent_sum))
ylim([-pi,pi])
title("angle(OTF)")
legend("Unshifted beam","new focal")
grid on

%% Incoherent sum, sum e field

PSFCoherent_E = zeros(N,N,N);
Pupil1 = squeeze(SWPupil(:,:,1));
Pupil2 = squeeze(SWPupil(:,:,2));
Pupil_sum = Pupil1 + Pupil2;
% propagation
for i = 1:length(y_exc)
    propagator_exc = exp(2*pi * 1i * ky_exc * y_exc(i));
    PSFCoherent_E(:,:,i) = fftshift( ifft2(ifftshift(Pupil_sum .* propagator_exc)) );
end
PSFCoherent_E = PSFCoherent_E/max(PSFCoherent_E,[],'all');

PSFCoherent = PSFCoherent_E;

yzPSFCoherent = squeeze(PSFCoherent(:,(N+1)/2,:));
yPSFCoherent = abs(yzPSFCoherent((N+1)/2,:));
figure
plot(Y_exc,yPSFCoherent)
title("abs(E field)")
findpeaks(yPSFCoherent,'MinPeakHeight',0.3);
[~,locs] = findpeaks(yPSFCoherent,'MinPeakHeight',0.3);
shift_period_pixel = round(mean(diff(locs))/2);
ShiftPSF_Coherent_sum = PSFCoherent + circshift(PSFCoherent,shift_period_pixel,3);
ShiftPSF_Coherent_sum = abs(ShiftPSF_Coherent_sum).^2;

FocalPlanexzPSF = PSFCoherent(:,:,(N+1)/2);
FocalPlanexzPSF = abs(FocalPlanexzPSF).^2;

FocalPlanexzOTF = fftshift(ifft2(ifftshift(FocalPlanexzPSF)));
FocalPlanexzOTF = FocalPlanexzOTF/max(FocalPlanexzOTF,[],'all');
FocalPlanezOTF = FocalPlanexzOTF(:,(N+1)/2);

Shifted_FocalPlanexzPSF = PSFCoherent(:,:,(N+1)/2-shift_period_pixel);
Shifted_FocalPlanexzPSF = abs(Shifted_FocalPlanexzPSF).^2;

Shifted_FocalPlanexzOTF = fftshift(ifft2(ifftshift(Shifted_FocalPlanexzPSF)));
Shifted_FocalPlanexzOTF = Shifted_FocalPlanexzOTF/max(Shifted_FocalPlanexzOTF,[],'all');
Shifted_FocalPlanezOTF = Shifted_FocalPlanexzOTF(:,(N+1)/2);

Sum_focalPlane_xzPSF = FocalPlanexzPSF + Shifted_FocalPlanexzPSF;
Sum_focalPlane_xzPSF = abs(Sum_focalPlane_xzPSF).^2;

ShiftxzOTF_Incoherent_sum = fftshift(ifft2(ifftshift(Sum_focalPlane_xzPSF)));
ShiftxzOTF_Incoherent_sum = ShiftxzOTF_Incoherent_sum/max(ShiftxzOTF_Incoherent_sum,[],'all');
ShiftzOTF_Incoherent_sum = ShiftxzOTF_Incoherent_sum(:,(N+1)/2);

figure
imagesc(Y_exc,Z_exc,squeeze(ShiftPSF_Coherent_sum(:,(N+1)/2,:)))
ylim([-40,40]);
xlim([-40,40]);

figure
imagesc(X_exc,Z_exc,Sum_focalPlane_xzPSF) %xzPSF

fig = figure;
plot(KZ_exc,abs(FocalPlanezOTF))
hold on
plot(KZ_exc,abs(Shifted_FocalPlanezOTF))
plot(KZ_exc,abs(ShiftzOTF_Incoherent_sum))
title("abs(OTF)")
legend("Unshifted beam","shifted Beam","Coherent Sum")
grid on

fig = figure;
plot(KZ_exc,real(FocalPlanezOTF))
hold on
plot(KZ_exc,real(Shifted_FocalPlanezOTF))
plot(KZ_exc,real(ShiftzOTF_Incoherent_sum))
title("real(OTF)")
legend("Unshifted beam","shifted Beam","Coherent Sum")
grid on

fig = figure;
plot(KZ_exc,angle(FocalPlanezOTF))
hold on
plot(KZ_exc,angle(Shifted_FocalPlanezOTF))
% plot(KZ_exc,angle(ShiftzOTF_Incoherent_sum))
ylim([-pi,pi])
title("angle(OTF)")
legend("Unshifted beam","shifted Beam","Coherent Sum")
grid on