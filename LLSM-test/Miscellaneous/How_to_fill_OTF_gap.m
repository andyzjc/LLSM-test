clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Simulation of SW Pair
hex1NA1 = 0.58;
hex1deltaNA1 = 0.04;
% Phase = GetSingleZmodePupil(2,0,0.65);
% Phase_factor = exp(1i.*Phase);

%% hex1
% [SWPupil,SWMask] = GetSWPairPupil(ProfileType,NA1Ideal,NA2Ideal,
%                                               deltaNA1,deltaNA2,...
%                                               WeightRatio=I(NA1)/I(NA2));
[hex1SWPupil,~,~] = GetSWPairPupil('tophat',hex1NA1,hex1NA1/2,...
                                           hex1deltaNA1,hex1deltaNA1*2,...
                                           1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[hex1PSFCoherent,~,~] = SimulateSWPair(hex1SWPupil);
% [hex1PSFCoherent_defocus,~,~] = SimulateSWPair(hex1SWPupil.*Phase_factor);

hex1PSFCoherent = hex1PSFCoherent./max(hex1PSFCoherent,[],'all');
CoherentHex1yzPSFexc = squeeze(hex1PSFCoherent(:,(N+1)/2,:));
CoherentHex1xzPSFexc = squeeze(hex1PSFCoherent(:,:,(N+1)/2));
CoherentHex1yPSFexc = squeeze(hex1PSFCoherent((N+1)/2,(N+1)/2,:));
CoherentHex1zPSFexc = squeeze(hex1PSFCoherent(:,(N+1)/2,(N+1)/2));

    % yindex = 1-(CoherentHex1yPSFexc(1:(N+1)/2) >= 0.5*max(squeeze(CoherentHex1yPSFexc(1:(N+1)/2))));
    % hexFWHMindex = find(yindex,1,'last');
    % hexyFWHN = abs(Y_exc(hexFWHMindex)*2)
    
% PrettyPlotSWPair(SWPupil,SWMask,SWPupilMetaData,PSFCoherent,PSFIncoherent,SWcenter);

%% hex2
hex2NA1 = 1.8/2 * hex1NA1 ;
hex2deltaNA1 = ((hex1NA1+hex1deltaNA1/2)-(hex1NA1-hex1deltaNA1/2)).*((hex1NA1+hex1deltaNA1/2)+(hex1NA1-hex1deltaNA1/2)) ./ (2*hex2NA1);

factor = 0.6;
% hex2NA2 = 0.75/2 * hex1NA1 ;
hex2NA2 = factor.*1/2 * hex2NA1 ;
hex2deltaNA2 = ((hex1NA1+hex1deltaNA1/2)-(hex1NA1-hex1deltaNA1/2)).*((hex1NA1+hex1deltaNA1/2)+(hex1NA1-hex1deltaNA1/2)) ./ (2*hex2NA2);

[hex2SWPupil,~,~] = GetSWPairPupil('tophat',hex2NA1,hex2NA2,...
                                           hex2deltaNA1,hex2deltaNA2,...
                                           1);

% Coherent/Incoherent Propagation of SW pairs
% [PSFCoherent,PSFIncoherent] = SimulateSWPair(SWPupil);
[hex2PSFCoherent,~,~] = SimulateSWPair(hex2SWPupil);
% [hex2PSFCoherent_defocus,~,~] = SimulateSWPair(hex2SWPupil.*Phase_factor);

hex2PSFCoherent = hex2PSFCoherent./max(hex2PSFCoherent,[],'all');
CoherentHex2yzPSFexc = squeeze(hex2PSFCoherent(:,(N+1)/2,:));
CoherentHex2xzPSFexc = squeeze(hex2PSFCoherent(:,:,(N+1)/2));
CoherentHex2yPSFexc = squeeze(hex2PSFCoherent((N+1)/2,(N+1)/2,:));
CoherentHex2zPSFexc = squeeze(hex2PSFCoherent(:,(N+1)/2,(N+1)/2));

%% adding two PSF with same T, Incoherently
CoherentHex1zOTFexc = fftshift(ifft(ifftshift(CoherentHex1zPSFexc)));CoherentHex1zOTFexc = CoherentHex1zOTFexc./max(CoherentHex1zOTFexc,[],'all');
CoherentHex2zOTFexc = fftshift(ifft(ifftshift(CoherentHex2zPSFexc))); CoherentHex2zOTFexc = CoherentHex2zOTFexc./max(CoherentHex2zOTFexc,[],'all');

shift_distance = 27; % remember to change when change NA/deltaNA
sumyzPSF = CoherentHex1yzPSFexc + circshift(CoherentHex2yzPSFexc,shift_distance,2); sumyzPSF = sumyzPSF/max(sumyzPSF,[],'all');
sumyPSF = sumyzPSF((N+1)/2,:);
sumzPSF = sumyzPSF(:,(N+1)/2);
sumzOTF = fftshift(ifft(ifftshift(sumzPSF))); sumzOTF = sumzOTF/max(sumzOTF,[],'all');

% close all
figure
imagesc(imfuse(hex2SWPupil(:,:,1)+hex2SWPupil(:,:,2), hex1SWPupil(:,:,1)+hex1SWPupil(:,:,2)))

figure
plot(Y_exc,CoherentHex1yPSFexc)
hold on
plot(Y_exc,CoherentHex2yPSFexc)
plot(Y_exc,sumyPSF)
title("yPSF")
legend("hex reference","hex fill","shift,sum")

%% incoherently, add up 4 PSFs, physically not impossible
CoherentHex1zOTFexc = fftshift(ifft(ifftshift(CoherentHex1zPSFexc)));CoherentHex1zOTFexc = CoherentHex1zOTFexc./max(CoherentHex1zOTFexc,[],'all');
CoherentHex2zOTFexc = fftshift(ifft(ifftshift(CoherentHex2zPSFexc))); CoherentHex2zOTFexc = CoherentHex2zOTFexc./max(CoherentHex2zOTFexc,[],'all');

sumyzPSF = CoherentHex1yzPSFexc + CoherentHex2yzPSFexc; sumyzPSF = sumyzPSF/max(sumyzPSF,[],'all');
sumyPSF = sumyzPSF((N+1)/2,:);
sumzPSF = sumyzPSF(:,(N+1)/2);
sumzOTF = fftshift(ifft(ifftshift(sumzPSF))); sumzOTF = sumzOTF/max(sumzOTF,[],'all');

shift_distance = 28; % remember to change when change NA/deltaNA
sumyzPSF_shift = sumyzPSF + circshift(sumyzPSF,shift_distance,2); sumyzPSF_shift = sumyzPSF_shift/max(sumyzPSF_shift,[],'all');
sumyPSF_shift = sumyzPSF_shift((N+1)/2,:);
sumzPSF_shift = sumyzPSF_shift(:,(N+1)/2+round(shift_distance/2)); % there is a new center/focal plane
sumzOTF_shift = fftshift(ifft(ifftshift(sumzPSF_shift))); sumzOTF_shift = sumzOTF_shift/max(sumzOTF_shift,[],'all');

close all
figure
imagesc(imfuse(hex2SWPupil(:,:,1)+hex2SWPupil(:,:,2), hex1SWPupil(:,:,1)+hex1SWPupil(:,:,2)))

figure
plot(Y_exc,CoherentHex1yPSFexc)
hold on
plot(Y_exc,CoherentHex2yPSFexc)
plot(Y_exc,sumyPSF)
plot(Y_exc,sumyPSF_shift)
title("yPSF")
legend("hex reference","hex fill","sum","sum and shift")

figure
plot(Z_exc,CoherentHex1zPSFexc)
hold on
plot(Z_exc,CoherentHex2zPSFexc)
plot(Z_exc,sumzPSF)
plot(Z_exc,sumzPSF_shift)
title("zPSF")
legend("hex reference","hex fill","sum","sum and shift")
xlim([-15,15])

figure
subplot(2,2,1)
imagesc(CoherentHex1yzPSFexc)
axis image
subplot(2,2,2)
imagesc(CoherentHex2yzPSFexc)
axis image
subplot(2,2,3)
imagesc(sumyzPSF)
axis image
subplot(2,2,4)
imagesc(sumyzPSF_shift)
axis image

figure
grid on
plot(KZ_exc,real(CoherentHex1zOTFexc))
hold on
plot(KZ_exc,real(CoherentHex2zOTFexc))
plot(KZ_exc,real(sumzOTF))
plot(KZ_exc,real(sumzOTF_shift))
title("zOTF")
legend("hex reference","hex fill","sum","sum and shift")

figure
title("xzPSF")
subplot(2,2,1)
imagesc(CoherentHex1xzPSFexc)
axis image
subplot(2,2,2)
imagesc(CoherentHex2xzPSFexc)
axis image

%% Coherently add up, option 1
CoherentHex1zOTFexc = fftshift(ifft(ifftshift(abs(CoherentHex1zPSFexc).^2)));CoherentHex1zOTFexc = CoherentHex1zOTFexc./max(CoherentHex1zOTFexc,[],'all');
CoherentHex2zOTFexc = fftshift(ifft(ifftshift(abs(CoherentHex2zPSFexc).^2))); CoherentHex2zOTFexc = CoherentHex2zOTFexc./max(CoherentHex2zOTFexc,[],'all');

hex1_shift_distance = 0;
hex1_shift_and_sum_yz_efield = CoherentHex1yzPSFexc + circshift(CoherentHex1yzPSFexc,hex1_shift_distance,2);
hex1_shift_and_sum_yz_PSF = abs(hex1_shift_and_sum_yz_efield).^2; hex1_shift_and_sum_yz_PSF = hex1_shift_and_sum_yz_PSF/max(hex1_shift_and_sum_yz_PSF,[],'all');
hex1_shift_and_sum_y_PSF = hex1_shift_and_sum_yz_PSF((N+1)/2,:);
hex1_shift_and_sum_z_PSF = hex1_shift_and_sum_yz_PSF(:,(N+1)/2);
hex1_shift_and_sum_z_OTF = fftshift(ifft2(ifftshift(hex1_shift_and_sum_yz_PSF))); 

hex2_shift_distance = 32;
hex2_shift_and_sum_yz_efield = CoherentHex2yzPSFexc + circshift(CoherentHex2yzPSFexc,hex2_shift_distance,2);
hex2_shift_and_sum_yz_PSF = abs(hex2_shift_and_sum_yz_efield).^2; hex2_shift_and_sum_yz_PSF = hex2_shift_and_sum_yz_PSF/max(hex2_shift_and_sum_yz_PSF,[],'all');
hex2_shift_and_sum_y_PSF = hex2_shift_and_sum_yz_PSF((N+1)/2,:);
hex2_shift_and_sum_z_PSF = hex2_shift_and_sum_yz_PSF(:,(N+1)/2);
hex2_shift_and_sum_z_OTF = 0;

sumyzPSF =  hex1_shift_and_sum_yz_PSF+ hex2_shift_and_sum_yz_PSF; 
sumyPSF = sumyzPSF((N+1)/2,:);
sumzPSF = sumyzPSF(:,(N+1)/2);
sumzOTF = fftshift(ifft(ifftshift(sumzPSF))); sumzOTF = sumzOTF/max(sumzOTF,[],'all');

close all
figure
imagesc(imfuse(hex2SWPupil(:,:,1)+hex2SWPupil(:,:,2), hex1SWPupil(:,:,1)+hex1SWPupil(:,:,2)))

figure
title("hex1")
plot(Y_exc,abs(CoherentHex1yPSFexc).^2)
hold on
plot(Y_exc,circshift(abs(CoherentHex1yzPSFexc((N+1)/2,:)).^2,hex1_shift_distance))
plot(Y_exc,hex1_shift_and_sum_y_PSF)
legend("hex reference E,I","hex reference shifted E, I","I,E field sum")

figure
title("hex1")
plot(Y_exc,abs(CoherentHex2yPSFexc).^2)
hold on
plot(Y_exc,circshift(abs(CoherentHex2yzPSFexc((N+1)/2,:)).^2,hex2_shift_distance))
plot(Y_exc,hex2_shift_and_sum_y_PSF)
legend("hex fill E,I","hex fill shifted E, I","I,E field sum")

figure
plot(Y_exc,hex1_shift_and_sum_y_PSF)
hold on
plot(Y_exc,hex2_shift_and_sum_y_PSF)
plot(Y_exc,sumyPSF)
title("yPSF")
legend("hex reference coherent sum","hex fill coherent sum","Incoherent sum")

figure
plot(Z_exc,hex1_shift_and_sum_z_PSF)
hold on
plot(Z_exc,hex2_shift_and_sum_z_PSF)
plot(Z_exc,sumzPSF)
title("zPSF")
legend("hex reference coherent sum","hex fill coherent sum","Incoherent sum")
xlim([-15,15])

figure
subplot(2,2,1)
imagesc(CoherentHex1yzPSFexc)
axis image
subplot(2,2,2)
imagesc(CoherentHex2yzPSFexc)
axis image
subplot(2,2,3)
imagesc(sumyzPSF)
axis image
subplot(2,2,4)
imagesc(sumyzPSF_shift)
axis image

figure
grid on
plot(KZ_exc,real(CoherentHex1zOTFexc))
hold on
plot(KZ_exc,real(CoherentHex2zOTFexc))
plot(KZ_exc,real(sumzOTF))
plot(KZ_exc,real(sumzOTF_shift))
title("zOTF")
legend("hex reference","hex fill","sum","sum and shift")

figure
title("xzPSF")
subplot(2,2,1)
imagesc(CoherentHex1xzPSFexc)
axis image
subplot(2,2,2)
imagesc(CoherentHex2xzPSFexc)
axis image



