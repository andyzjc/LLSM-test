clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Simulation of diffraction, trying to understand how much energy is lost during diffraction 

% Start with normal pupil, here use a sinc function, do a x-bounding
[SWPupil,SWMask,SWPupilMetaData] = GetSWPairPupil('gaussian',0.3/2,0.3/2,...
                                                           0.3,0.3,...
                                                           1);

xzPSF = abs(fftshift(ifft2(ifftshift(SWPupil(:,:,1))))).^2;

% bound gaussian in x
xFWHM = 50; % pixel 
gaussian_bound = exp(-((-(N-1)/2:(N-1)/2).^2)/ ( (xFWHM).^2));
[xbound,~] = meshgrid(gaussian_bound,gaussian_bound);

% Inverse fourier transform 
xzPSFbound = xzPSF .* xbound; % this is the default energy distribution, i.e unclipped one 

% reverse back to pupil 
boundPupil = fftshift(fft2(ifftshift(xzPSFbound)));
xzPSFbound = abs(fftshift(ifft2(ifftshift(boundPupil)))).^2;

figure(1)
subplot(2,2,1)
imagesc(xzPSF)
axis image
subplot(2,2,2)
imagesc(xzPSFbound)
axis image
subplot(2,2,3)
imagesc(abs(boundPupil))
axis image
subplot(2,2,4)
[IncoherentSumI,absIncoherentSumI,IncoherentIFWHM_unclip,half] = IFWHM(boundPupil);
totalNRG_unclip = sum(sum(xzPSFbound));
peakNRG_unclip = max(max(xzPSFbound));
    plot(Z_exc+0.125,IncoherentSumI,'Color','r','LineWidth',2)
    hold on
    grid on
    plot(Z_exc,absIncoherentSumI,'LineWidth',2)
    axis square
    ylim([0,1])
    xlim([-50,50])
    yline(half)
    title("IFWHM=" + num2str(IncoherentIFWHM_unclip) + "pixel")

%% Real clipping simulation start here:)

% clip the pupil with mask
mask_half_width = [1,2,4,8,16,32]; %pixel 

figure(2)
for ii = 1:length(mask_half_width)
    newpupil = zeros(N,N);
    IncoherentSumI = [];
    absIncoherentSumI = [];

    newpupil(:,(N+1)/2-mask_half_width(1,ii):(N+1)/2+mask_half_width(1,ii)) = boundPupil(:,(N+1)/2-mask_half_width(1,ii):(N+1)/2+mask_half_width(1,ii));
    subplot(3,6,ii)
    imagesc(abs(newpupil))
    axis image
    title(num2str(mask_half_width(1,ii)) + "pixel")

    subplot(3,6,ii+6)
    newPSF = abs(fftshift(ifft2(ifftshift(newpupil)))).^2;
    imagesc(newPSF)
    [IncoherentSumI,absIncoherentSumI,IncoherentIFWHM,half] = IFWHM(newpupil);
    axis image

    subplot(3,6,ii+12)
    plot(Z_exc+0.125,IncoherentSumI,'Color','r','LineWidth',2)
    hold on
    grid on
    plot(Z_exc,absIncoherentSumI,'LineWidth',2)
    axis square
    ylim([0,1])
    yline(half)
    xlim([-50,50])
    title(num2str(IncoherentIFWHM) + "pixel")
end

% just iterate more to get more points
mask_half_width = 1:1:35; %pixel 
IncoherentIFWHM = [];
for ii = 1:length(mask_half_width)
    newpupil = zeros(N,N);
    newpupil(:,(N+1)/2-mask_half_width(1,ii):(N+1)/2+mask_half_width(1,ii)) = boundPupil(:,(N+1)/2-mask_half_width(1,ii):(N+1)/2+mask_half_width(1,ii));

    newPSF = abs(fftshift(ifft2(ifftshift(newpupil)))).^2;
    [~,~,IncoherentIFWHM(1,ii),half] = IFWHM(newpupil);
    % curious about how much energy pass through the mask by looking at the
    % sum of sample
    totalNRG_aftermask(1,ii) = sum(sum(newPSF));
    peakNRG_aftermask(1,ii) = max(max(newPSF));
end

%% Let's do some ratio analysis, to see what is the ratio 
TotalNRG_ratio = totalNRG_unclip./totalNRG_aftermask;
peakNRG_ratio = peakNRG_unclip./peakNRG_aftermask;
IncoherentIFWHM_ratio = IncoherentIFWHM./IncoherentIFWHM_unclip;
figure(3)
subplot(2,2,1)
pupilx_cut = abs(boundPupil((N+1)/2,:))/max(abs(boundPupil((N+1)/2,:)));
index_array = ax(1,:);
plot(index_array,pupilx_cut)
hold on
for ii = 1:length(mask_half_width) % generate block
    block = abs(index_array) <= mask_half_width(1,ii);
    plot(index_array,block)
end
% title("mask half width = "+num2str(mask_half_width) + "pixels")
hold off
xlim([-50,50])
grid on
axis square

% calculate the width of the airy unit
airy_width = max(index_array(pupilx_cut>=0.001))+1;
diffraction_spot_mask_ratio = mask_half_width./airy_width;

subplot(2,2,2)
scatter(diffraction_spot_mask_ratio,TotalNRG_ratio)
title("total NRG")
axis square
grid on
xlabel("width ratio = mask size/diffraction spot")
ylabel("total energy/pass energy")
ylim([0,10])
xlim([0,1])

subplot(2,2,3)
scatter(diffraction_spot_mask_ratio,peakNRG_ratio)
title("peak intensity")
axis square
grid on
xlabel("width ratio = mask size/diffraction spot")
ylabel("peak intensity ratio=unclip/clip")
ylim([0,60])
xlim([0,1])

subplot(2,2,4)
scatter(diffraction_spot_mask_ratio,IncoherentIFWHM_ratio)
title("I_x FWHM ratio")
axis square
grid on
xlabel("width ratio = mask size/diffraction spot")
ylabel("Ix_FWHM (pixel)")
% ylim([0,5])
xlim([0,1])

