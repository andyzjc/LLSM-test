%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% detection 
% detection
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));
xzPSFdet = PSFdet(:,:,(N+1)/2);

%% Unaberrated 
clc
NA1 = 0.58; %0.6 for gaussian
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 4/3; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1.9,1,1,1.9,1]; % 1.9 for V2 LLS
PupilNA = 0.65;

SWPupil = zeros(N,N);
LatticePupil = zeros(N,N);
k_apertureNA = NA1 * k_wave / n;
if isequal(LatticeType,'hex')
    [SWPupil,SWmask,SWPupilMetaData] = GetSWPairPupil(ProfileType,NA1,NA1/2,...
    deltaNA,2*deltaNA,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'SW')
    [SWPupil,SWmask,SWPupilMetaData] = GetSWPairPupil(ProfileType,NA1,NA1,...
    deltaNA,deltaNA,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil('hex',ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'NC')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,NA1*3/4,...
    deltaNA,deltaNA*4/3,...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil('hex',ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
elseif isequal(LatticeType,'sinc')
    [SWPupil,~,~] = GetSWPairPupil('tophat',0,0,...
    NA1*2,NA1*2,...
    SWweighting);

    LatticePupil(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));

    % LatticePupil = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p6.mat');
    % temp = LatticePupil.AiryPupil;
    % LatticePupil = zeros(size(temp));
    % LatticePupil(:,(N+1)/2) = temp(:,(N+1)/2);
elseif isequal(LatticeType,'1dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);
    
    gaussian_mask = zeros(N,N);
    gaussian_mask(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
    LatticePupil(:,(N+1)/2) = exp( -(kz_exc(:,1).^2)/ ((k_apertureNA).^2) );
    % LatticePupil = LatticePupil.* gaussian_mask;

    % LatticePupil = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p6.mat');
    % LatticePupil = LatticePupil.GaussianPupil;
elseif isequal(LatticeType,'2dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);

    gaussian_mask = zeros(N,N);
    gaussian_mask= (k_apertureNA).^2 >= abs((kx_exc.^2 + kz_exc.^2));
    LatticePupil = exp( -(kx_exc.^2 + kz_exc.^2)/ ((k_apertureNA).^2) );
    LatticePupil = LatticePupil.* gaussian_mask;
    % LatticePupil = load('LLSM-test/Aberration/Zernike/2dGaussianPupil513_NA_0p6.mat');
    % LatticePupil = LatticePupil.GaussianPupil;
elseif isequal(LatticeType,'bessel')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',NA1,NA1,...
    deltaNA,deltaNA,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/1dBesselPupil513_NA_0p6_0p05.mat');
    LatticePupil = LatticePupil.BesselPupil;
else % square
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,sqrt(2*NA1*deltaNA)/2,...
    deltaNA,sqrt(2*NA1*deltaNA),...
    SWweighting);

    [LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA1,deltaNA, ...
    0.8,0,...
    Latticeweighting);
end
%% SW
[~,PSFIncoherent,~] = SimulateSWPair(SWPupil);
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');

%% Lattice & and other
[~,LatticePSFDithered,~] = SimulateLattice(LatticePupil);
LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');

%% Making video
outputVideo = VideoWriter('outputVideo.mp4', 'MPEG-4'); % Create a video writer object
outputVideo.FrameRate = 20; % Set the frame rate (frames per second)
open(outputVideo); % Open the video writer object
firemap = fire(256);

yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
SWyFWHM1 = find(yindex,1,'first') ;
SWyFWHM2 = find(yindex,1,'last');

fig = figure;
fig.WindowState = 'maximized';
hf=colordef(fig,'black'); %Set color scheme
hf.Color='k'; %Set background color of figure window

for jj = 1:N
PSF1 = PSFIncoherent(:,:,jj);
overallPSF1 = PSF1 .* xzPSFdet;

PSF2 = LatticePSFDithered(:,:,jj);
overallPSF2 = PSF2 .* xzPSFdet;

s1 = subplot(3,3,4);
imagesc(X_exc,Z_exc,PSF1);
title("PEARLS excitation xzPSF")
axis image;
xlim([-20,20])
ylim([-20,20])
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
colormap(firemap)
clim([0,1])
s1Pos = get(s1,'position');
bar = colorbar;
bar.Location = 'eastoutside';
set(s1,'position',s1Pos)
s1.XAxis.FontSize = 15;
s1.YAxis.FontSize = 15;
s1.Title.FontSize = 15;

s2 = subplot(3,3,5);
imagesc(X_exc,Z_exc,overallPSF1);
title("PEARLS Overall PSF")
axis image;
xlim([-5,5])
ylim([-5,5])
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
colormap(firemap)
clim([0,1])
s2.XAxis.FontSize = 15;
s2.YAxis.FontSize = 15;
s2.Title.FontSize = 15;

s3 = subplot(3,3,6);
plot(Z_exc,PSF1(:,(N+1)/2),'Color','magenta','LineWidth',2);
hold on
plot(Z_exc,overallPSF1(:,(N+1)/2),'Color','g','LineWidth',2);
title("PEARLS Axial linecut, X=0")
axis square
xlim([-10,10])
ylim([0,1])
xlabel("z/(\lambda_{exc}/n)")
ylabel("Intensity (a.u.)")
lgd = legend("Excitation","Overall");
lgd.FontSize = 15;
lgd.Box = 'off';
lgd.Orientation = 'vertical';
lgd.Location = 'north';
grid on
hold off
s3.XAxis.FontSize = 15;
s3.YAxis.FontSize = 15;
s3.Title.FontSize = 15;

s4 = subplot(3,3,1);
imagesc(X_exc,Z_exc,PSF2);
title("LLS excitation xzPSF")
axis image;
xlim([-20,20])
ylim([-20,20])
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
colormap(firemap)
clim([0,1])
s4.XAxis.FontSize = 15;
s4.YAxis.FontSize = 15;
s4.Title.FontSize = 15;

s5 = subplot(3,3,2);
imagesc(X_exc,Z_exc,overallPSF2);
title("LLS Overall PSF")
axis image;
xlim([-5,5])
ylim([-5,5])
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
colormap(firemap)
clim([0,1])
s5.XAxis.FontSize = 15;
s5.YAxis.FontSize = 15;
s5.Title.FontSize = 15;

s6 = subplot(3,3,3);
plot(Z_exc,PSF2(:,(N+1)/2),'Color','magenta','LineWidth',2);
hold on
plot(Z_exc,overallPSF2(:,(N+1)/2),'Color','g','LineWidth',2);
axis square
title("LLS Axial linecut, X=0")
xlim([-10,10])
ylim([0,1])
xlabel("z/(\lambda_{exc}/n)")
ylabel("Intensity (a.u.)")
lgd = legend("Excitation","Overall");
lgd.FontSize = 15;
lgd.Box = 'off';
lgd.Orientation = 'vertical';
lgd.Location = 'North';
grid on
hold off
s6.XAxis.FontSize = 15;
s6.YAxis.FontSize = 15;
s6.Title.FontSize = 15;

s7 = subplot(3,3,7:9);
imagesc(Y_exc,Z_exc,squeeze(PSFIncoherent(:,(N+1)/2,:)));
hold on
xline(Y_exc(jj),'LineWidth',2,'Color','g')
x1 = xline(Y_exc(SWyFWHM1),'Color','cyan','LineWidth',2);
x2 = xline(Y_exc(SWyFWHM2),'Color','cyan','LineWidth',2);
title("X=0, YZ excitation PSF, " + "Y=" + num2str(Y_exc(jj),'%.2f') + "\lambda_{exc}/n")
axis image;
xlim([-107,107])
ylim([-20,20])
xlabel("y/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
lgd = legend(x1,'Non-diffractive region');
lgd.FontSize = 15;
hold off
s7.XAxis.FontSize = 15;
s7.YAxis.FontSize = 15;
s7.Title.FontSize = 15;

% Capture the frame
frame = getframe(gcf);
writeVideo(outputVideo, frame); % Write the frame to the video
% close(gcf); % Close the figure

% exportgraphics(fig4,'Profile.gif','Append',true,'Resolution',500)
% clear fig4
% delete fig4

end

% Close the video file
close(outputVideo);
disp('Video creation complete.');