%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% detection PSF
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

%% Phase

PupilNA = 1.1;
MinRadialOrder = 2;
MaxRadialOrder = 6;
RadioOrderArray = [];
AngularFrequencyArray =[];
counter =1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(1,counter) = i;
        AngularFrequencyArray(1,counter) = AngularFrequency(k);
        counter = counter +1;
    end
end

ModeWeighting = [-0.169329073
-0.006389776
-0.255591054
-0.297124601
0.329073482
-0.121405751
-0.019169329
-0.051118211
0.07028754
0.156549521
1
0.680511182
-0.111821086
0.10543131
0.220447284
0.162939297
0.146964856
-0.14057508
0.063897764
0.044728435
0.134185304
-0.204472843
-0.571884984
-0.504792332
0.54313099
];

tempPhase = zeros(N,N);
Phase = zeros(N,N);
for i = 1:length(ModeWeighting)
    tempPhase = GetSingleZmodePupil(RadioOrderArray(i),AngularFrequencyArray(i),PupilNA);
    Phase = Phase + ModeWeighting(i).*tempPhase;
end
Phase = Phase .* 0.254;
ComplexPhase = exp(1i .* Phase);

fig1 = figure;
    imagesc(KZ_exc,KX_exc,Phase/(2*pi))
    title("Unit: lambda/n")
    colorbar
    colormap(turbo)
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    axis image


%% Pupil
clc
NA1 = 0.58;
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 7/10; %4/3 for equal OTF V2 LLS, 7/10 for V1 LLS
Latticeweighting = 1; % 1.9 for V2 LLS
SNR = 10;
Iter = 10; 
OTFthreshold = 0.001;
wavefrontNA = 1.1;

savingdir = ['Eric_Betzig_Zebrafish_Aberration_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '_SWWeighting_' num2str(SWweighting) '_LLSWeighting_' num2str(Latticeweighting) '_SNR_' num2str(SNR) '/'];
mkdir(savingdir)
print(fig1, '-dsvg', [  savingdir 'Wavefront' '.SVG'],'-r300')
print(fig1, '-dpng', [  savingdir 'Wavefront' '.PNG'],'-r300')   

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

    LatticePupil = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p6.mat');
    temp = LatticePupil.AiryPupil;
    LatticePupil = zeros(size(temp));
    LatticePupil(:,(N+1)/2) = temp(:,(N+1)/2);
elseif isequal(LatticeType,'1dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p6.mat');
    LatticePupil = LatticePupil.GaussianPupil;
elseif isequal(LatticeType,'2dgaussian')
    [SWPupil,~,~] = GetSWPairPupil('gaussian',0,0,...
    0,NA1,...
    SWweighting);

    LatticePupil = load('LLSM-test/Aberration/Zernike/2dGaussianPupil513_NA_0p6.mat');
    LatticePupil = LatticePupil.GaussianPupil;
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
%% Correction 
[~,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');

Pupil1 = SWPupil(:,:,1) .* ComplexPhase;
Pupil2 = SWPupil(:,:,2) .* ComplexPhase;

[~,AberratedBeam1,Beam1Center] = SimulateLattice(Pupil1);
AberratedBeam1 = AberratedBeam1./SWcenter(2,1);
yzAberratedBeam1 = squeeze(AberratedBeam1(:,(N+1)/2,:));
[~,column1] = max(max(yzAberratedBeam1));
[~,row1] = max(yzAberratedBeam1(:,column1));

[~,AberratedBeam2,Beam2Center] = SimulateLattice(Pupil2);
AberratedBeam2 = AberratedBeam2./SWcenter(2,1);
yzAberratedBeam2 = squeeze(AberratedBeam2(:,(N+1)/2,:));
[~,column2] = max(max(yzAberratedBeam2));
[~,row2] = max(yzAberratedBeam2(:,column2));

AberratedPSFIncoherent = AberratedBeam1+AberratedBeam2;
OverallAberratedPSFIncoherent = AberratedPSFIncoherent.*PSFdet;

% y translation, shift NA1 beam to NA2 beam
yshift = column2 - column1;
yCorrectedBeam1 = circshift(AberratedBeam1,yshift,3);
yCorrectedAberratedPSFIncoherent = yCorrectedBeam1 + AberratedBeam2;

% z translation, shift NA1 beam to NA2 beam
zshift = row2 - row1;
yzCorrectedBeam1 = circshift(yCorrectedBeam1,zshift,1);
yzCorrectedAberratedPSFIncoherent = yzCorrectedBeam1 + AberratedBeam2;

%%
xzviewSize = 20;
yViewSize = 100;

fig1 = figure;
imagesc(Y_exc,Z_exc,squeeze(AberratedBeam1(:,(N+1)/2,:)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
ylim([-xzviewSize,xzviewSize])
xlim([-yViewSize,yViewSize])
title("Aberrated yz, NA1")
    print(fig1, '-dsvg', [ savingdir 'AberratedNA1' '.SVG'],'-r300')
    print(fig1, '-dpng', [ savingdir 'AberratedNA1' '.PNG'],'-r300')

fig1 = figure;
imagesc(Y_exc,Z_exc,squeeze(AberratedBeam2(:,(N+1)/2,:)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
ylim([-xzviewSize,xzviewSize])
xlim([-yViewSize,yViewSize])
title("Aberrated yz, NA2")
    print(fig1, '-dsvg', [ savingdir 'AberratedNA2' '.SVG'],'-r300')
    print(fig1, '-dpng', [ savingdir 'AberratedNA2' '.PNG'],'-r300')

fig1 = figure;
imagesc(Y_exc,Z_exc,squeeze(AberratedPSFIncoherent(:,(N+1)/2,:)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
ylim([-xzviewSize,xzviewSize])
xlim([-yViewSize,yViewSize])
title("Aberrated yz")
    print(fig1, '-dsvg', [ savingdir 'Aberratedyz' '.SVG'],'-r300')
    print(fig1, '-dpng', [ savingdir 'Aberratedyz' '.PNG'],'-r300')

fig2 = figure;
imagesc(Y_exc,Z_exc,squeeze(AberratedPSFIncoherent(:,:,(N+1)/2)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
xlim([-xzviewSize,xzviewSize])
ylim([-xzviewSize,xzviewSize])
title("Aberrated xz")
    print(fig2, '-dsvg', [ savingdir 'Aberratedxz' '.SVG'],'-r300')
    print(fig2, '-dpng', [ savingdir 'Aberratedxz' '.PNG'],'-r300')

fig2 = figure;
imagesc(Y_exc,Z_exc,squeeze(yCorrectedAberratedPSFIncoherent(:,(N+1)/2,:)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
ylim([-xzviewSize,xzviewSize])
xlim([-yViewSize,yViewSize])
title("y-correct yz")
    print(fig2, '-dsvg', [ savingdir 'y_corrected_yz' '.SVG'],'-r300')
    print(fig2, '-dpng', [ savingdir 'y_corrected_yz' '.PNG'],'-r300')

fig2 = figure;
imagesc(Y_exc,Z_exc,squeeze(yCorrectedAberratedPSFIncoherent(:,:,column2)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
xlim([-xzviewSize,xzviewSize])
ylim([-xzviewSize,xzviewSize])
title("y-correct xz, at NA2 max (Not the same focal)")
    print(fig2, '-dsvg', [ savingdir 'y_corrected_xz' '.SVG'],'-r300')
    print(fig2, '-dpng', [ savingdir 'y_corrected_xz' '.PNG'],'-r300')

fig2 = figure;
imagesc(Y_exc,Z_exc,squeeze(yzCorrectedAberratedPSFIncoherent(:,(N+1)/2,:)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
ylim([-xzviewSize,xzviewSize])
xlim([-yViewSize,yViewSize])
title("yz-correct yz")
    print(fig2, '-dsvg', [ savingdir 'yz_corrected_yz' '.SVG'],'-r300')
    print(fig2, '-dpng', [ savingdir 'yz_corrected_yz' '.PNG'],'-r300')

fig2 = figure;
imagesc(Y_exc,Z_exc,squeeze(yzCorrectedAberratedPSFIncoherent(:,:,column2)))
xlabel("y(\lambda_{exc}/n)");
ylabel("z(\lambda_{exc}/n)");
colormap(fire(256))
axis image
xlim([-xzviewSize,xzviewSize])
ylim([-xzviewSize,xzviewSize])
title("yz-correct xz, at NA2 max (Not the same focal)")
    print(fig2, '-dsvg', [ savingdir 'yz_corrected_xz' '.SVG'],'-r300')
    print(fig2, '-dpng', [ savingdir 'yz_corrected_xz' '.PNG'],'-r300')

close all

%% PSF
% [~,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
% PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
%     yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
%     SWyFWHM1 = find(yindex,1,'first') ;
%     SWyFWHM2 = find(yindex,1,'last');
%     overallOTF = fftshift(ifftn(ifftshift(PSFIncoherent.* PSFdet)));
%     xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
%     SWOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;
% 
% [~,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
% LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
%     yindex = 1-(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:))));
%     LatticeyFWHM1 = find(yindex,1,'first');
%     LatticeyFWHM2 = find(yindex,1,'last');
%     overallOTF = fftshift(ifftn(ifftshift(LatticePSFDithered.* PSFdet)));
%     xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
%     LatticeOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;
% 
% AberratedSWPupil = SWPupil.*ComplexPhase;
% AberratedLatticePupil= LatticePupil.*ComplexPhase;
% 
% [AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);
% AberratedPSFIncoherent = AberratedPSFIncoherent./SWcenter(2,1);
% AberratedPSFCoherent = AberratedPSFCoherent./SWcenter(1,1);
% [AberratedLatticePSF,AberratedLatticePSFDithered,~] = SimulateLattice(AberratedLatticePupil);
% AberratedLatticePSFDithered = AberratedLatticePSFDithered./Latticecenter(2,1);
% AberratedLatticePSF = AberratedLatticePSF/Latticecenter(1,1);
% 
% OverallAberratedPSFIncoherent = AberratedPSFIncoherent.* PSFdet;
% OverallAberratedLatticePSFDithered = AberratedLatticePSFDithered.* PSFdet;
% 
% %% Analysis
% SRatioSW = max(OverallAberratedPSFIncoherent(:,:,(N+1)/2),[],'all');
% SRatioSW_corrected = (AberratedSWcenter(3,1) + AberratedSWcenter(4,1))/SWcenter(2,1);
% SRatioLattice = max(OverallAberratedLatticePSFDithered(:,:,(N+1)/2),[],'all');
% 
% [AberratedSWAveragefc3,AberratedSWAveragefc3FWHM,~,~] = RunFC3(AberratedPSFIncoherent,PSFdet,SWOTFmask,Iter,SNR,SWyFWHM2);
% [AberratedLatticeAveragefc3,AberratedLatticeAveragefc3FWHM,~,~] = RunFC3(AberratedLatticePSFDithered,PSFdet,LatticeOTFmask,Iter,SNR,LatticeyFWHM2);
% 
% [SWconvLines,~,~,~] = ConvRes(PSFIncoherent,PSFdet,SNR);
% [LatticeconvLines,~,~,~] = ConvRes(LatticePSFDithered,PSFdet,SNR);
% [AberratedSWconvLines,lineSpot,Spacing,LineZ] = ConvRes(AberratedPSFIncoherent,PSFdet,SNR);
% [AberratedLatticeconvLines,~,~,~] = ConvRes(AberratedLatticePSFDithered,PSFdet,SNR);
% 
% GTBeads = load('/Users/andyzjc/Dropbox (Princeton)/Tian-Ming for Andy/Manuscripts/Adaptive Polarization Controlled Incoherent SW Light Sheet/Figures/Figure 5/hex_0.58_0.04_tophat_V1/AberrationAnalysis/BeadsImage/GroundTruth.mat');
% Beads = GTBeads.GTallbeads;
% Vol = zeros(N,N,N);
% Vol(:,:,(N+1)/2) = Beads;
% 
% fftvol = fftshift(fftn(ifftshift(Vol))) .* fftshift(fftn(ifftshift(OverallAberratedPSFIncoherent)));
% temp = abs(fftshift(ifftn(ifftshift(fftvol))));
% AberratedSWBeadsVol = temp + poissrnd(temp) .* 1/SNR;
% 
% fftvol = fftshift(fftn(ifftshift(Vol))) .* fftshift(fftn(ifftshift(OverallAberratedLatticePSFDithered)));
% temp = abs(fftshift(ifftn(ifftshift(fftvol))));
% AberratedLatticeBeadsVol = temp + poissrnd(temp) .* 1/SNR;
% 
% %% Plot
% grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticePupil,PSFdet,savingdir)
% grapthSW(NA1,deltaNA,LatticeType,SWweighting,AberratedPSFIncoherent,AberratedPSFCoherent,AberratedSWPupil,PSFdet,savingdir)
% 
% fig1 = figure;
%     imagesc(Z_exc,X_exc,AberratedSWBeadsVol(:,:,(N+1)/2))
%     axis image
%     xlim([-40,40])
%     ylim([-40,40])
%     clim([0,1])
%     colormap(hot)
%     colorbar
%     print(fig1, '-dsvg', [ savingdir 'AberratedSWBeads' '.SVG'],'-r300')
%     print(fig1, '-dpng', [ savingdir 'AberratedSWBeads' '.PNG'],'-r300')
% 
% fig2 = figure;
%     imagesc(Z_exc,X_exc,AberratedLatticeBeadsVol(:,:,(N+1)/2))
%     axis image
%     xlim([-40,40])
%     ylim([-40,40])
%     colormap(hot)
%     clim([0,1])
%     colorbar
%     print(fig2, '-dsvg', [  savingdir 'LatticeBeads' '.SVG'],'-r300')
%     print(fig2, '-dpng', [  savingdir 'LatticeBeads' '.PNG'],'-r300')
% 
% fig3 = figure;
% h1 = subplot(2,2,1,'Parent',fig3);
%     imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3)
%     title("SW")
%     xlabel("k_r/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     % colormap(h1,fire(256))
%     colormap(h1,turbo)
%     colorbar
%     clim([0,1])
%     set(gca,'YDir','normal')
%     axis image
% 
% h1 = subplot(2,2,2,'Parent',fig3);
%     imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3)
%     title("LLS")
%     xlabel("k_r/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     % colormap(h1,fire(256))
%     colormap(h1,turbo)
%     colorbar
%     clim([0,1])
%     set(gca,'YDir','normal')
%     axis image
% 
% h1 = subplot(2,2,3,'Parent',fig3);
%     imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3FWHM)
%     title("SW FWHM")
%     xlabel("k_r/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     % colormap(h1,fire(256))
%     colormap(h1,turbo)
%     colorbar
%     clim([0,1])
%     set(gca,'YDir','normal')
%     axis image
% 
% h1 = subplot(2,2,4,'Parent',fig3);
%     imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3FWHM)
%     title("LLS FWHM")
%     xlabel("k_r/(4\pin/\lambda_{exc})")
%     ylabel("k_z/(4\pin/\lambda_{exc})")
%     % colormap(h1,fire(256))
%     colormap(h1,turbo)
%     colorbar
%     clim([0,1])
%     set(gca,'YDir','normal')
%     axis image
% 
% print(fig3, '-dsvg', [  savingdir 'FC3' '.SVG'],'-r300')
% print(fig3, '-dpng', [  savingdir 'FC3' '.PNG'],'-r300')
% 
% fig4 = figure;
%     plot(LineZ,SWconvLines((N+1)/2,:),'magenta','LineWidth',0.5)
%     hold on
%     plot(LineZ,LatticeconvLines((N+1)/2,:),'Color','b','LineWidth',0.5)
%     plot(LineZ,AberratedSWconvLines((N+1)/2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5)
%     plot(LineZ,AberratedLatticeconvLines((N+1)/2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',0.5)
%     xline(LineZ(lineSpot),'k')
%     xlim([0,50])
%     ylim([0,50])
%     legend("pSW","LLS","aberrated iSW","aberrated LLS")
%     xlabel("z(lambda/n)")
%     ylabel("Intensity (a.u)")
%     hold off 
%     grid on
%     pbaspect([7 2 1])
% print(fig4, '-dsvg', [  savingdir 'LineGrating' '.SVG'],'-r300')
% print(fig4, '-dpng', [  savingdir 'LineGrating' '.PNG'],'-r300')
% 
% fig5 = figure;
%     imagesc(LineZ,LineZ,AberratedSWconvLines)
%     colormap(hot)
%     hold on
%     axis image
%     xline(LineZ(lineSpot),'k','LineWidth',1)
%     xlim([0,50])
%     ylim([0,20])
%     clim([0,50])
%     colorbar
%     xlabel("z(um)")
%     ylabel("um")
%     hold off
% print(fig5, '-dsvg', [  savingdir 'SWLines' '.SVG'],'-r300')
% print(fig5, '-dpng', [  savingdir 'SWLines' '.PNG'],'-r300')    
% 
% fig6 = figure;
%     imagesc(LineZ,LineZ,AberratedLatticeconvLines)
%     colormap(hot)
%     hold on
%     xline(LineZ(lineSpot),'k','LineWidth',1)
%     axis image
%     xlim([0,50])
%     ylim([0,20])
%     clim([0,50])
%     colorbar
%     xlabel("z(um)")
%     ylabel("um")
%     hold off
% print(fig6, '-dsvg', [  savingdir 'LLSLines' '.SVG'],'-r300')
% print(fig6, '-dpng', [  savingdir 'LLSLines' '.PNG'],'-r300')  
% close all