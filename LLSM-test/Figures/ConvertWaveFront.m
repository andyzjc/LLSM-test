clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

getParameters; %modify image parameter here
CalculatePhysics;

%% detection PSF
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

%% mask
image = Tiff('/Users/andyzjc/Dropbox (Princeton)/Tian-Ming for Andy/Manuscripts/Adaptive Polarization Controlled Incoherent SW Light Sheet/Figures/Figure 5/Wavefront.tif');
wavefront = read(image);
% wavefront(wavefront>=220) = 0;
wavefront = double(wavefront);
[row,col] = size(wavefront);
[x,y] = meshgrid(-(row-1)/2:(row-1)/2,-(col-1)/2:(col-1)/2);
mask = (row/2-3).^2 >= (x.^2+y.^2); 
wavefront = wavefront.*mask;

%% scale to 0.65 NA
PupilNA = 0.65;
[theta,r] = cart2pol(kx_exc./(PupilNA./n*k_wave),kz_exc./(PupilNA./n*k_wave));
idx = r<=1;
k_NA = (PupilNA * k_wave/n);
deltak_wavefront = k_NA/(row/2-5);

wavefront_pad = padarray(wavefront,[round(k_bound/deltak_wavefront-(row/2-5)),round(k_bound/deltak_wavefront-(row/2-5))],0);

%% scale to N
rescale_factor = deltak_wavefront/deltak;
up_Image_size = round(size(wavefront_pad,1) * rescale_factor);
Image_center = (up_Image_size+1)/2;
Scaledwavefront = imresize(wavefront_pad,[up_Image_size,up_Image_size]);
Scaledwavefront = Scaledwavefront(Image_center-(N+1)/2+1:Image_center+(N+1)/2-1,...
                          Image_center-(N+1)/2+1:Image_center+(N+1)/2-1);
%% scale intensity
maxWFE = 2; %um
IntScaledwavefront = Scaledwavefront/max(Scaledwavefront,[],'all');
IntScaledwavefront = maxWFE.*(IntScaledwavefront');
IntScaledwavefront(idx) = IntScaledwavefront(idx) - mode(mode(IntScaledwavefront(idx)));
phase = IntScaledwavefront * 2*pi/(wavelength_exc);
ComplexPhase = exp(1i.*phase);

fig1 = figure;
    imagesc(KZ_exc,KX_exc,IntScaledwavefront)
    title("Unit: um")
    colorbar
    colormap(turbo)
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    axis image
    clim([-0.5,0.5])

fig2 = figure;
    imagesc(KZ_exc,KX_exc,IntScaledwavefront/wavelength_exc)
    title("Unit: lambda/n")
    colorbar
    colormap(turbo)
    xlabel("k_x/(4\pin/\lambda_{exc})");
    ylabel("k_z/(4\pin/\lambda_{exc})");
    axis image
    % clim([min(IntScaledwavefront/wavelength_exc,[],'all'),max(IntScaledwavefront/wavelength_exc,[],'all')])
    clim([-2,2])
%% Pupil
clc
NA1 = 0.58;
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 7/10; %4/3 for equal OTF V2 LLS, 7/10 for V1 LLS
Latticeweighting = [1 1 1 1 1 1]; % 1.9 for V2 LLS
SNR = 10;
Iter = 10;
OTFthreshold = 0.001;

savingdir = ['NaJi_Zebrafish_motor_neuron_Aberration_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '_SWWeighting_' num2str(SWweighting) '_LLSWeighting_' num2str(Latticeweighting) '_SNR_' num2str(SNR) '/'];
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




%% PSF
[~,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
    yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
    SWyFWHM1 = find(yindex,1,'first') ;
    SWyFWHM2 = find(yindex,1,'last');
    overallOTF = fftshift(ifftn(ifftshift(PSFIncoherent.* PSFdet)));
    xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
    SWOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;

[~,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
    yindex = 1-(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:))));
    LatticeyFWHM1 = find(yindex,1,'first');
    LatticeyFWHM2 = find(yindex,1,'last');
    overallOTF = fftshift(ifftn(ifftshift(LatticePSFDithered.* PSFdet)));
    xzOTF = abs(squeeze(overallOTF(:,(N+1)/2,:)))/max(abs(squeeze(overallOTF(:,(N+1)/2,:))),[],'all');
    LatticeOTFmask = xzOTF((N+1)/2:end,(N+1)/2:end) >= OTFthreshold;

AberratedSWPupil = SWPupil.*ComplexPhase;
AberratedLatticePupil= LatticePupil.*ComplexPhase;

[AberratedPSFCoherent,AberratedPSFIncoherent,AberratedSWcenter] = SimulateSWPair(AberratedSWPupil);
AberratedPSFIncoherent = AberratedPSFIncoherent./SWcenter(2,1);
AberratedPSFCoherent = AberratedPSFCoherent./SWcenter(1,1);
[AberratedLatticePSF,AberratedLatticePSFDithered,~] = SimulateLattice(AberratedLatticePupil);
AberratedLatticePSFDithered = AberratedLatticePSFDithered./Latticecenter(2,1);
AberratedLatticePSF = AberratedLatticePSF/Latticecenter(1,1);

OverallAberratedPSFIncoherent = AberratedPSFIncoherent.* PSFdet;
OverallAberratedLatticePSFDithered = AberratedLatticePSFDithered.* PSFdet;

%% Analysis
SRatioSW = max(OverallAberratedPSFIncoherent(:,:,(N+1)/2),[],'all');
SRatioSW_corrected = (AberratedSWcenter(3,1) + AberratedSWcenter(4,1))/SWcenter(2,1);
SRatioLattice = max(OverallAberratedLatticePSFDithered(:,:,(N+1)/2),[],'all');

[AberratedSWAveragefc3,AberratedSWAveragefc3FWHM,~,~] = RunFC3(AberratedPSFIncoherent,PSFdet,SWOTFmask,Iter,SNR,SWyFWHM2);
[AberratedLatticeAveragefc3,AberratedLatticeAveragefc3FWHM,~,~] = RunFC3(AberratedLatticePSFDithered,PSFdet,LatticeOTFmask,Iter,SNR,LatticeyFWHM2);

[SWconvLines,~,~,~] = ConvRes(PSFIncoherent,PSFdet,SNR);
[LatticeconvLines,~,~,~] = ConvRes(LatticePSFDithered,PSFdet,SNR);
[AberratedSWconvLines,lineSpot,Spacing,LineZ] = ConvRes(AberratedPSFIncoherent,PSFdet,SNR);
[AberratedLatticeconvLines,~,~,~] = ConvRes(AberratedLatticePSFDithered,PSFdet,SNR);

GTBeads = load('/Users/andyzjc/Dropbox (Princeton)/Tian-Ming for Andy/Manuscripts/Adaptive Polarization Controlled Incoherent SW Light Sheet/Figures/Figure 5/hex_0.58_0.04_tophat_V1/AberrationAnalysis/BeadsImage/GroundTruth.mat');
Beads = GTBeads.GTallbeads;
Vol = zeros(N,N,N);
Vol(:,:,(N+1)/2) = Beads;

fftvol = fftshift(fftn(ifftshift(Vol))) .* fftshift(fftn(ifftshift(OverallAberratedPSFIncoherent)));
temp = abs(fftshift(ifftn(ifftshift(fftvol))));
AberratedSWBeadsVol = temp + poissrnd(temp) .* 1/SNR;

fftvol = fftshift(fftn(ifftshift(Vol))) .* fftshift(fftn(ifftshift(OverallAberratedLatticePSFDithered)));
temp = abs(fftshift(ifftn(ifftshift(fftvol))));
AberratedLatticeBeadsVol = temp + poissrnd(temp) .* 1/SNR;

%% Plot

grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,AberratedLatticePSF,AberratedLatticePSFDithered,AberratedLatticePupil,PSFdet,savingdir)
grapthSW(NA1,deltaNA,LatticeType,SWweighting,AberratedPSFIncoherent,AberratedPSFCoherent,AberratedSWPupil,PSFdet,savingdir)

fig1 = figure;
    imagesc(Z_exc,X_exc,AberratedSWBeadsVol(:,:,(N+1)/2))
    axis image
    xlim([-40,40])
    ylim([-40,40])
    clim([0,1])
    colormap(hot)
    colorbar
    print(fig1, '-dsvg', [ savingdir 'AberratedSWBeads' '.SVG'],'-r300')
    print(fig1, '-dpng', [ savingdir 'AberratedSWBeads' '.PNG'],'-r300')

fig2 = figure;
    imagesc(Z_exc,X_exc,AberratedLatticeBeadsVol(:,:,(N+1)/2))
    axis image
    xlim([-40,40])
    ylim([-40,40])
    colormap(hot)
    clim([0,1])
    colorbar
    print(fig2, '-dsvg', [  savingdir 'LatticeBeads' '.SVG'],'-r300')
    print(fig2, '-dpng', [  savingdir 'LatticeBeads' '.PNG'],'-r300')

fig3 = figure;
h1 = subplot(2,2,1,'Parent',fig3);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3)
    title("SW")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    % colormap(h1,fire(256))
    colormap(h1,turbo)
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

h1 = subplot(2,2,2,'Parent',fig3);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3)
    title("LLS")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    % colormap(h1,fire(256))
    colormap(h1,turbo)
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

h1 = subplot(2,2,3,'Parent',fig3);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3FWHM)
    title("SW FWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    % colormap(h1,fire(256))
    colormap(h1,turbo)
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

h1 = subplot(2,2,4,'Parent',fig3);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3FWHM)
    title("LLS FWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    % colormap(h1,fire(256))
    colormap(h1,turbo)
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

print(fig3, '-dsvg', [  savingdir 'FC3' '.SVG'],'-r300')
print(fig3, '-dpng', [  savingdir 'FC3' '.PNG'],'-r300')

fig4 = figure;
    plot(LineZ,SWconvLines((N+1)/2,:),'magenta','LineWidth',0.5)
    hold on
    plot(LineZ,LatticeconvLines((N+1)/2,:),'Color','b','LineWidth',0.5)
    plot(LineZ,AberratedSWconvLines((N+1)/2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',0.5)
    plot(LineZ,AberratedLatticeconvLines((N+1)/2,:),'Color',[0.8500 0.3250 0.0980],'LineWidth',0.5)
    xline(LineZ(lineSpot),'k')
    xlim([0,50])
    ylim([0,50])
    legend("pSW","LLS","aberrated iSW","aberrated LLS")
    xlabel("z(lambda/n)")
    ylabel("Intensity (a.u)")
    hold off 
    grid on
    pbaspect([7 2 1])
print(fig4, '-dsvg', [  savingdir 'LineGrating' '.SVG'],'-r300')
print(fig4, '-dpng', [  savingdir 'LineGrating' '.PNG'],'-r300')

fig5 = figure;
    imagesc(LineZ,LineZ,AberratedSWconvLines)
    colormap(hot)
    hold on
    axis image
    xline(LineZ(lineSpot),'k','LineWidth',1)
    xlim([0,50])
    ylim([0,20])
    clim([0,50])
    colorbar
    xlabel("z(um)")
    ylabel("um")
    hold off
print(fig5, '-dsvg', [  savingdir 'SWLines' '.SVG'],'-r300')
print(fig5, '-dpng', [  savingdir 'SWLines' '.PNG'],'-r300')    

fig6 = figure;
    imagesc(LineZ,LineZ,AberratedLatticeconvLines)
    colormap(hot)
    hold on
    xline(LineZ(lineSpot),'k','LineWidth',1)
    axis image
    xlim([0,50])
    ylim([0,20])
    clim([0,50])
    colorbar
    xlabel("z(um)")
    ylabel("um")
    hold off
print(fig6, '-dsvg', [  savingdir 'LLSLines' '.SVG'],'-r300')
print(fig6, '-dpng', [  savingdir 'LLSLines' '.PNG'],'-r300')  
close all

%% decomposition
MinRadialOrder = 1;
MaxRadialOrder = 6;

LabelArray = cell(length(RadioOrderArray),1);
for i = 1:length(RadioOrderArray)
    LabelArray{i,1} = "Z^{" + num2str(AngularFrequencyArray(i)) + "}_{" + num2str(RadioOrderArray(i)) + "}";
end

counter = 1;
RadioOrderArray = [];
AngularFrequencyArray =[];
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(1,counter) = i;
        AngularFrequencyArray(1,counter) = AngularFrequency(k);

        temp = zeros(N,N);
        temp(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx));
        factor_top(:,counter) = max(temp(idx)) - min(temp(idx));

        Z(:,counter) = temp(idx);
        counter = counter +1;
    end
end

temp = IntScaledwavefront(idx);
a_mn = Z\temp(:); % nm
b_mn_A = a_mn.*factor_top';
ZCoeff = b_mn_A;

counter = 1;
composite_wavefront = zeros(N,N);
Z_iter = zeros(N,N);
for i = 1:3
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        Z_iter(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx));
        composite_wavefront(idx) = composite_wavefront(idx) +  (ZCoeff(counter)) .* Z_iter(idx);
        counter = counter + 1;
    end
end
figure
imagesc(KZ_exc,KX_exc,composite_wavefront)
xlabel("k_x/(4\pin/\lambda_{exc})");
ylabel("k_z/(4\pin/\lambda_{exc})");
axis image
title("n=1-3")
colormap(turbo)
colorbar
clim([-0.5 0.5])

for i = 4:6
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        Z_iter(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx));
        composite_wavefront(idx) = composite_wavefront(idx) +  (ZCoeff(counter)) .* Z_iter(idx);
        counter = counter + 1;
    end
end
figure
imagesc(KZ_exc,KX_exc,composite_wavefront)
xlabel("k_x/(4\pin/\lambda_{exc})");
ylabel("k_z/(4\pin/\lambda_{exc})");
axis image
clim([-0.5 0.5])
title("n=1-3+n=4-6")
colormap(turbo)
colorbar

fig1 = figure;
    h1 = subplot(1,1,1,'Parent',fig1);
    bar(1:length(ZCoeff),ZCoeff,0.75)
    xlabel("Zernike Mode")
    ylabel("Mode Amplitude (um)")
    title("Naji science advance 2020, fig6")
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 8;
    pbaspect([6 2 1])
    hold off
    ylim([-0.5,0.5])
    grid on
