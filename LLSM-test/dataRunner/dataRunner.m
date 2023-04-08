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

% xzPSFdet = PSFdet(:,:,(N+1)/2);
% yzPSFdet = squeeze(PSFdet(:,(N+1)/2,:)); 
% xzOTFdet = fftshift(fft2(ifftshift(xzPSFdet)));
% zOTFdet = real(xzOTFdet(:,(N+1)/2));

%% Unaberrated
clc
NA1 = 0.42;
deltaNA = 0.08;
LatticeType = 'hex';
ProfileType = 'gaussian';
SWweighting = 7/10; %7/5 for equal OTF
Latticeweighting = 1;
SNR = 10;
Iter = 10;

savingdir = [LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_' ProfileType '/'];
mkdir(savingdir) 
SWsavingdir = [savingdir 'SW_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_WR' num2str(SWweighting) '/'];
mkdir(SWsavingdir) 
Latticesavingdir = [savingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_WR' num2str(Latticeweighting) '/'];
mkdir(Latticesavingdir) 
SWDatasavingdir = [savingdir 'SWdata/'];
mkdir(SWDatasavingdir) 
LLSDatasavingdir = [savingdir 'LLSdata/'];
mkdir(LLSDatasavingdir) 

if isequal(LatticeType,'hex')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,NA1/2,...
    deltaNA,2*deltaNA,...
    SWweighting);
elseif isequal(LatticeType,'SW')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,0,...
    deltaNA,0,...
    SWweighting);
    SWPupil = SWPupil .* Phase_factor;
elseif isequal(LatticeType,'NC')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,NA1*3/4,...
    deltaNA,deltaNA*4/3,...
    SWweighting);
elseif isequal(LatticeType,'Airy')
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,0,0,...
    0,NA1,...
    SWweighting);
else
    [SWPupil,~,~] = GetSWPairPupil(ProfileType,NA1,sqrt(2*NA1*deltaNA)/2,...
    deltaNA,sqrt(2*NA1*deltaNA),...
    SWweighting);
end
[PSFCoherent,PSFIncoherent,SWcenter] = SimulateSWPair(SWPupil);
PSFCoherent = PSFCoherent/max(PSFCoherent,[],'all');
PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
mkdir([SWsavingdir 'Unaberrated/'])
grapthSW(NA1,deltaNA,LatticeType,SWweighting,PSFIncoherent,PSFCoherent,SWPupil,PSFdet,[SWsavingdir 'Unaberrated/'])
    yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
    SWyFWHM1 = find(yindex,1,'first') ;
    SWyFWHM2 = find(yindex,1,'last');
[SWAveragefc3,SWAveragefc3FWHM,~,~] = RunFC3(PSFIncoherent,PSFdet,Iter,SNR,SWyFWHM2);
[SWconvLines,lineSpot,~,Line_Z] = ConvRes(PSFIncoherent,PSFdet);
save([SWDatasavingdir '/PSFIncoherent.mat'],'PSFIncoherent')
save([ SWDatasavingdir '/PSFIncoherent_Center.mat'], 'SWcenter')
save([ SWDatasavingdir '/PSFIncoherent_FC3.mat'], 'SWAveragefc3')
save([ SWDatasavingdir '/PSFIncoherent_FC3FWHM.mat'], 'SWAveragefc3FWHM')
save([ SWDatasavingdir '/PSFIncoherent_LineGrating.mat'], 'SWconvLines')

% Lattice
[LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
NA1,deltaNA, ...
0.6,0,...
Latticeweighting);
[LatticePSF,LatticePSFDithered,Latticecenter] = SimulateLattice(LatticePupil);
LatticePSF = LatticePSF/max(LatticePSF,[],'all');
LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
mkdir([Latticesavingdir 'Unaberrated/'])
grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,LatticePSF,LatticePSFDithered,LatticePupil,PSFdet,[Latticesavingdir 'Unaberrated/'])
    yindex = 1-(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:))));
    LatticeyFWHM1 = find(yindex,1,'first') ;
    LatticeyFWHM2 = find(yindex,1,'last');
[LatticeAveragefc3,LatticeAveragefc3FWHM,~,~] = RunFC3(LatticePSFDithered,PSFdet,Iter,SNR,LatticeyFWHM2);
[LatticeconvLines,~,~,~] = ConvRes(LatticePSFDithered,PSFdet);
save([ LLSDatasavingdir '/LatticePSF.mat'], 'LatticePSF')
save([ LLSDatasavingdir '/LatticePSFDithered.mat'], 'LatticePSFDithered')
save([ LLSDatasavingdir '/Latticecenter.mat'], 'Latticecenter')
save([ LLSDatasavingdir '/LatticePSFDithered_FC3.mat'], 'LatticeAveragefc3')
save([ LLSDatasavingdir '/LatticePSFDithered_FC3FWHM.mat'], 'LatticeAveragefc3FWHM')
save([ LLSDatasavingdir '/LatticePSFDithered_LineGrating.mat'], 'LatticeconvLines')

%% Aberrated
clc
[theta,r] = cart2pol(kx_exc./(0.6./n*k_wave),kz_exc./(0.6./n*k_wave));
idx = r<=1;

MinRadialOrder = 0;
MaxRadialOrder = 6;
PhaseAmplitude = 4; %default 4 
counter = 1;

RadioOrderArray = zeros(28,1);
AngularFrequencyArray = RadioOrderArray;

% AberratedPSFCoherent = cell(28,1);
% AberratedSWPupil = AberratedPSFCoherent;
AberratedPSFIncoherent = cell(28,1);
AberratedSWcenter = AberratedPSFIncoherent;
% AberratedLatticePupil = AberratedPSFCoherent;
AberratedLatticePSF = AberratedPSFIncoherent;
AberratedLatticePSFDithered = AberratedPSFIncoherent;
% AberratedLatticecenter = AberratedPSFCoherent;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        RadioOrderArray(counter) = i;
        AngularFrequencyArray(counter) = AngularFrequency(k);
        phase = zeros(size(kx_exc));
        phase(idx) = zernfun(i,AngularFrequency(k),r(idx),theta(idx),'norm');

        %SW 
        AberratedPupil1 = zeros(size(phase));
        AberratedPupil2 = zeros(size(phase));
        SWPupil1 = squeeze(SWPupil(:,:,1));
        SWPupil2 = squeeze(SWPupil(:,:,2));

        AberratedPupil1(idx) = SWPupil1(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
        AberratedPupil2(idx) = SWPupil2(idx) .* exp(PhaseAmplitude.* 1i.*phase(idx));
        temp1(:,:,1) = AberratedPupil1;
        temp1(:,:,2) = AberratedPupil2;
        [temp3,temp4,temp5] = SimulateSWPair(temp1);
        % AberratedSWPupil{counter,1} = temp1;
        % AberratedPSFCoherent{counter,1} = temp3 /SWcenter(1,1);
        AberratedPSFIncoherent{counter,1} = temp4 /SWcenter(2,1);
        AberratedSWcenter{counter,1} = temp5;
        SWAberratedsavingdir = [SWsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_Amplitude_4/'];
        mkdir(SWAberratedsavingdir)
        grapthSW(NA1,deltaNA,LatticeType,SWweighting,temp4/SWcenter(2,1),temp3/SWcenter(1,1),temp1,...
           PSFdet,[SWAberratedsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_']);

        %Lattice
        temp2 = zeros(size(phase));
        temp2(idx) = LatticePupil(idx) .* exp(PhaseAmplitude.* 1i .*phase(idx));
        [temp6,temp7,~] = SimulateLattice(temp2);
        % AberratedLatticePupil{counter,1} = temp2;
        AberratedLatticePSF{counter,1} = temp6/Latticecenter(1,1);
        AberratedLatticePSFDithered{counter,1} = temp7/Latticecenter(2,1);
        % AberratedLatticecenter{counter,1} = temp8;
        LatticeAberratedsavingdir = [Latticesavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_Amplitude_4/'];
        mkdir(LatticeAberratedsavingdir)
        grapthLattice(NA1,deltaNA,LatticeType,Latticeweighting,temp6/Latticecenter(1,1),temp7/Latticecenter(2,1),temp2,...
            PSFdet,[LatticeAberratedsavingdir 'Z_' num2str(i) '_' num2str(AngularFrequency(k)) '_'])
        
        counter = counter +1;
    end
end
% save([SWDatasavingdir '/AberratedPSFCoherent.mat'],'AberratedPSFCoherent')
save([SWDatasavingdir '/AberratedPSFIncoherent.mat'],'AberratedPSFIncoherent')
save([SWDatasavingdir '/AberratedSWcenter.mat'],'AberratedSWcenter')
% save([SWDatasavingdir '/AberratedSWPupil.mat'],'AberratedSWPupil')

% save([ LLSDatasavingdir '/AberratedLatticePupil.mat'], 'AberratedLatticePupil')
save([ LLSDatasavingdir '/AberratedLatticePSF.mat'], 'AberratedLatticePSF')
save([ LLSDatasavingdir '/AberratedLatticePSFDithered.mat'], 'AberratedLatticePSFDithered')

clear LatticeAberratedsavingdir SWAberratedsavingdir
clear temp1 temp2 temp3 temp4 temp5 temp6 temp7 temp8 phase
clear AberratedPupil1 AberratedPupil2
clear SWsavingdir Latticesavingdir

%% Analysis
SW_SRatio_Focal = zeros(size(AberratedPSFIncoherent));
SW_SRatio_FWHM = SW_SRatio_Focal;
SW_SRatio_3D = SW_SRatio_Focal;
SW_SRatio_corrected = SW_SRatio_Focal;

AberratedSWAveragefc3 = cell(28,1);
AberratedSWAveragefc3FWHM = AberratedSWAveragefc3;
AberratedSWconvLines = AberratedSWAveragefc3;
for i = 1:length(AberratedPSFIncoherent)
    tempPSF = AberratedPSFIncoherent{i,1};
    tempSWcenter = AberratedSWcenter{i,1};

    % Srethl ratio
    SW_SRatio_Focal(i,1) = tempPSF((N+1)/2,(N+1)/2,(N+1)/2);
    SW_SRatio_FWHM(i,1) = tempPSF((N+1)/2,(N+1)/2,SWyFWHM2);
    SW_SRatio_3D(i,1) = max(tempPSF,[],'all');

    % aberration correction (quick guess) 
    SW_SRatio_corrected(i,1) = (tempSWcenter(3,1) + tempSWcenter(4,1))/SWcenter(2,1);

    % FC3 
    [AberratedSWAveragefc3{i,1},AberratedSWAveragefc3FWHM{i,1},~,~] = RunFC3(tempPSF,PSFdet,Iter,SNR,SWyFWHM2);

    % line grating
    [AberratedSWconvLines{i,1},~,~,~] = ConvRes(tempPSF,PSFdet);
end
save([SWDatasavingdir '/SW_SRatio_Focal.mat'],'SW_SRatio_Focal')
save([SWDatasavingdir '/SW_SRatio_FWHM.mat'],'SW_SRatio_FWHM')
save([SWDatasavingdir '/SW_SRatio_3D.mat'],'SW_SRatio_3D')
save([SWDatasavingdir '/SW_SRatio_corrected.mat'],'SW_SRatio_corrected')
save([SWDatasavingdir '/AberratedSWAveragefc3.mat'],'AberratedSWAveragefc3')
save([SWDatasavingdir '/AberratedSWAveragefc3FWHM.mat'],'AberratedSWAveragefc3FWHM')
save([SWDatasavingdir '/AberratedSWconvLines.mat'],'AberratedSWconvLines')

Lattice_SRatio_Focal = zeros(size(AberratedLatticePSFDithered));
Lattice_SRatio_FWHM = Lattice_SRatio_Focal;
Lattice_SRatio_3D = Lattice_SRatio_Focal;

AberratedLatticeAveragefc3 = AberratedSWAveragefc3;
AberratedLatticeAveragefc3FWHM = AberratedLatticeAveragefc3;
AberratedLatticeconvLines = AberratedLatticeAveragefc3;
for i = 1:length(AberratedLatticePSFDithered)
    tempPSF = AberratedLatticePSFDithered{i,1};

    % Srethl ratio
    Lattice_SRatio_Focal(i,1) = tempPSF((N+1)/2,(N+1)/2,(N+1)/2);
    Lattice_SRatio_FWHM(i,1) = tempPSF((N+1)/2,(N+1)/2,LatticeyFWHM2);
    Lattice_SRatio_3D(i,1) = max(tempPSF,[],'all');

    % FC3 
    [AberratedLatticeAveragefc3{i,1},AberratedLatticeAveragefc3FWHM{i,1},~,~] = RunFC3(tempPSF,PSFdet,Iter,SNR,LatticeyFWHM2);

    % line grating
    [AberratedLatticeconvLines{i,1},~,~,~] = ConvRes(tempPSF,PSFdet);
end
save([ LLSDatasavingdir '/Lattice_SRatio_Focal.mat'], 'Lattice_SRatio_Focal')
save([ LLSDatasavingdir '/Lattice_SRatio_FWHM.mat'], 'Lattice_SRatio_FWHM')
save([ LLSDatasavingdir '/Lattice_SRatio_3D.mat'], 'Lattice_SRatio_3D')
save([ LLSDatasavingdir '/AberratedLatticeAveragefc3.mat'], 'AberratedLatticeAveragefc3')
save([ LLSDatasavingdir '/AberratedLatticeAveragefc3FWHM.mat'], 'AberratedLatticeAveragefc3FWHM')
save([ LLSDatasavingdir '/AberratedLatticeconvLines.mat'], 'AberratedLatticeconvLines')

clear tempPSF tempSWcenter

%% Pretty plots, Strehl
clc
LabelArray = cell(28,1);
for i = 1:length(RadioOrderArray)
    LabelArray{i,1} = "Z^{" + num2str(AngularFrequencyArray(i)) + "}_{" + num2str(RadioOrderArray(i)) + "}";
end

Analysis_savingdir = [savingdir 'AberrationAnalysis/'];
mkdir(Analysis_savingdir) 

Strehlsavingdir = [Analysis_savingdir 'StrehlRatio/'];
mkdir(Strehlsavingdir)

fig1 = figure;
    h1 = subplot(1,1,1,'Parent',fig1);
    plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal,1,'Parent',h1,'LineStyle','-o');
    hold on
    plot(1:length(Lattice_SRatio_Focal),Lattice_SRatio_Focal,1,'Parent',h1,'LineStyle','-o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("Focal Plane")
    pbaspect([5 1 1])
    hold off
    print(fig1, '-dsvg', [Strehlsavingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_FocalStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [Strehlsavingdir LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_FocalStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig2 = figure;
    h1 = subplot(1,1,1,'Parent',fig2);
    plot(1:length(SW_SRatio_FWHM),SW_SRatio_FWHM,1,'Parent',h1,'LineStyle','-o');
    hold on
    plot(1:length(Lattice_SRatio_FWHM),Lattice_SRatio_FWHM,1,'Parent',h1,'LineStyle','-o');
    % p1 = plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal./Lattice_SRatio_Focal,'magenta');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("yFWHM Plane")
    pbaspect([5 1 1])
    hold off
    print(fig2, '-dsvg', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yFWHMStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig2, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_yFWHMStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig3 = figure;
    h1 = subplot(1,1,1,'Parent',fig3);
    plot(1:length(SW_SRatio_3D),SW_SRatio_3D,1,'Parent',h1,'LineStyle','-o');
    hold on
    plot(1:length(Lattice_SRatio_3D),Lattice_SRatio_3D,1,'Parent',h1,'LineStyle','-o');
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("3D")
    pbaspect([5 1 1])
    hold off
    print(fig3, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_3DStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig3, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_3DStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')

fig4 = figure;
    h1 = subplot(1,1,1,'Parent',fig4);
    plot(1:length(SW_SRatio_Focal),SW_SRatio_Focal,1,'Parent',h1,'LineStyle','-o');
    hold on
    plot(1:length(SW_SRatio_corrected),SW_SRatio_corrected,1,'Parent',h1,'LineStyle','-o');
    lgd = legend(b1,"Aberrated","Corrected");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("Correction")
    pbaspect([5 1 1])
    hold off
    print(fig4, '-dsvg', [ Strehlsavingdir   LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_CorrectedStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig4, '-dpng', [ Strehlsavingdir  LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_CorrectedStrehl_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
close all
clear h1 b1 lgd fig1 fig2 fig3 fig4 Strehlsavingdir LabelArray

%% FC3 
clc
SWFC3savingdir = [Analysis_savingdir 'FC3/iSW/'];
mkdir(SWFC3savingdir)
fig1 = figure;
for i = 1:length(AberratedSWAveragefc3)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),SWAveragefc3)
    title("Unaberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = SWAveragefc3./AberratedSWAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("Unaberrated/Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),SWAveragefc3FWHM)
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3FWHM{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = SWAveragefc3FWHM./AberratedSWAveragefc3FWHM{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ SWFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ SWFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_SWFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all

LatticeFC3savingdir = [Analysis_savingdir 'FC3/LLS/'];
mkdir(LatticeFC3savingdir)
fig1 = figure;
for i = 1:length(AberratedLatticeAveragefc3)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),LatticeAveragefc3)
    title("Unaberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    title("Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = LatticeAveragefc3./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("Unaberrated/Aberrated")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),LatticeAveragefc3FWHM)
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3FWHM{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = LatticeAveragefc3FWHM./AberratedLatticeAveragefc3FWHM{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ LatticeFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ LatticeFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LatticeFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all

RatioMapFC3savingdir = [Analysis_savingdir 'FC3/RatioMap/'];
mkdir(RatioMapFC3savingdir)
fig1 = figure;
for i = 1:length(AberratedLatticeAveragefc3)
    h1 = subplot(2,3,1,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("SW")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,2,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    title("LLS")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = AberratedSWAveragefc3{i,1}./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,3,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    title("SW/LLS")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,4,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedSWAveragefc3{i,1})
    title("yFWHM")
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    h1 = subplot(2,3,5,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),AberratedLatticeAveragefc3{i,1})
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,fire(256))
    colorbar
    clim([0,1])
    set(gca,'YDir','normal')
    axis image

    Ratiomap = AberratedSWAveragefc3{i,1}./AberratedLatticeAveragefc3{i,1};
    h1 = subplot(2,3,6,'Parent',fig1);
    imagesc(KX_exc((N+1)/2:N),KZ_exc((N+1)/2:N),Ratiomap)
    xlabel("k_r/(4\pin/\lambda_{exc})")
    ylabel("k_z/(4\pin/\lambda_{exc})")
    colormap(h1,hot(256))
    colorbar
    clim([0,10])
    set(gca,'YDir','normal')
    axis image

    print(fig1, '-dsvg', [ RatioMapFC3savingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_RatioMapFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ RatioMapFC3savingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_RatioMapFC3_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all
clear SWFC3savingdir LatticeFC3savingdir RatioMapFC3savingdir
clear Ratiomap 

%% line grating
clc
lineGratesavingdir = [Analysis_savingdir 'ConvRes/'];
mkdir(lineGratesavingdir)
fig1 = figure;
for i = 1:length(AberratedSWconvLines)
    temp1 = AberratedSWconvLines{i,1};
    temp2 = AberratedLatticeconvLines{i,1};
    plot(Line_Z,SWconvLines((N+1)/2,:),'magenta','LineWidth',2)
    hold on
    plot(Line_Z,LatticeconvLines((N+1)/2,:),'Color','b','LineWidth',2,'LineStyle','-.')
    plot(Line_Z,temp1((N+1)/2,:),'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
    plot(Line_Z,temp2((N+1)/2,:),'Color',[0.8500 0.3250 0.0980],'LineStyle','-.','LineWidth',2)
    xline(Line_Z(lineSpot),'k')
    xlim([0,20])
    legend("iSW","LLS","aberrated iSW","aberrated LLS")
    title("Spacing increment=" + num2str(deltax*1000) + "nm")
    xlabel("z(um)")
    ylabel("Intensity (a.u)")
    hold off 
    grid on
    pbaspect([7 2 1])

    print(fig1, '-dsvg', [ lineGratesavingdir  'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LineGrating_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.SVG'],'-r300')
    print(fig1, '-dpng', [ lineGratesavingdir 'Z_' num2str(RadioOrderArray(i)) '_' num2str(AngularFrequencyArray(i)) '_' LatticeType '_' num2str(NA1) '_' num2str(deltaNA) '_LineGrating_SWweight' num2str(SWweighting) '_LLSweight' num2str(Latticeweighting) '.PNG'],'-r300')
end
close all
clear lineGratesavingdir temp1 temp2

%%
% figure
% plot(KZ_exc,log(zOTFexc2),'g','LineWidth',2,'LineStyle','-')
%     hold on
%     % axis image
% plot(KZ_exc,log(zOTFexc1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
%     xlabel("k_z/(4\pin/\lambda_{exc})")
%     ylabel("ln(OTF Strength)")
%     xticks([-0.5,-0.25,0,0.25,0.5]);
%     yticks([])
%     % yticks([0,0.5,1]);
%     xlim([-0.5,0.5])
%     ylim([-7,0])
%     grid on
% 
% figure
% plot(Z_exc,xzPSFexc2(:,(N+1)/2),'g','LineWidth',2,'LineStyle','-')
%     hold on
%     % axis image
% plot(Z_exc,xzPSFexc1(:,(N+1)/2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1)
%     xlabel("z/(\lambda_{exc}/n)")
%     ylabel("Intensity (a.u)")
%     xticks([-10,-5,0,5,10]);
%     yticks([0,0.25,0.5,0.75,1]);
%     xlim([-10,10])
%     ylim([0,1])
%     grid on