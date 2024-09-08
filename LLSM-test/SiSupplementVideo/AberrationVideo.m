%% 
clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%%
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

%% Unaberrated 
clc
NA1 = 0.58; %0.6 for gaussian
deltaNA = 0.04;
LatticeType = 'hex';
ProfileType = 'tophat';
SWweighting = 7/10; %4/3 for equal OTF V2 LLS, 1/sqrt(2) for V1 LLS
Latticeweighting = [1,1,1,1,1,1]; % 1.9 for V2 LLS
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

PSFIncoherent1 = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,1) )) ) ).^2;
PSFIncoherent2 = abs( fftshift( ifft2(ifftshift(SWPupil(:,:,2) )) ) ).^2;
PSFIncoherent = PSFIncoherent2 + PSFIncoherent1;
PEARLSmax = max(PSFIncoherent,[],'all');
LatticePSF = abs( fftshift( ifft2(ifftshift(LatticePupil )) ) ).^2;
LatticePSFDithered = meshgrid(mean(LatticePSF,2))';
LLSmax = max(LatticePSFDithered,[],'all');

%% Aberration
[theta,r] = cart2pol(kx_exc./(PupilNA./n*k_wave),kz_exc./(PupilNA./n*k_wave));
idx = r<=1;

MinRadialOrder = 2;
MaxRadialOrder = 6;
PhaseAmplitude = 6*wavelength_exc/(2*pi); 

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

LabelArray = cell(length(RadioOrderArray),1);
for i = 1:length(RadioOrderArray)
    LabelArray{i,1} = "Z^{" + num2str(AngularFrequencyArray(i)) + "}_{" + num2str(RadioOrderArray(i)) + "}";
end

%% Making video
outputVideo = VideoWriter('outputVideo.mp4', 'MPEG-4'); % Create a video writer object
outputVideo.FrameRate = 1; % Set the frame rate (frames per second)
open(outputVideo); % Open the video writer object
firemap = fire(256);

counter = 1;
for i = MinRadialOrder:MaxRadialOrder
    AngularFrequency = -i:2:i;
    for k = 1:length(AngularFrequency)
        AberratedPSFIncoherent = [];
        AberratedLatticePSFDithered = [];
        fig = figure;
        
        hf=colordef(fig,'black'); %Set color scheme
        hf.Color='k'; %Set background color of figure window

        fig.WindowState = 'maximized';
        phase = GetSingleZmodePupil(i,AngularFrequency(k),PupilNA);
        ComplexPhase = exp(PhaseAmplitude .* 1i .* phase);

        AberratedSWPupil = zeros(size(phase));
        AberratedLLSPupil = zeros(size(phase));

        AberratedPupil1 = zeros(size(phase));
        AberratedPupil2 = zeros(size(phase));
        SWPupil1 = squeeze(SWPupil(:,:,1));
        SWPupil2 = squeeze(SWPupil(:,:,2));

        AberratedPupil1 = SWPupil1 .* ComplexPhase;
        AberratedPupil2 = SWPupil2 .* ComplexPhase;
        
        Profile_pupil1 = abs( fftshift( ifft2(ifftshift(AberratedPupil1 )) ) ).^2;
        Profile_pupil2 = abs( fftshift( ifft2(ifftshift(AberratedPupil2 )) ) ).^2;
        AberratedPSFIncoherent = Profile_pupil1 + Profile_pupil2;
        AberratedPSFIncoherent = AberratedPSFIncoherent./PEARLSmax;
        % PEARLSStrehl(counter,1) = AberratedPSFIncoherent((N+1)/2,(N+1)/2);
        
        AberratedLLSPupil = LatticePupil .* ComplexPhase;
        AberratedLatticePSF = abs( fftshift( ifft2(ifftshift(AberratedLLSPupil )) ) ).^2;
        AberratedLatticePSFDithered = meshgrid(mean(AberratedLatticePSF,2))';
        AberratedLatticePSFDithered = AberratedLatticePSFDithered./LLSmax;
        % LLSStrehl(counter,1) = AberratedLatticePSFDithered((N+1)/2,(N+1)/2);

        PEARLSStrehl(counter,1) = max(AberratedPSFIncoherent.*PSFdet(:,:,(N+1)/2),[],'all'); % only at focal point/ on optical axis
        LLSStrehl(counter,1) = max(AberratedLatticePSFDithered.*PSFdet(:,:,(N+1)/2),[],'all');
        
        s1 = subplot(2,3,1);
        imagesc(KX_exc,KZ_exc,phase/2/pi/3)
        axis image
        title("Zernike mode, Z=" + num2str(i) + "," + num2str(AngularFrequency(k)))
        xlim([-0.5,0.5])
        ylim([-0.5,0.5])
        xlabel("k_x/(4\pin/\lambda_{exc})");
        ylabel("k_z/(4\pin/\lambda_{exc})");
        s1.Colormap = jet;
        s1.XAxis.FontSize = 15;
        s1.YAxis.FontSize = 15;
        s1.Title.FontSize = 15;
        s1Pos = get(s1,'position');
        bar = colorbar;
        bar.Location = 'eastoutside';
        bar.Limits = [-2.5,2.5];
        set(s1,'position',s1Pos)

        s2 = subplot(2,3,2);
        imagesc(X_exc,Z_exc,AberratedLatticePSFDithered);
        title("LLS xzPSF")
        axis image;
        xlim([-20,20])
        ylim([-20,20])
        xlabel("x/(\lambda_{exc}/n)")
        ylabel("z/(\lambda_{exc}/n)")
        s2.Colormap = firemap;
        clim([0,1])
        s2.XAxis.FontSize = 15;
        s2.YAxis.FontSize = 15;
        s2.Title.FontSize = 15;

        s3 = subplot(2,3,3);
        imagesc(X_exc,Z_exc,AberratedPSFIncoherent);
        title("PEARLS xzPSF")
        axis image;
        xlim([-20,20])
        ylim([-20,20])
        xlabel("x/(\lambda_{exc}/n)")
        ylabel("z/(\lambda_{exc}/n)")
        s3.Colormap = firemap;
        clim([0,1])
        s3Pos = get(s3,'position');
        bar = colorbar;
        bar.Location = 'eastoutside';
        set(s3,'position',s3Pos)
        s3.XAxis.FontSize = 15;
        s3.YAxis.FontSize = 15;
        s3.Title.FontSize = 15; 

        h1 = subplot(2,3,4:6);
        scatter(1:counter,LLSStrehl',200,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980]);
        hold on
        scatter(1:counter,PEAzRLSStrehl',200,'filled','MarkerFaceColor','blue');
        lgd = legend("LLS","PEARLS");
        lgd.Location = 'northoutside';
        lgd.Orientation = 'horizontal';
        lgd.FontSize = 15;
        h1.XAxis.TickValues = 1:length(RadioOrderArray);
        h1.XAxis.TickLabels = LabelArray;
        h1.XAxis.FontSize = 15;
        h1.YAxis.FontSize = 15;
        ylim([0,1])
        xlim([1,25])
        grid on
        xlabel("Aberration Mode")
        ylabel("Strehl Ratio")
        pbaspect([5 1 1])
        hold off

        pause(2)
        counter = counter +1;
        % Capture the frame
        frame = getframe(gcf);
        writeVideo(outputVideo, frame); % Write the frame to the video
        close(gcf); % Close the figure
    end
end


% Close the video file
close(outputVideo);
disp('Video creation complete.');