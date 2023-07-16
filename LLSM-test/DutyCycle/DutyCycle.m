clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%%
LatticeType = 'square';
ProfileType = 'tophat';
NA1 = 0.3:0.01:0.6;
deltaNA1 = 0.02:0.02:0.1;
ROI_range = 10;
%% LLS 

fig1 = figure;
    hold on
for j = 1:length(deltaNA1)
    for i = 1:length(NA1)
        [LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,ProfileType, ...
                                                                     NA1(i),deltaNA1(j), ...
                                                                     0.8,0.0,...
                                                                     [1 1 1 1 1 1]);
        xzLatticePSF = abs(fftshift(ifft2(ifftshift(LatticePupil)))).^2;
        xzLatticePSF = xzLatticePSF./max(xzLatticePSF,[],'all');

        xzLatticePSFDither = meshgrid(max(xzLatticePSF,[],2))';

        % real dither
        xLatticePSF = xzLatticePSF((N+1)/2,:);
        [~,locs] = findpeaks(xLatticePSF,'MinPeakHeight',0.8);
        xzLatticePSF = xzLatticePSF(:,locs(2):locs(3));
        xzLatticePSFDither = xzLatticePSFDither(:,locs(2):locs(3));
        % dither_period_pixel(j,i) = mode(diff(locs))-1;
        % xzLatticePSFDither = zeros(size(xzLatticePSF));
        % for k = 0:dither_period_pixel(j,i)
        %     xzLatticePSFDither = xzLatticePSFDither + circshift(xzLatticePSF,k,2);
        % end
        % xzLatticePSFDither = xzLatticePSFDither/max(xzLatticePSFDither,[],'all');
        
        % fig2 = figure;imagesc(xzLatticePSFDither)
        % drawnow
        % pause(1)
        % close(fig2)
        
        Duty_cycle_LLS(j,i) = sum(sum(xzLatticePSF)) ./ sum(sum(xzLatticePSFDither)); % sum over FOV
    end
    plot(NA1,Duty_cycle_LLS(j,:),'-o');
end
    hold off
    grid on
    ylim([0,1])
    xlabel("NA")
    ylabel("Duty Cycle")
    title([LatticeType ', ' ProfileType])
    legend(num2str(deltaNA1'))
%% duty cycle, gaussian
gaussianNA = 0.1:0.01:0.6;

fig1 = figure;
    hold on
for i = 1:length(gaussianNA)
    k_apertureNA = gaussianNA(i) * k_wave / n;
    gaussian_mask = zeros(N,N);
    gaussian_mask= (k_apertureNA).^2 >= abs((kx_exc.^2 + kz_exc.^2));
    gaussianPupil = exp( -(kx_exc.^2 + kz_exc.^2)/ ((k_apertureNA).^2) );
    gaussianPupil = gaussianPupil.* gaussian_mask;

    xzgaussianPSF = abs(fftshift(ifft2(ifftshift(gaussianPupil)))).^2;
    xzgaussianPSF = xzgaussianPSF./max(xzgaussianPSF,[],'all');

    xzGaussianPSFDither = meshgrid(max(xzgaussianPSF,[],2))';

    %real dither
    % xLatticePSF = xzgaussianPSF((N+1)/2,:);
    % [~,locs] = findpeaks(xLatticePSF,'MinPeakHeight',0.8);
    % xzGaussianPSFDither = zeros(size(xzgaussianPSF));
    % for k = -(N+1)/2+1:(N+1)/2-1
    %     xzGaussianPSFDither = xzGaussianPSFDither + circshift(xzgaussianPSF,k,2);
    % end
    % xzGaussianPSFDither = xzGaussianPSFDither/max(xzGaussianPSFDither,[],'all');
    % fig2 = figure;imagesc(xzLatticePSFDither)
    % drawnow
    % pause(1)
    % close(fig2)
    
    Duty_cycle_gaussian(i) = sum(sum(xzgaussianPSF)) ./ sum(sum(xzGaussianPSFDither)); % sum over FOV
end
    plot(gaussianNA,Duty_cycle_gaussian,'-o');
    hold off
    grid on
    ylim([0,1])
    xlabel("NA")
    ylabel("Duty Cycle")
    title('gaussian')
    % legend(num2str(deltaNA'))

%% Field Synthesis
NA = 0.5;
deltaNA = 0.08;
if LatticeType == "hex"
    FSWeighting1 = [1 0 0 0 0 1];
    FSWeighting2 = [0 1 0 0 1 0]; % 1.9 for V2 LLS
    FSWeighting3 = [0 0 1 1 0 0]; 
else
    FSWeighting1 = [1 0 0 0];
    FSWeighting2 = [0 1 0 1]; % 1.9 for V2 LLS
    FSWeighting3 = [0 0 1 0]; 
end
SNR = [1,3,10,40,10^50];
z_Duty_cycle_FS = [];

for k = 1:length(SNR)

    FSPupil = zeros(N,N,3);
    [FSPupil(:,:,1),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA,deltaNA, ...
    0.8,0,...
    FSWeighting1);
    
    [FSPupil(:,:,2),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA,deltaNA, ...
    0.8,0,...
    FSWeighting2);
    
    [FSPupil(:,:,3),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
    NA,deltaNA, ...
    0.8,0,...
    FSWeighting3);
    
    Beam1PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,1))))).^2;
    Beam2PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,2))))).^2;
    Beam3PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,3))))).^2;

    FSPSF = Beam1PSF + Beam2PSF + Beam3PSF;
    Beam1PSF = Beam1PSF/max(FSPSF,[],'all');
    Beam2PSF = Beam2PSF/max(FSPSF,[],'all');
    Beam3PSF = Beam3PSF/max(FSPSF,[],'all');
    FSPSF = FSPSF/max(FSPSF,[],'all');
    
    ROI = Z_exc > -ROI_range & Z_exc < ROI_range;
    zFSPSF = FSPSF(ROI,(N+1)/2);
    zBeam1PSF = Beam1PSF(ROI,(N+1)/2);
    zBeam2PSF = Beam2PSF(ROI,(N+1)/2);
    zBeam3PSF = Beam3PSF(ROI,(N+1)/2);

    zBeam1PSF(zBeam1PSF < 10^-20) = 0;
    zBeam2PSF(zBeam2PSF < 10^-20) = 0;
    zBeam3PSF(zBeam3PSF < 10^-20) = 0;

    % generate threshold 
    noise = zFSPSF./(3*SNR(k)); 
    Beam1_Duty_cycle_FS = double(zBeam1PSF-noise > 10^-20);
    Beam1_Duty_cycle_FS(zBeam1PSF==0) = 1;
    Beam2_Duty_cycle_FS = double(zBeam2PSF-noise > 10^-20);
    Beam2_Duty_cycle_FS(zBeam2PSF==0) = 1;
    Beam3_Duty_cycle_FS = double(zBeam3PSF-noise > 10^-20);
    Beam3_Duty_cycle_FS(zBeam3PSF==0) = 1;
    z_Duty_cycle_FS(:,k) = (Beam1_Duty_cycle_FS + Beam2_Duty_cycle_FS + Beam3_Duty_cycle_FS)/3;
    Dutycycle(1,k) = sum(z_Duty_cycle_FS(:,k)) ./ (sum(ROI));

    fig = figure;
    plot(Z_exc(ROI),zFSPSF)
    hold on
    plot(Z_exc(ROI),zBeam1PSF,'LineWidth',2)
    plot(Z_exc(ROI),zBeam2PSF,'LineWidth',2)
    plot(Z_exc(ROI),zBeam3PSF,'LineWidth',1)
    plot(Z_exc(ROI),noise,'LineWidth',2)
    plot(Z_exc(ROI),z_Duty_cycle_FS(:,k),'-o')
    ylim([0,1])
    xlim([-ROI_range,ROI_range])
    legend("t1+t2+t3,Incoherent","+kx,t1","kx0,t2","-kx,t3","noise","Duty Cycle(z)")
    grid on
    title("SNR=" + num2str(SNR(k)) + ", DutyCycle=" + num2str(Dutycycle(1,k)))
    xlabel("Z")
    hold off
    print(fig, '-dsvg', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
    print(fig, '-dpng', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
    % count(k,:) = [sum(z_Duty_cycle_FS(:,k)==0) sum(z_Duty_cycle_FS(:,k)==1/3) sum(z_Duty_cycle_FS(:,k)==2/3) sum(z_Duty_cycle_FS(:,k)==1)];

end

% fig = figure;
% hold on
% for j = 1:size(z_Duty_cycle_FS,2)
% end
% ylim([0,1])
% xlim([-30,30])
% xlabel("Z")
% ylabel("Duty Cycle")
% legend(num2str(SNR'))
% title([LatticeType ', ' ProfileType])
% print(fig, '-dsvg', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_DutyCycle' '.SVG'],'-r300')
% print(fig, '-dpng', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_DutyCycle' '.PNG'],'-r300')


% fig = figure;
% X = categorical(SNR);
% b = bar(X ,count./N,"stacked");
% legend("0","1/3","2/3","1")
% xlabel("SNR")
% ylabel("count (%)")
% print(fig, '-dsvg', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_Count' '.SVG'],'-r300')
% print(fig, '-dpng', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_Count' '.PNG'],'-r300')
    
%% FS, all NA

if LatticeType == "hex"
    FSWeighting1 = [1 0 0 0 0 1];
    FSWeighting2 = [0 1 0 0 1 0]; % 1.9 for V2 LLS
    FSWeighting3 = [0 0 1 1 0 0]; 
else
    FSWeighting1 = [1 0 0 0];
    FSWeighting2 = [0 1 0 1]; % 1.9 for V2 LLS
    FSWeighting3 = [0 0 1 0]; 
end
SNR = [1,3,10,40,10^50];
z_Duty_cycle_FS = [];
for k = 1:length(SNR)
    fig = figure;
    hold on
    for j = 1:length(deltaNA1)
        for i = 1:length(NA1)
            FSPupil = zeros(N,N,3);
            [FSPupil(:,:,1),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
            NA1(i),deltaNA1(j), ...
            0.8,0,...
            FSWeighting1);
            
            [FSPupil(:,:,2),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
            NA1(i),deltaNA1(j), ...
            0.8,0,...
            FSWeighting2);
            
            [FSPupil(:,:,3),~,~] = GetLatticePupil(LatticeType,ProfileType, ...
            NA1(i),deltaNA1(j), ...
            0.8,0,...
            FSWeighting3);
            
            Beam1PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,1))))).^2;
            Beam2PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,2))))).^2;
            Beam3PSF= abs(fftshift(ifft2(ifftshift(FSPupil(:,:,3))))).^2;
        
            FSPSF = Beam1PSF + Beam2PSF + Beam3PSF;
            Beam1PSF = Beam1PSF/max(FSPSF,[],'all');
            Beam2PSF = Beam2PSF/max(FSPSF,[],'all');
            Beam3PSF = Beam3PSF/max(FSPSF,[],'all');
            FSPSF = FSPSF/max(FSPSF,[],'all');
            
            ROI = Z_exc > -ROI_range & Z_exc < ROI_range;
            zFSPSF = FSPSF(ROI,(N+1)/2);
            zBeam1PSF = Beam1PSF(ROI,(N+1)/2);
            zBeam2PSF = Beam2PSF(ROI,(N+1)/2);
            zBeam3PSF = Beam3PSF(ROI,(N+1)/2);
        
            zBeam1PSF(zBeam1PSF < 10^-20) = 0;
            zBeam2PSF(zBeam2PSF < 10^-20) = 0;
            zBeam3PSF(zBeam3PSF < 10^-20) = 0;
        
            % generate threshold 
            noise = zFSPSF./(3*SNR(k)); 
            Beam1_Duty_cycle_FS = double(zBeam1PSF-noise > 10^-20);
            Beam1_Duty_cycle_FS(zBeam1PSF==0) = 1;
            Beam2_Duty_cycle_FS = double(zBeam2PSF-noise > 10^-20);
            Beam2_Duty_cycle_FS(zBeam2PSF==0) = 1;
            Beam3_Duty_cycle_FS = double(zBeam3PSF-noise > 10^-20);
            Beam3_Duty_cycle_FS(zBeam3PSF==0) = 1;

            z_Duty_cycle_FS = (Beam1_Duty_cycle_FS + Beam2_Duty_cycle_FS + Beam3_Duty_cycle_FS)/3;
            Dutycycle(j,i) = sum(z_Duty_cycle_FS) ./ (sum(ROI));
        
            % fig = figure;
            % plot(Z_exc,zFSPSF)
            % hold on
            % plot(Z_exc,zBeam1PSF,'LineWidth',2)
            % plot(Z_exc,zBeam2PSF,'LineWidth',2)
            % plot(Z_exc,zBeam3PSF,'LineWidth',1)
            % plot(Z_exc,noise,'LineWidth',2)
            % plot(Z_exc,z_Duty_cycle_FS(:,k),'-o')
            % ylim([0,1])
            % xlim([-60,60])
            % legend("Incoherent","+kx","kx0","-kx","noise","Duty Cycle(z)")
            % grid on
            % title("SNR=" + num2str(SNR(k)) + ", DutyCycle=" + num2str(Dutycycle(1,k)))
            % xlabel("Z")
            % hold off
            % print(fig, '-dsvg', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
            % print(fig, '-dpng', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
            % count(k,:) = [sum(z_Duty_cycle_FS(:,k)==0) sum(z_Duty_cycle_FS(:,k)==1/3) sum(z_Duty_cycle_FS(:,k)==2/3) sum(z_Duty_cycle_FS(:,k)==1)];
        end
        plot(NA1, Dutycycle(j,:),'-o')
    end
    hold off
    grid on
    ylim([0,1])
    xlabel("NA")
    ylabel("Duty Cycle")
    title([LatticeType ', ' ProfileType ', SNR=' num2str(SNR(k))])
    legend(num2str(deltaNA1'))   
    print(fig, '-dsvg', [ LatticeType '_' 'DutyCycle_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
    print(fig, '-dpng', [ LatticeType '_' 'DutyCycle_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
end

%% PEARLS
NA = 0.5;
deltaNA = 0.08;
SNR = [1,3,10,40,10^50];
z_Duty_cycle_SW = [];

for k = 1:length(SNR)

    [SWPupil,~,~] = GetSWPairPupil('tophat',NA,NA/2,...
                                           deltaNA,deltaNA*2,...
                                           1);

    Beam1= abs(fftshift(ifft2(ifftshift(SWPupil(:,:,1))))).^2;
    Beam2= abs(fftshift(ifft2(ifftshift(SWPupil(:,:,2))))).^2;

    IncoherentPSF = Beam1 + Beam2;
    Beam1PSF = IncoherentPSF/max(IncoherentPSF,[],'all');
    Beam2PSF = Beam1PSF;
    Beam3PSF = Beam1PSF;
    IncoherentPSF = IncoherentPSF/max(IncoherentPSF,[],'all');
    
    ROI = Z_exc > -ROI_range & Z_exc < ROI_range;
    zIncoherentPSF = IncoherentPSF(ROI,(N+1)/2);
    zBeam1PSF = Beam1PSF(ROI,(N+1)/2);
    zBeam2PSF = Beam2PSF(ROI,(N+1)/2);
    zBeam3PSF = Beam3PSF(ROI,(N+1)/2);

    zBeam1PSF(zBeam1PSF < 10^-20) = 0;
    zBeam2PSF(zBeam2PSF < 10^-20) = 0;
    zBeam3PSF(zBeam3PSF < 10^-20) = 0;

    % generate threshold 
    noise = zIncoherentPSF./(3*SNR(k)); 
    Beam1_Duty_cycle_SW = double(zBeam1PSF-noise > 10^-20);
    Beam1_Duty_cycle_SW(zBeam1PSF==0) = 1;
    Beam2_Duty_cycle_SW = double(zBeam2PSF-noise > 10^-20);
    Beam2_Duty_cycle_SW(zBeam2PSF==0) = 1;
    Beam3_Duty_cycle_SW = double(zBeam3PSF-noise > 10^-20);
    Beam3_Duty_cycle_SW(zBeam3PSF==0) = 1;

    z_Duty_cycle_SW(:,k) = (Beam1_Duty_cycle_SW + Beam2_Duty_cycle_SW + Beam3_Duty_cycle_SW)/3;
    Dutycycle(1,k) = sum(z_Duty_cycle_SW(:,k)) ./ (sum(ROI));

    fig = figure;
    plot(Z_exc(ROI),zFSPSF)
    hold on
    plot(Z_exc(ROI),zBeam1PSF,'LineWidth',2)
    plot(Z_exc(ROI),zBeam2PSF,'LineWidth',2)
    plot(Z_exc(ROI),zBeam3PSF,'LineWidth',2)
    plot(Z_exc(ROI),noise,'LineWidth',2)
    plot(Z_exc(ROI),z_Duty_cycle_SW(:,k),'-o')
    ylim([0,1])
    xlim([-ROI_range,ROI_range])
    legend("t1+t2+t3,Incoherent","t1","t2","t3","noise","Duty Cycle(z)")
    grid on
    title("SNR=" + num2str(SNR(k)) + ", DutyCycle=" + num2str(Dutycycle(1,k)))
    xlabel("Z")
    hold off
    print(fig, '-dsvg', [ 'SW_' LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
    print(fig, '-dpng', [ 'SW_' LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
    % count(k,:) = [sum(z_Duty_cycle_FS(:,k)==0) sum(z_Duty_cycle_FS(:,k)==1/3) sum(z_Duty_cycle_FS(:,k)==2/3) sum(z_Duty_cycle_FS(:,k)==1)];

end

%% PEARLS, all NA

SNR = [1,3,10,40,10^50];
z_Duty_cycle_SW = [];
for k = 1:length(SNR)
    fig = figure;
    hold on
    for j = 1:length(deltaNA1)
        NA1_array = [];
        Dutycycle = [];
        for i = 1:length(NA1)
            NA1min = NA1(i) - deltaNA1(j)/2;
            if LatticeType == "hex"
                NA2max = NA1(i)/2 + deltaNA1(j);
                NA2 = NA1(i)/2;
                deltaNA2 = deltaNA1(j)*2;
            else
                NA2max = sqrt(2 * NA1(i) * deltaNA1(j));
                deltaNA2 = NA2max;
                NA2 = deltaNA2/2;
            end

            if NA1min > NA2max
                NA1_array = [NA1_array NA1(i)];
                [SWPupil,~,~] = GetSWPairPupil('tophat',NA,NA/2,...
                                           deltaNA,deltaNA*2,...
                                           1);

    Beam1= abs(fftshift(ifft2(ifftshift(SWPupil(:,:,1))))).^2;
    Beam2= abs(fftshift(ifft2(ifftshift(SWPupil(:,:,2))))).^2;

    IncoherentPSF = Beam1 + Beam2;
    Beam1PSF = IncoherentPSF/max(IncoherentPSF,[],'all');
    Beam2PSF = Beam1PSF;
    Beam3PSF = Beam1PSF;
    IncoherentPSF = IncoherentPSF/max(IncoherentPSF,[],'all');
    
    ROI = Z_exc > -ROI_range & Z_exc < ROI_range;
    zIncoherentPSF = IncoherentPSF(ROI,(N+1)/2);
    zBeam1PSF = Beam1PSF(ROI,(N+1)/2);
    zBeam2PSF = Beam2PSF(ROI,(N+1)/2);
    zBeam3PSF = Beam3PSF(ROI,(N+1)/2);

    zBeam1PSF(zBeam1PSF < 10^-20) = 0;
    zBeam2PSF(zBeam2PSF < 10^-20) = 0;
    zBeam3PSF(zBeam3PSF < 10^-20) = 0;

    % generate threshold 
    noise = zIncoherentPSF./(3*SNR(k)); 
    Beam1_Duty_cycle_SW = double(zBeam1PSF-noise > 10^-20);
    Beam1_Duty_cycle_SW(zBeam1PSF==0) = 1;
    Beam2_Duty_cycle_SW = double(zBeam2PSF-noise > 10^-20);
    Beam2_Duty_cycle_SW(zBeam2PSF==0) = 1;
    Beam3_Duty_cycle_SW = double(zBeam3PSF-noise > 10^-20);
    Beam3_Duty_cycle_SW(zBeam3PSF==0) = 1;
          
    z_Duty_cycle_SW = (Beam1_Duty_cycle_SW + Beam2_Duty_cycle_SW + Beam3_Duty_cycle_SW)/3;
    Dutycycle = sum(z_Duty_cycle_SW) ./ (sum(ROI));

                % fig = figure;
                % plot(Z_exc,zFSPSF)
                % hold on
                % plot(Z_exc,zBeam1PSF,'LineWidth',2)
                % plot(Z_exc,zBeam2PSF,'LineWidth',2)
                % plot(Z_exc,zBeam3PSF,'LineWidth',1)
                % plot(Z_exc,noise,'LineWidth',2)
                % plot(Z_exc,z_Duty_cycle_FS(:,k),'-o')
                % ylim([0,1])
                % xlim([-60,60])
                % legend("Incoherent","+kx","kx0","-kx","noise","Duty Cycle(z)")
                % grid on
                % title("SNR=" + num2str(SNR(k)) + ", DutyCycle=" + num2str(Dutycycle(1,k)))
                % xlabel("Z")
                % hold off
                % print(fig, '-dsvg', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
                % print(fig, '-dpng', [ LatticeType '_' num2str(NA) '_' num2str(deltaNA) '_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
                % count(k,:) = [sum(z_Duty_cycle_FS(:,k)==0) sum(z_Duty_cycle_FS(:,k)==1/3) sum(z_Duty_cycle_FS(:,k)==2/3) sum(z_Duty_cycle_FS(:,k)==1)];
            end
        end
        plot(NA1_array, Dutycycle,'-o')
    end
    hold off
    grid on
    ylim([0,1])
    xlabel("NA")
    ylabel("Duty Cycle")
    title([LatticeType ', ' ProfileType ', SNR=' num2str(SNR(k))])
    legend(num2str(deltaNA1'))   
    print(fig, '-dsvg', [ 'SW_' LatticeType '_' 'DutyCycle_SNR_' num2str(SNR(k)) '.SVG'],'-r300')
    print(fig, '-dpng', [ 'SW_' LatticeType '_' 'DutyCycle_SNR_' num2str(SNR(k)) '.PNG'],'-r300')
end