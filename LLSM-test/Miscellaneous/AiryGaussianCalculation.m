clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

% detection
PSFdet = getDetectionPSF;
PSFdet = PSFdet./(max(max(max(PSFdet))));

%%
GaussianPupil = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p6.mat');
GaussianPupil = GaussianPupil.GaussianPupil;

Airy024 = load('LLSM-test/Aberration/Zernike/sincPupil513_0p24.mat');
AiryPupil024 = Airy024.sincPupil;

AiryPupil = AiryPupil024;

firemap = fire(256);

Phase = GetSingleZmodePupil(2,-2,4);
AberratedAiryPupil = AiryPupil .* Phase;
AberratedGaussianPupil = GaussianPupil .* Phase;

fig1 = figure;
subplot(1,3,1)
imagesc(KX_exc,KZ_exc,AiryPupil)
axis image
subplot(1,3,2)
imagesc(KX_exc,KZ_exc,GaussianPupil)
axis image
subplot(1,3,3)
plot(KZ_exc,GaussianPupil(:,(N+1)/2))
hold on
plot(KZ_exc,AiryPupil(:,(N+1)/2))
axis square
subplot(1,3,3)
legend("Gaussian","Airy")

%% 

[GaussianPSF,GaussianPSFDithered,GaussianPSFCenter] = SimulateLattice(GaussianPupil);
% GaussianPSF = PSFIncoherent/SWcenter(2,1);

GaussianPSF = GaussianPSF/max(GaussianPSF,[],'all');
GaussianOTF = real(fftshift(ifft2(ifftshift(GaussianPSF(:,:,(N+1)/2)))));
GaussianOTF = GaussianOTF/max(GaussianOTF,[],'all');
yGaussian = squeeze(GaussianPSF((N+1)/2,(N+1)/2,:));
zGaussian = squeeze(GaussianPSF(:,(N+1)/2,(N+1)/2));
GaussianOverall = GaussianPSF.* PSFdet;

[AiryPSF,AiryPSFDithered,AiryPSFCenter] = SimulateLattice(AiryPupil);
% AiryPSF = PSFCoherent/SWcenter(2,1);

AiryPSF = AiryPSF/max(AiryPSF,[],'all');
AiryOTF = real(fftshift(ifft2(ifftshift(AiryPSF(:,:,(N+1)/2)))));
AiryOTF = AiryOTF/max(AiryOTF,[],'all');
yAiry = squeeze(AiryPSF((N+1)/2,(N+1)/2,:));
zAiry = squeeze(AiryPSF(:,(N+1)/2,(N+1)/2));
AiryOverall = AiryPSF .* PSFdet;

[AberratedGaussianPSF,~,~] = SimulateLattice(AberratedGaussianPupil);
AberratedGaussianPSF = AberratedGaussianPSF/max(AberratedGaussianPSF,[],'all');
AberratedGaussianOverall = AberratedGaussianPSF.* PSFdet;

[AberratedAiryPSF,~,~] = SimulateLattice(AberratedAiryPupil);
AberratedAiryPSF = AberratedAiryPSF/max(AberratedAiryPSF,[],'all');
AberratedAiryOverall = AberratedAiryPSF .* PSFdet;

fig2 = figure;
subplot(3,2,1)
imagesc(X_exc,Z_exc,squeeze(AiryPSF(:,:,(N+1)/2)));
title("Airy, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])

subplot(3,2,2)
imagesc(X_exc,Z_exc,squeeze(GaussianPSF(:,:,(N+1)/2)));
title("Gaussian, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])

subplot(3,2,3)
imagesc(KX_exc,KZ_exc,AiryOTF);
hold on
plot(AiryOTF(:,(N+1)/2)-0.5,KZ_exc)
title("Airy real(zOTF), K_Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("k_x/(4\pin/\lambda_{exc})")
ylabel("k_z/(4\pin/\lambda_{exc})")
xlim([-0.5,0.5])
ylim([-0.5,0.5])

subplot(3,2,4)
imagesc(KX_exc,KZ_exc,GaussianOTF);
hold on
plot(GaussianOTF(:,(N+1)/2)-0.5,KZ_exc)
title("Gaussian real(zOTF), K_Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("k_x/(4\pin/\lambda_{exc})")
ylabel("k_z/(4\pin/\lambda_{exc})")
xlim([-0.5,0.5])
ylim([-0.5,0.5])

[~,maxindex] = max(yAiry);
index = 1-(yAiry) <= 0.5*max(yAiry);
if ~isempty(index)
yFWHM1 = Y_exc(find(index,1,'first')) ;
yFWHM2 = Y_exc(find(index,1,'last'));
if abs(yFWHM1) == abs(yFWHM2)
yFWHM = abs(yFWHM1) + abs(yFWHM2);
elseif abs(Y_exc(maxindex) - yFWHM1) > abs(Y_exc(maxindex) - yFWHM2)
yFWHM = abs(Y_exc(maxindex) - yFWHM1)*2;
else
yFWHM = abs(Y_exc(maxindex) - yFWHM2)*2;
end
else
yFWHM = "N/A";
end
AiryyFWHM = yFWHM;

[~,maxindex] = max(zAiry);
index = 1-(zAiry) <= 0.5*max(zAiry);
if ~isempty(index)
zFWHM1 = Z_exc(find(index,1,'first')) ;
zFWHM2 = Z_exc(find(index,1,'last'));
if abs(zFWHM1) == abs(zFWHM2)
zFWHM = abs(zFWHM1) + abs(zFWHM2);
elseif abs(Z_exc(maxindex) - zFWHM1) > abs(Z_exc(maxindex) - zFWHM2)
zFWHM = abs(Z_exc(maxindex) - zFWHM1)*2;
else
zFWHM = abs(Z_exc(maxindex) - zFWHM2)*2;
end
else
zFWHM = "N/A";
end
AiryzFWHM = zFWHM;

subplot(3,2,5)
imagesc(Y_exc,Z_exc,squeeze(AiryPSF(:,(N+1)/2,:)));
title("Airy, X=0, yFWHM=" + num2str(AiryyFWHM) + ", zFWHM=" + num2str(AiryzFWHM))
colorbar;
axis image;
colormap(firemap)
xlabel("y/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
ylim([-20,20])
xlim([-60,60])


[~,maxindex] = max(yGaussian);
index = 1-(yGaussian) <= 0.5*max(yGaussian);
if ~isempty(index)
yFWHM1 = Y_exc(find(index,1,'first')) ;
yFWHM2 = Y_exc(find(index,1,'last'));
if abs(yFWHM1) == abs(yFWHM2)
yFWHM = abs(yFWHM1) + abs(yFWHM2);
elseif abs(Y_exc(maxindex) - yFWHM1) > abs(Y_exc(maxindex) - yFWHM2)
yFWHM = abs(Y_exc(maxindex) - yFWHM1)*2;
else
yFWHM = abs(Y_exc(maxindex) - yFWHM2)*2;
end
else
yFWHM = "N/A";
end
GaussianyFWHM = yFWHM;

[~,maxindex] = max(zGaussian);
index = 1-(zGaussian) <= 0.5*max(zGaussian);
if ~isempty(index)
zFWHM1 = Z_exc(find(index,1,'first')) ;
zFWHM2 = Z_exc(find(index,1,'last'));
if abs(zFWHM1) == abs(zFWHM2)
zFWHM = abs(zFWHM1) + abs(zFWHM2);
elseif abs(Z_exc(maxindex) - zFWHM1) > abs(Z_exc(maxindex) - zFWHM2)
zFWHM = abs(Z_exc(maxindex) - zFWHM1)*2;
else
zFWHM = abs(Z_exc(maxindex) - zFWHM2)*2;
end
else
zFWHM = "N/A";
end
GaussianzFWHM = zFWHM;


subplot(3,2,6)
imagesc(Y_exc,Z_exc,squeeze(GaussianPSF(:,(N+1)/2,:)));
title("Gaussian, X=0, yFWHM=" + num2str(GaussianyFWHM) + ", zFWHM=" + num2str(GaussianzFWHM))
colorbar;
axis image;
colormap(firemap)
xlabel("y/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
ylim([-20,20])
xlim([-60,60])

fig = figure;
subplot(1,2,2)
imagesc(Y_exc,X_exc,squeeze(GaussianPSF((N+1)/2,:,:)));
title("Gaussian, Z=0")
colorbar;
axis image;
colormap(firemap)
xlabel("y/(\lambda_{exc}/n)")
ylabel("x/(\lambda_{exc}/n)")
ylim([-20,20])
xlim([-60,60])

subplot(1,2,1)
imagesc(Y_exc,X_exc,squeeze(AiryPSF((N+1)/2,:,:)));
title("Airy, Z=0")
colorbar;
axis image;
colormap(firemap)
xlabel("y/(\lambda_{exc}/n)")
ylabel("x/(\lambda_{exc}/n)")
ylim([-20,20])
xlim([-60,60])

fig = figure;
plot(KZ_exc,AiryOTF(:,(N+1)/2))
hold on
plot(KZ_exc,GaussianOTF(:,(N+1)/2))
legend("Airy","Gaussian")
xlim([-0.5,0.5])
xlabel("z/(\lambda_{exc}/n)")
grid on
title("zOTF")

fig = figure; 
subplot(2,2,1)
imagesc(X_exc,Z_exc,squeeze(AiryPSF(:,:,(N+1)/2)));
title("Airy, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])

subplot(2,2,2)
imagesc(X_exc,Z_exc,squeeze(GaussianPSF(:,:,(N+1)/2)));
title("Gaussian, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])

subplot(2,2,3)
imagesc(X_exc,Z_exc,squeeze(AberratedAiryPSF(:,:,(N+1)/2)));
title("Aberrated Airy, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])

subplot(2,2,4)
imagesc(X_exc,Z_exc,squeeze(AberratedGaussianPSF(:,:,(N+1)/2)));
title("Aberrated Gaussian, Y=0")
colorbar;
axis image;
colormap(firemap)
xlabel("x/(\lambda_{exc}/n)")
ylabel("z/(\lambda_{exc}/n)")
xlim([-10,10])
ylim([-10,10])


GaussianSumI = cumsum(zGaussian((N+1)/2:end));
GaussianSumI = GaussianSumI/max(GaussianSumI,[],'all');

AirySumI = cumsum(zAiry((N+1)/2:end));
AirySumI = AirySumI/max(AirySumI,[],'all');

fig = figure;
plot(X_exc((N+1)/2:end),GaussianSumI)
hold on
plot(X_exc((N+1)/2:end),AirySumI)
legend("0.4","0.3")




