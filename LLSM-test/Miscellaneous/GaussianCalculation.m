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
GaussianPupil06 = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p6.mat');
GaussianPupil06 = GaussianPupil06.GaussianPupil;
GaussianPupil04 = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p4.mat');
GaussianPupil04 = GaussianPupil04.GaussianPupil;
GaussianPupil03 = load('LLSM-test/Aberration/Zernike/1dGaussianPupil513_NA_0p3.mat');
GaussianPupil03 = GaussianPupil03.GaussianPupil;

Airy06 = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p6.mat');
AiryPupil06 = Airy06.AiryPupil;
Airy04 = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p4.mat');
AiryPupil04 = Airy04.AiryPupil;
Airy03 = load('LLSM-test/Aberration/Zernike/AiryPupil513_0p3.mat');
AiryPupil03 = Airy03.AiryPupil;
Airy024 = load('LLSM-test/Aberration/Zernike/sincPupil513_0p24.mat');
AiryPupil024 = Airy024.sincPupil;

GaussianPupil04 = AiryPupil04 .* GaussianPupil04;
GaussianPupil03 = AiryPupil03 .* GaussianPupil03;

firemap = fire(256);

fig1 = figure;
plot(KZ_exc,GaussianPupil06(:,(N+1)/2))
hold on
plot(KZ_exc,GaussianPupil04(:,(N+1)/2))
plot(KZ_exc,GaussianPupil03(:,(N+1)/2))
axis square
legend("0.6","0.4","0.3")

%%
[GaussianPSF03,~,~] = SimulateLattice(GaussianPupil03);
% GaussianPSF = PSFIncoherent/SWcenter(2,1);

GaussianPSF03 = GaussianPSF03/max(GaussianPSF03,[],'all');
yGaussian03 = squeeze(GaussianPSF03((N+1)/2,(N+1)/2,:));
zGaussian03 = squeeze(GaussianPSF03(:,(N+1)/2,(N+1)/2));

[GaussianPSF04,~,~] = SimulateLattice(GaussianPupil04);
% GaussianPSF = PSFIncoherent/SWcenter(2,1);

GaussianPSF04 = GaussianPSF04/max(GaussianPSF04,[],'all');
yGaussian04 = squeeze(GaussianPSF04((N+1)/2,(N+1)/2,:));
zGaussian04 = squeeze(GaussianPSF04(:,(N+1)/2,(N+1)/2));

[GaussianPSF06,~,~] = SimulateLattice(GaussianPupil06);
% GaussianPSF = PSFIncoherent/SWcenter(2,1);

GaussianPSF06 = GaussianPSF06/max(GaussianPSF06,[],'all');
yGaussian06 = squeeze(GaussianPSF06((N+1)/2,(N+1)/2,:));
zGaussian06 = squeeze(GaussianPSF06(:,(N+1)/2,(N+1)/2));

fig1 = figure;
plot(Z_exc,zGaussian06)
hold on
plot(Z_exc,zGaussian04)
plot(Z_exc,zGaussian03)
axis square
legend("0.6","0.4","0.3")

GaussianSumI06 = [0; cumsum(zGaussian06((N+1)/2:end))];
GaussianSumI06 = GaussianSumI06/max(GaussianSumI06,[],'all');

GaussianSumI04 = [0; cumsum(zGaussian04((N+1)/2:end))];
GaussianSumI04 = GaussianSumI04/max(GaussianSumI04,[],'all');

GaussianSumI03 = [0; cumsum(zGaussian03((N+1)/2:end))];
GaussianSumI03 = GaussianSumI03/max(GaussianSumI03,[],'all');

fig = figure;
plot(Z_exc((N+1)/2:end),GaussianSumI06(1:end-1))
hold on
plot(Z_exc((N+1)/2:end),GaussianSumI04(1:end-1))
plot(Z_exc((N+1)/2:end),GaussianSumI03(1:end-1))
legend("0.6","0.4","0.3")
xlim([0,5])
ylim([0,1])