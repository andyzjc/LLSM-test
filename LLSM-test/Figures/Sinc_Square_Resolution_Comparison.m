clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;

%% Approximation 
NA1 = 0.5;
deltaNA1 = 0.02:0.005:0.17; 
sinc_deltaNA2 = 0.16:0.005:0.4; 
weighting = 1;
% figure;
% scatter(deltaNA1,deltaNA2)
% grid on
% xlim([0,0.5])
% ylim([0,0.5])
% hold on

%% Numeric way

Square_yFWHM = [];
for i = 1:length(deltaNA1)
    Square_deltaNA2 = sqrt(2*NA1*deltaNA1(i));
    Square_NA2 = Square_deltaNA2/2;
        [SWPupil,~,~] = GetSWPairPupil('tophat',NA1,Square_NA2,...
                                                deltaNA1(i),Square_deltaNA2,...
                                                weighting);
        [~,PSFIncoherent,~] = SimulateSWPair(SWPupil);
        PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
        yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
        SWyFWHM1 = find(yindex,1,'first') ;
        Square_yFWHM(i) = abs(Y_exc(SWyFWHM1)).*2;
end

Sinc_yFWHM = [];
for i = 1:length(sinc_deltaNA2)
    k_apertureNA = sinc_deltaNA2(i) * k_wave / n;
    sincPupil = zeros(N,N);
    sincPupil(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
    [~,sincPSFDithered,~] = SimulateLattice(sincPupil);
    yindex = 1-(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:))));
    SincyFWHM1index = find(yindex,1,'first');
    Sinc_yFWHM(i) = abs(Y_exc(SincyFWHM1index)).*2;
end

figure
plot(deltaNA1,Square_yFWHM)
hold on
plot(sinc_deltaNA2,Sinc_yFWHM)
grid on
legend("Square","sinc")

beam1minyFWHM = min(Square_yFWHM);
 beam1maxyFWHM = max(Square_yFWHM);
 beam2minyFWHM = min(Sinc_yFWHM);
 beam2maxyFWHM = max(Sinc_yFWHM);

 if beam1minyFWHM >= beam2minyFWHM
    minyFWHM_equal = beam1minyFWHM;
 else
    minyFWHM_equal = beam2minyFWHM;
 end

 if beam1maxyFWHM >= beam2maxyFWHM
    maxyFWHM_equal = beam2maxyFWHM;
 else
    maxyFWHM_equal = beam1maxyFWHM;
 end

yFWHM_equal = (minyFWHM_equal:0.1*wavelength_exc:maxyFWHM_equal)';

% find intersection
deltaNA_reference = NaN(size(yFWHM_equal));
deltaNA_equivalence = deltaNA_reference;
%%
for jj = 2:length(yFWHM_equal)
    [x1,~,~,~] = intersections(deltaNA1,Square_yFWHM,...
        0:0.01:0.4,yFWHM_equal(jj)*ones(1,length(0:0.01:0.4)),0);
    [x2,~,~,~] = intersections(sinc_deltaNA2,Sinc_yFWHM,...
        0:0.01:0.4,yFWHM_equal(jj)*ones(1,length(0:0.01:0.4)),0);
    deltaNA_reference(jj,1) = x1(1);
    deltaNA_equivalence(jj,1) = x2(1);
end
 
plot(deltaNA_reference,deltaNA_equivalence)
xlabel("deltaNA Sqaure")
ylabel("deltaNA sinc")
grid on


%% Strict Condition

%% Analytical curve
% x = 0.01:0.01:2;
% y = 1/2*(sqrt(x)+sqrt(1./x));
% plot(x,y)
% grid on
% xlim([0,2])
% % ylim([0,1])
% hold on