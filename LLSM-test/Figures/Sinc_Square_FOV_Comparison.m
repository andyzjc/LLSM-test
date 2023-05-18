clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

% get predefine variables 
getParameters; %modify image parameter here
CalculatePhysics;
%% Get analytical curve 
x = 0.01:0.01:2;
y = 1/4*(sqrt(x)+sqrt(1./x)).^2;
plot(x,y)
grid on
xlim([0,2])
% ylim([0,1])
hold on

%%
clc
NA1 = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7];
deltaNA1 = 0.04:0.02:0.2;
weighting = 1;

NA1Array = [];
deltaNA1Array = [];
NA2Array = [];
deltaNA2Array = [];
NA1_diff_sum = [];
yFWHM_Ratio = [];

counter = 1;
for i = 1:length(NA1)
    for j = 1:length(deltaNA1)
        deltaNA2 = sqrt(2*NA1(i)*deltaNA1(j));
        NA1min = NA1(i) - deltaNA1(j)/2;
        NA1max = NA1(i) + deltaNA1(j)/2;
        if deltaNA2 < NA1min
            NA1Array(counter,1) = NA1(i)
            deltaNA1Array(counter,1) = deltaNA1(j)
            NA2Array(counter,1) = deltaNA2/2;
            deltaNA2Array(counter,1) = deltaNA2;
            NA1_diff_sum(counter,1) = (NA1max-NA1min)/(NA1max+NA1min);

            [SWPupil,~,~] = GetSWPairPupil('tophat',NA1(i),deltaNA2/2,...
                                                    deltaNA1(j),deltaNA2,...
                                                    weighting);
            [~,PSFIncoherent,~] = SimulateSWPair(SWPupil);
            PSFIncoherent = PSFIncoherent/max(PSFIncoherent,[],'all');
            yindex = 1-(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(PSFIncoherent((N+1)/2,(N+1)/2,:))));
            SWyFWHM1 = find(yindex,1,'first') ;
            SWyFWHM2 = find(yindex,1,'last');

            k_apertureNA = NA1max * k_wave / n;
            sincPupil = zeros(N,N);
            sincPupil(:,(N+1)/2) = (k_apertureNA) >= abs(kz_exc(:,1));
            [~,sincPSFDithered,~] = SimulateLattice(sincPupil);
            sincPSFDithered = sincPSFDithered/max(sincPSFDithered,[],'all');
            yindex = 1-(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(sincPSFDithered((N+1)/2,(N+1)/2,:))));
            SincyFWHM1 = find(yindex,1,'first');
            SincyFWHM2 = find(yindex,1,'last');

            yFWHM_Ratio(counter,1) = Y_exc(SWyFWHM1)/Y_exc(SincyFWHM1);
            counter = counter + 1;
        end
    end
end

scatter(NA1_diff_sum,yFWHM_Ratio)
