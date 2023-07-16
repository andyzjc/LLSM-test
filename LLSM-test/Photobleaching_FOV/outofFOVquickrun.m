clear all
close all

addpath([pwd '/' addpath(genpath("LLSM-test/"))])

getParameters; %modify image parameter here
CalculatePhysics;

%% Simulation of out of FOV photobleaching
w = [50 100 250 500] ./ wavelength_exc; % lambda/n
NA = 0.3:0.01:0.6;
FWHM = 10:5:150; % lambda/n
[NAgrid,FWHMgrid] = meshgrid(NA,FWHM);
faceColorArray = [0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];
figure
for i = 1:length(w)
    Ratio = [];
    Ratio = 3./( 2 .* (w(i) ./ (NAgrid./n .* FWHMgrid) -1) );
    surf(NAgrid,FWHMgrid,Ratio,'FaceColor',faceColorArray(i,:))
    hold on
end
xlabel("NA")
ylabel("FWHM(\lambda_{exc}/n)")
legend(num2str(w' .* wavelength_exc) + "\lambda/n")
xlim([0.3,0.6])
hold off

%% Isoplanic patch 
w = 50; %lambda/n
FWHM = 50; %lambda/n
NA = 0.3:0.01:0.6;
figure
Ratio = [];
Ratio = 3./( 2 .* (w ./ (NAgrid./n .* FWHM) -1) );
plot(NA,Ratio)
xlabel("NA")
ylabel("Ratio")
grid on
title("w=" + num2str(w) + "\lambda/n, FHWM=" + num2str(FWHM) + " \lambda/n")

%% Simulation 
% w = [100 250 500 1000] ./ wavelength_exc; % lambda/n
% NA_simulation = 0.3:0.05:0.6;
% deltaNA_simulation = 0.02:0.02:0.14;
% [NAgrid,deltaNAgrid] = meshgrid(NA_simulation,deltaNA_simulation);
% LatticeType = 'hex';
% ProfileType = 'tophat';
% 
% % Approximated Analytical solution of FWHM 
% Ratio = [];
% for k = 1:length(w)
%     for i = 1:length(NA_simulation)
%         for j = 1:length(deltaNA_simulation)
%             if LatticeType == "hex"
%                 NAmax = NA_simulation(i)/2 + deltaNA_simulation(j);
%                 NAmin = NA_simulation(i) - deltaNA_simulation(j)/2;
%             else
%                 NAmax = sqrt(2 * NA_simulation(i) * deltaNA_simulation(j));
%                 NAmin = -NAmax;
%             end
%             [LatticePupil,~,~] = GetLatticePupil(LatticeType,ProfileType, ...
%                                                              NA_simulation(i),deltaNA_simulation(j), ...
%                                                              0.8,0.0,...
%                                                              [1 1 1 1 1 1]);
%             [~,LatticePSFDithered,~] = SimulateLattice(LatticePupil);
%             LatticePSFDithered = LatticePSFDithered/max(LatticePSFDithered,[],'all');
%             yindex = 1-(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:)) <= 0.5*max(squeeze(LatticePSFDithered((N+1)/2,(N+1)/2,:))));
%             LatticeyFWHM1 = find(yindex,1,'first');
%             yFWHM = abs(Y_exc(LatticeyFWHM1)).*2;
%             Ratio(i,j,k) = 3./( 2 .* (w(k) ./ (tand(asind(NAmax./n)) .* yFWHM) -1) );
%         end
%     end
%     hold on
%     surf(NAgrid,deltaNAgrid,Ratio(:,:,k))
% end