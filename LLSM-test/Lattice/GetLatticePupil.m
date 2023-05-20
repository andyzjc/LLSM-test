function [LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,ProfileType,NAIdeal,deltaNA,MaskNAmax,MaskNAmin,weighting)
    % generates pupil function of hex/square lattice
    % LatticeType = 'hex','square'
    % ProfileType = 'gaussian','tophat'

if contains(LatticeType,'hex')
       theta = [30,90,150,210,270,330];
%       theta = [30,150,210,330];
%          theta = [90, 270];
elseif contains(LatticeType,'square')
     theta = [0,90,180,270];
%           theta = [0,180];
%          theta = [90, 270];
else
    error("Incorrect Lattice type")
end

getParameters;
CalculatePhysics;

NAmax = NAIdeal + deltaNA/2;
NAmin = NAIdeal - deltaNA/2;
k_ideal = k_wave * NAIdeal / n;
k_deltaNA = deltaNA / n * k_wave;
k_NAmax = NAmax /n * k_wave; % k
k_NAmin = NAmin /n * k_wave; 
k_MaskNAmax = MaskNAmax /n * k_wave; % k
k_MaskNAmin = MaskNAmin /n * k_wave; 

kxposition = k_ideal * cosd(theta) /deltak; % pixel
kzposition = k_ideal * sind(theta) /deltak; % pixel
deltaNApixels = k_deltaNA / deltak;

LatticePupil_SW = zeros(N,N);
LatticePupil_Lattice = zeros(N,N);
LatticeMask = ((k_MaskNAmax >= sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_MaskNAmin <= sqrt(kx_exc.^2 + kz_exc.^2)));
if contains(ProfileType,'gaussian')
    gaussian_SW = exp( -(kz_exc(:,1).^2)./ ( (k_deltaNA/2).^2) ); %deltaNA/2 at 1/e 
    for j = 1:length(kxposition)
        if theta(j) == 90 || theta(j) == 270
            LatticePupil_SW( ...
                (N+1)/2 + round(kzposition(j)),...
                (N+1)/2 + round(kxposition(j)) ) = weighting(j);
        end
    end
    LatticePupil_SW = conv2(LatticePupil_SW,gaussian_SW,'same');

    LatticePupil = zeros(N,N);
    if contains(LatticeType,'hex')
        gaussian_lattice = exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA).^2) );
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
                LatticePupil_Lattice( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = weighting(j);
            end
        end
        LatticePupil_Lattice = conv2(LatticePupil_Lattice,gaussian_lattice,'same');
%         LatticePupil =  LatticePupil .* LatticeMask;
    else
        deltaNA_Square = sqrt(2*deltaNA*NAIdeal);
        k_deltaNA_Square = deltaNA_Square ./ n * k_wave;
        gaussian_lattice = exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA_Square).^2) );
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
                LatticePupil_Lattice( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = weighting(j);
            end
        end
        LatticePupil_Lattice = conv2(LatticePupil_Lattice,gaussian_lattice,'same');
%         LatticePupil = LatticePupil_SW.*LatticeMask + LatticePupil_Lattice;
    end 
elseif contains(ProfileType,'tophat') %tophat doesnt need mask -> beam width define by deltaNA, following the non-diffractive condition
    
    for j = 1:length(kxposition)
        if theta(j) == 90 || theta(j) == 270
        LatticePupil_SW( ...
            (N+1)/2 + round(kzposition(j) - deltaNApixels/2) : (N+1)/2 + round(kzposition(j)+ deltaNApixels/2),...
            (N+1)/2 + round(kxposition(j)) ) = weighting(j);
        end
    end

    if contains(LatticeType,'square')
        deltaNA_Square = sqrt(2*deltaNA*NAIdeal);
        k_deltaNA_Square = deltaNA_Square ./ n * k_wave;
        deltaNApixels = k_deltaNA_Square / deltak;
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
            LatticePupil_Lattice( ...
                (N+1)/2 + round(kzposition(j) - deltaNApixels) : (N+1)/2 + round(kzposition(j)+ deltaNApixels),...
                (N+1)/2 + round(kxposition(j)) ) = weighting(j);
            end
        end
    else
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
            LatticePupil_Lattice( ...
                (N+1)/2 + round(kzposition(j) - deltaNApixels) : (N+1)/2 + round(kzposition(j)+ deltaNApixels),...
                (N+1)/2 + round(kxposition(j)) ) = weighting(j);
            end
        end
    end
    LatticeMask = zeros(N,N);
else
    error("Incorrect Intensity Profile")
end

% apply sigmoid bounding to x and z sample based on pupil

% temp = zeros(1,N);
% temp(1:(N+1)/2-1) = 1./(1+exp(-(-(N+1)/4+1:(N+1)/4-1)));
% temp((N+1)/2:end) = flip(1./(1+exp(-(-(N+1)/4:(N+1)/4-1))));
% [boundx, ~] = meshgrid(temp);
% Sample = fftshift(fft2(ifftshift(LatticePupil_SW)));
% Sample = Sample.*boundx;
% LatticePupil_SW = fftshift(ifft2(ifftshift((Sample))));
% LatticePupil_SW = LatticePupil_SW/max(LatticePupil_SW,[],'all');
% LatticePupil_SW(:,(N+1)/2) = WeightingRatio.* LatticePupil_SW(:,(N+1)/2);
% 
% Sample = fftshift(fft2(ifftshift(LatticePupil_Lattice)));
% Sample = Sample.*boundx;
% LatticePupil_Lattice = fftshift(ifft2(ifftshift((Sample))));
% LatticePupil_Lattice = LatticePupil_Lattice/max(LatticePupil_Lattice,[],'all');

LatticePupil = LatticePupil_SW + LatticePupil_Lattice;
LatticePupil = LatticePupil/max(LatticePupil,[],'all');

LatticePupil(LatticePupil == inf) = 0;
LatticePupil = fillmissing(LatticePupil,'constant',0);

LatticeMetaData.NA = NAIdeal;
LatticeMetaData.deltaNA = deltaNA;
LatticeMetaData.NAmax = NAmax;
LatticeMetaData.NAmin = NAmin;
LatticeMetaData.MaskNAmax = MaskNAmax;
LatticeMetaData.MaskNAmin = MaskNAmin;
LatticeMetaData.NAWeighting = weighting;

