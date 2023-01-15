function [LatticePupil,LatticeMask,LatticeMetaData] = GetLatticePupil(LatticeType,ProfileType,NAIdeal,deltaNA,Weighting)
    % generates pupil function of hex/square lattice
    % LatticeType = 'hex','square'
    % ProfileType = 'gaussian','tophat'
if contains(LatticeType,'hex')
    theta = [30,90,150,210,270,330];
elseif contains(LatticeType,'square')
    theta = [0,90,180,270];
else
    error("Incorrect Lattice type")
end

getParameters;
CalculatePhysics;

NAmax = NAIdeal + deltaNA/2;
NAmin = NAIdeal - deltaNA/2;
k_ideal = k_wave * NAIdeal / n;
k_deltaNA = deltaNA ./ n * k_wave;
k_NAmax = NAmax /n * k_wave; % k
k_NAmin = NAmin /n * k_wave; 

kxposition = k_ideal * cosd(theta) /deltak; % pixel
kzposition = k_ideal * sind(theta) /deltak; % pixel
deltaNApixels = k_deltaNA / deltak;

LatticePupil = zeros(N,N);

if contains(ProfileType,'gaussian')

    gaussian = Weighting .* exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA/2).^2) );

    if contains(LatticeType,'square')
        for j = 1:length(kxposition)
        LatticePupil( ...
            (N+1)/2 + round(kzposition(j)),...
            (N+1)/2 + round(kxposition(j)) ) = 1;
        end
        LatticePupil = conv2(LatticePupil,gaussian,'same');
    else
        for j = 1:length(kxposition)
            if theta(j) == 90 || theta(j) == 270
                LatticePupil( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = 1;
            end
        end
        LatticePupil_SW = conv2(LatticePupil,gaussian,'same');

        LatticePupil = zeros(N,N);
        gaussian = Weighting .* exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA).^2) );
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
                LatticePupil( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = 1;
            end
        end
        LatticePupil_Lattice = conv2(LatticePupil,gaussian,'same');
        
        LatticePupil = LatticePupil_SW + LatticePupil_Lattice;
    end

    LatticeMask = ((k_NAmax*2 > sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_NAmin/2 < sqrt(kx_exc.^2 + kz_exc.^2)));

elseif contains(ProfileType,'tophat')
    
    for j = 1:length(kxposition)
    LatticePupil( ...
        (N+1)/2 + round(kzposition(j) - deltaNApixels-20) : (N+1)/2 + round(kzposition(j)+ deltaNApixels+20),...
        (N+1)/2 + round(kxposition(j)) ) = Weighting;
    end

    LatticeMask = ((k_NAmax > sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_NAmin < sqrt(kx_exc.^2 + kz_exc.^2)));
else
    error("Incorrect Intensity Profile")
end

LatticePupil = LatticePupil .* LatticeMask .* k_wave./ky_exc;
LatticePupil(LatticePupil == inf) = 0;
LatticePupil = fillmissing(LatticePupil,'constant',0);

LatticeMetaData.NA = NAIdeal;
LatticeMetaData.deltaNA = deltaNA;
LatticeMetaData.NAmax = NAmax;
LatticeMetaData.NAmin = NAmin;
LatticeMetaData.NAWeighting = Weighting;
