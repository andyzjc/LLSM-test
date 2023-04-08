function [FSPupil,FSMask,FSMetaData] = GetFSPupil(LatticeType,ProfileType,NAIdeal,deltaNA,MaskNAmax,MaskNAmin,WeightingRatio)
    % generates pupil function of hex/square lattice
    % LatticeType = 'hex','square'
    % ProfileType = 'gaussian','tophat'
if contains(LatticeType,'hex')
     theta = [30,90,150,210,270,330];
%       theta = [30,150,210,330];
%        theta = [90, 270];
elseif contains(LatticeType,'square')
     theta = [0,90,180,270];
%         theta = [0,180];
%        theta = [90, 270];
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

LatticePupil = zeros(N,N);
if contains(ProfileType,'gaussian')
    gaussian_SW = 1 .* exp( -(kz_exc(:,1).^2)./ ( (k_deltaNA/2).^2) ); %deltaNA/2 at 1/e 
    for j = 1:length(kxposition)
        if theta(j) == 90 || theta(j) == 270
            LatticePupil( ...
                (N+1)/2 + round(kzposition(j)),...
                (N+1)/2 + round(kxposition(j)) ) = 1;
        end
    end
    LatticePupil_SW = conv2(LatticePupil,gaussian_SW,'same');

    LatticePupil = zeros(N,N);
    if contains(LatticeType,'hex')
        gaussian_lattice = (1 ./ WeightingRatio) .* exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA).^2) );
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
                LatticePupil( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = 1;
            end
        end
        LatticePupil_Lattice = conv2(LatticePupil,gaussian_lattice,'same');
        LatticePupil = LatticePupil_SW + LatticePupil_Lattice;
        LatticeMask = ((k_MaskNAmax >= sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_MaskNAmin <= sqrt(kx_exc.^2 + kz_exc.^2)));
%         LatticePupil =  LatticePupil .* LatticeMask.*k_wave./ky_exc;
        LatticePupil =  LatticePupil .* LatticeMask;
    else
        deltaNA_Square = sqrt(2*deltaNA*NAIdeal);
        k_deltaNA_Square = deltaNA_Square ./ n * k_wave;
        gaussian_lattice = (1 ./ WeightingRatio) .* exp( -(kz_exc(:,1).^2)/ ( (k_deltaNA_Square).^2) );
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
                LatticePupil( ...
                    (N+1)/2 + round(kzposition(j)),...
                    (N+1)/2 + round(kxposition(j)) ) = 1;
            end
        end
        LatticePupil_Lattice = conv2(LatticePupil,gaussian_lattice,'same');
        LatticePupil = LatticePupil_SW + LatticePupil_Lattice;
        LatticeMask = ((k_MaskNAmax >= sqrt(kx_exc.^2 + kz_exc.^2)) .* (k_MaskNAmin <= sqrt(kx_exc.^2 + kz_exc.^2)));
%         LatticePupil =  LatticePupil.*LatticeMask.*k_wave./ky_exc;
        LatticePupil =  LatticePupil.*LatticeMask;
    end 
elseif contains(ProfileType,'tophat') %tophat doesnt need mask -> beam width define by deltaNA, following the non-diffractive condition
    
    for j = 1:length(kxposition)
        if theta(j) == 90 || theta(j) == 270
        LatticePupil( ...
            (N+1)/2 + round(kzposition(j) - deltaNApixels/2) : (N+1)/2 + round(kzposition(j)+ deltaNApixels/2),...
            (N+1)/2 + round(kxposition(j)) ) = 1;
        end
    end

    if contains(LatticeType,'square')
        deltaNA_Square = sqrt(2*deltaNA*NAIdeal);
        k_deltaNA_Square = deltaNA_Square ./ n * k_wave;
        deltaNApixels = k_deltaNA_Square / deltak;
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
            LatticePupil( ...
                (N+1)/2 + round(kzposition(j) - deltaNApixels) : (N+1)/2 + round(kzposition(j)+ deltaNApixels),...
                (N+1)/2 + round(kxposition(j)) ) = (1 ./ WeightingRatio);
            end
        end
    else
        for j = 1:length(kxposition)
            if theta(j) ~= 90 && theta(j) ~= 270
            LatticePupil( ...
                (N+1)/2 + round(kzposition(j) - deltaNApixels) : (N+1)/2 + round(kzposition(j)+ deltaNApixels),...
                (N+1)/2 + round(kxposition(j)) ) = (1 ./ WeightingRatio);
            end
        end
        LatticePupil = LatticePupil .*k_wave./ky_exc;
    end
    LatticeMask = zeros(N,N);
else
    error("Incorrect Intensity Profile")
end

LatticePupil(LatticePupil == inf) = 0;
LatticePupil = fillmissing(LatticePupil,'constant',0);

% convert to non-coherent (brutally) 
counter = 1;
for i = 1:length(kxposition)/2
    temp = zeros(N,N);
    temp(:,(N+1)/2 + round(kxposition(i))-1:(N+1)/2 + round(kxposition(i))+1) = ...
        LatticePupil(:,(N+1)/2 + round(kxposition(i))-1:(N+1)/2 + round(kxposition(i))+1);
    FSPupil(:,:,counter) = temp;
    counter = counter + 1;
end
FSMask = LatticeMask;

FSMetaData.NA = NAIdeal;
FSMetaData.deltaNA = deltaNA;
FSMetaData.NAmax = NAmax;
FSMetaData.NAmin = NAmin;
FSMetaData.MaskNAmax = MaskNAmax;
FSMetaData.MaskNAmin = MaskNAmin;
FSMetaData.NAWeighting = WeightingRatio;
