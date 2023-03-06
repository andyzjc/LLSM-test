function [Vol1,Vol2,GTallbeads] = beadsSimulation3d_FWHM(PSFexc,PSFdet,SNR,FWHMIndex)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % determine random beads image size and generate random beads
    beadPercent = 0.005;
    beads = rand(N-12,N-12);
    beads(beads>beadPercent)=0;
    beads(beads>0)=1;
    GTallbeads = beads;
    missingCols = N - size(GTallbeads,2);
    GTallbeads = padarray(GTallbeads,[0,round(missingCols/2)],'post');
    GTallbeads = padarray(GTallbeads,[0,missingCols - round(missingCols/2)],'pre');
    GTallbeads = padarray(GTallbeads,[round(missingCols/2),0],'post');
    GTallbeads = padarray(GTallbeads,[missingCols - round(missingCols/2),0],'pre');

    Vol = zeros(N,N,N);
    Vol(:,:,(N+1)/2) = GTallbeads;

    if (N+1)/2-FWHMIndex > 0
        ShiftedPSFexc = PSFexc(:,:,1:((N+1)/2+(FWHMIndex-1))); 
        ShiftedPSFexc = padarray(ShiftedPSFexc,[0 0 (N+1)/2-FWHMIndex],0,'pre');
    elseif (N+1)/2-FWHMIndex < 0 
        ShiftedPSFexc = PSFexc(:,:,(FWHMIndex-(N+1)/2):N);
        ShiftedPSFexc = padarray(ShiftedPSFexc,[0 0 (FWHMIndex-1)-(N+1)/2],0,'post');
    else
        ShiftedPSFexc = PSFexc;
    end

    PSFoverall = ShiftedPSFexc.* PSFdet;
    fftvol = fftshift(fftn(ifftshift(Vol))) .* fftshift(fftn(ifftshift(PSFoverall)));
    temp = abs(fftshift(ifftn(ifftshift(fftvol))));

    window = hamming_3d(N,N,N);
    temp = temp .* window;

    % add poiisson noise after convolution as shot noise
    Vol1 = temp + poissrnd(temp) .* 1/SNR;
    Vol2 = temp + poissrnd(temp) .* 1/SNR;

    Vol1(Vol1<0) = 0;
    Vol2(Vol2<0) = 0;