function [Vol1,Vol2,GTallbeads] = beadsSimulation3d(PSFexc, PSFdet,SNR)
    getParameters; %modify image parameter here
    CalculatePhysics;

    % determine start and end index
    range = 50.1; % -50 to 50 lambda/n
    startindex = find(Y_exc >= -range,1);
    endindex = find(Y_exc <= range,1,'last');
    step = 20;
    Yindex = startindex:step:endindex;

    % determine random beads image size and generate random beads
    beadsStep = round(N/length(Yindex))-1;
    beadPercent = 0.0001;
    beads = rand(length(startindex:step:endindex)*beadsStep,...
                 length(startindex:step:endindex)*beadsStep,...
                 length(startindex:step:endindex)*beadsStep);
    beads(beads>beadPercent)=0;
    beads(beads>0)=1;
    GTallbeads = beads;
    missingCols = N - size(GTallbeads,2);
    GTallbeads = padarray(GTallbeads,[0,round(missingCols/2),0],'post');
    GTallbeads = padarray(GTallbeads,[0,missingCols - round(missingCols/2),0],'pre');
    GTallbeads = padarray(GTallbeads,[round(missingCols/2),0,0],'post');
    GTallbeads = padarray(GTallbeads,[missingCols - round(missingCols/2),0,0],'pre');
    GTallbeads = padarray(GTallbeads,[0,0,round(missingCols/2)],'post');
    GTallbeads = padarray(GTallbeads,[0,0,missingCols - round(missingCols/2)],'pre');

    Vol = zeros(N,N,N);
    counter = 0;   
    for i = startindex:step:endindex
        ImageIter = GTallbeads(:,:,1+beadsStep*counter:beadsStep+beadsStep*counter);
        ImageIter = padarray(ImageIter,[0,0,beadsStep*counter],'pre');
        ImageIter = padarray(ImageIter,[0,0,size(GTallbeads,2)-(beadsStep+beadsStep*counter)],'post');

        % shift detection PSF to correct slice
        if (N+1)/2-i > 0
            ShiftedPSFexc = PSFexc(:,:,1:((N+1)/2+(i-1))); 
            ShiftedPSFexc = padarray(ShiftedPSFexc,[0 0 (N+1)/2-i],0,'pre');
        elseif (N+1)/2-i < 0 
            ShiftedPSFexc = PSFexc(:,:,(i-(N+1)/2):N);
            ShiftedPSFexc = padarray(ShiftedPSFexc,[0 0 (i-1)-(N+1)/2],0,'post');
        else
            ShiftedPSFexc = PSFexc;
        end
        overallPSF = ShiftedPSFexc.* PSFdet;
       
        % convolve beads with overallxzPSF
        temp = fftshift(fftn(ifftshift(ImageIter))) .* fftshift(fftn(ifftshift(overallPSF)));
        ImageIter3d = abs(fftshift(ifftn(ifftshift(temp))));
        Vol = Vol+ImageIter3d;

        counter = counter + 1;
    end

        % add poiisson noise after convolution as shot noise
    temp = Vol;
    Vol1 = temp + poissrnd(temp) .* 1/SNR;
    Vol2 = temp + poissrnd(temp) .* 1/SNR;

    Vol1(Vol1<0) = 0;
    Vol2(Vol2<0) = 0;
