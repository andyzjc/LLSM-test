function [Image1,Image2] = beadsSimulation(PSFexc,PSFdet,SNR)
    % Generate random beads and convolve with overallPSF(y)
%      GTallbeads = importdata("GTbeads,rho=0.001.mat");

    getParameters; %modify image parameter here
    CalculatePhysics;
    xzPSFdet = squeeze(PSFdet(:,:,(N+1)/2));

    % determine start and end index
    range = 50.1; % -50 to 50 lambda/n
    startindex = find(Y_exc >= -range,1);
    endindex = find(Y_exc <= range,1,'last');
    step = 10;
    Yindex = startindex:step:endindex;

    % determine random beads image size and generate random beads
    beadsStep = round(N/length(Yindex))-1;
    beadPercent = 0.01;
    beads = rand(length(startindex:step:endindex)*beadsStep,length(startindex:step:endindex)*beadsStep);
    beads(beads>beadPercent)=0;
    beads(beads>0)=1;
    GTallbeads = beads;
    missingCols = N - size(GTallbeads,2);
    GTallbeads = padarray(GTallbeads,[0,round(missingCols/2)],'post');
    GTallbeads = padarray(GTallbeads,[0,missingCols - round(missingCols/2)],'pre');
    GTallbeads = padarray(GTallbeads,[round(missingCols/2),0],'post');
    GTallbeads = padarray(GTallbeads,[missingCols - round(missingCols/2),0],'pre');

    % loop through -50 to 50 lambda 
    Image1 = zeros(N,N);
    counter = 0;   
    for i = startindex:step:endindex
        ImageIter = GTallbeads(:,1+beadsStep*counter:beadsStep+beadsStep*counter);
        ImageIter = padarray(ImageIter,[0,beadsStep*counter],'pre');
        ImageIter = padarray(ImageIter,[0,size(GTallbeads,2)-(beadsStep+beadsStep*counter)],'post');
        missingCols = N - size(ImageIter,2);
        ImageIter = padarray(ImageIter,[0,round(missingCols/2)],'post');
        ImageIter = padarray(ImageIter,[0,missingCols - round(missingCols/2)],'pre');

        % get PSF
        xzPSFexc = PSFexc(:,:,i);
        xzPSFoverall = xzPSFexc .* xzPSFdet;
       
        % convolve beads with overallxzPSF
        ImageIter1 = conv2(ImageIter,xzPSFoverall,'same');
        Image1 = Image1+ImageIter1;

        counter = counter + 1;
    end
    % add poiisson noise after convolution as shot noise
    temp = Image1;
    Image1 = temp + poissrnd(temp) * 1/SNR;
    Image2 = temp + poissrnd(temp) * 1/SNR;

    Image1 = Image1/max(max(Image1));
    Image2 = Image2/max(max(Image2));

    % apply a hamming window 
    hmwindow = hamming(N)*hamming(N)';
    Image1 = Image1.*hmwindow;
    Image2 = Image2.*hmwindow;

    Image1(Image1<0) = 0;
    Image2(Image2<0) = 0;