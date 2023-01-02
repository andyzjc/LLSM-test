function [Image1,Image2,ImageY,GTallbeads] = beadsSimulation(PSFexc,PSFdet)
    % Image1-> no poisson shot noise

    % Generate random beads and convolve with overallPSF(y) 
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    % determine start and end index
    range = 50.1;
    startindex = find(Y_exc >= -range,1);
    endindex = find(Y_exc <= range,1,'last');
    step = 10;

    % determine random beads image size and generate random beads
    beadsStep = round(N/step);
    beadPercent = 0.001;
    beads = rand(N,length(startindex:step:endindex)*beadsStep);
    beads(beads>beadPercent)=0;
    beads(beads>0)=1;
    GTallbeads(:,:,1) = beads;
    % add gaussian noise for autofluorescence 
    GTallbeads(:,:,2) = imnoise(beads,'gaussian');

    % loop through -50 to 50 lambda 
    Image1 = [];
    Image2 = [];
    ImageY = [];
    counter = 0;
    for i = startindex:step:endindex
        ImageIter = GTallbeads(:,1+beadsStep*counter:beadsStep+beadsStep*counter,1);
        % get PSF
        xzPSFexc = PSFexc(:,:,i);
        xzPSFdet = squeeze(PSFdet(:,i,:))'; 
        overallxzPSF = xzPSFexc .* xzPSFdet;
       
        % convolve beads with overallxzPSF
%         overallxzPSF = imresize(overallxzPSF,[size(ImageIter)]);
        ImageIter = conv2(ImageIter,overallxzPSF,'same');
        ImageIter1 = ImageIter;
        ImageIter2 = imnoise(ImageIter1,'Poisson');

        Image1 = [Image1 ImageIter1];
        Image2 = [Image2 ImageIter2];

        % add poiisson noise after convolution as shot noise
        ImageY = [ImageY Y_exc(i)];
        counter = counter + 1;
    end

%     Image = imresize(Image, [N,N]); 
    Image1 = Image1/max(max(Image1));
    Image2 = Image2/max(max(Image2));