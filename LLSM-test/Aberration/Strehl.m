function [SW_SRatio,Lattice_SRatio,RadioOrderArray,AngularFrequencyArray] = Strehl(SWPupil,LatticePupil,MinRadialOrder,MaxRadialOrder,PhaseAmplitude,PSFdet)
    getParameters; %modify image parameter here
    CalculatePhysics;

    SWPSF = zeros(N,N);
    for j = 1:size(SWPupil,3)
        PupilIter = SWPupil(:,:,j);
        SWPSF = SWPSF + abs(fftshift(ifft2(fftshift(PupilIter)))).^2;
    end
    [SWvalue,~] = max(SWPSF,[],'all'); % value, index
    SWPSF = SWPSF./SWvalue;


    LatticePSF = abs(fftshift(ifft2(fftshift(LatticePupil)))).^2;
    LatticePSF = meshgrid(mean(LatticePSF,2))';
    [Latticevalue,~] = max(LatticePSF,[],'all'); % value, index
    LatticePSF = LatticePSF./Latticevalue;
    
    counter = 1;
    for i = MinRadialOrder:MaxRadialOrder
        AngularFrequency = -i:2:i;
        for k = 1:length(AngularFrequency)
            RadioOrderArray(1,counter) = i;
            AngularFrequencyArray(counter) = AngularFrequency(k);
            
            [ComplexPhase,phase] = GetSingleZmodePupil(i,AngularFrequency(k),PhaseAmplitude);
    
            AberratedSWPSF = zeros(size(phase));
            %SW 
            for j = 1:size(SWPupil,3)
                PupilIter = SWPupil(:,:,j);
                AberratedSWPupil = PupilIter .* ComplexPhase;
                temp = abs(fftshift(ifft2(fftshift(AberratedSWPupil)))).^2;
                AberratedSWPSF = AberratedSWPSF + temp;
            end
            AberratedSWPSF = AberratedSWPSF./SWvalue; 
            SW_SRatio(counter,1) = max(AberratedSWPSF.*PSFdet(:,:,(N+1)/2),[],'all'); % only at focal point/ on optical axis

            AberratedLatticePupil = zeros(size(phase));
            AberratedLatticePupil = LatticePupil .* exp( 1i.*phase);
            AberratedLatticePSF = abs(fftshift(ifft2(fftshift(AberratedLatticePupil)))).^2;
            AberratedLatticePSF = meshgrid(mean(AberratedLatticePSF,2))';
            AberratedLatticePSF = AberratedLatticePSF/Latticevalue;
            Lattice_SRatio(counter,1) = max(AberratedLatticePSF.*PSFdet(:,:,(N+1)/2),[],'all');
            counter = counter + 1;
        end
    end

    LabelArray = cell(length(RadioOrderArray),1);
    for i = 1:length(RadioOrderArray)
        LabelArray{i,1} = "Z^{" + num2str(AngularFrequencyArray(i)) + "}_{" + num2str(RadioOrderArray(i)) + "}";
    end

    fig1 = figure;
    h1 = subplot(1,1,1,'Parent',fig1);
    plot(1:length(SW_SRatio),SW_SRatio,'Parent',h1,'LineStyle','-','Marker','o','MarkerSize',15);
    hold on
    plot(1:length(Lattice_SRatio),Lattice_SRatio,'Parent',h1,'LineStyle','-.','Marker','o','MarkerSize',15);
    lgd = legend("iSW","LLS");
    lgd.Location = 'northoutside';
    h1.XAxis.TickValues = 1:length(RadioOrderArray);
    h1.XAxis.TickLabels = LabelArray;
    h1.XAxis.FontSize = 6;
    grid on
    xlabel("Aberration Mode")
    ylabel("Strehl Ratio")
    title("Focal Plane")
    pbaspect([5 1 1])
    hold off
