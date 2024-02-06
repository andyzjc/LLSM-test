function zFWHM = CalculateResolution(excyzPSF)
    getParameters; %modify image parameter here
    CalculatePhysics;

    excyzPSF = excyzPSF./max(excyzPSF,[],'all');
    %  calculate yFWHM 
    yindex = 1-(squeeze(excyzPSF((N+1)/2,:)) <= 0.5*max(squeeze(excyzPSF((N+1)/2,:))));
    yFWHM1 = find(yindex,1,'first') ;
    yFWHM2 = find(yindex,1,'last');
    
    % input a xz excitation PSF, get only one z excitation PSF
    maxZ = max(Z_exc);
    for i = 1:size(excyzPSF,2)
        zPSF = excyzPSF(:,i);
        % re-sample to a higher resolution and the Z_exc
        sample_rate = 4;
        zPSF = resample(zPSF,sample_rate,1);
        zPSF = zPSF(2:end-4);
       
        M = length(zPSF);
        Z_exc_new = linspace(-maxZ,maxZ,M);
        % [~,maxindex] = max(overallzPSF);        
        % index1 = (overallzPSF(maxindex:end) <= 0.5*max(overallzPSF) );
        % index2 = (overallzPSF(1:maxindex) <= 0.5*max(overallzPSF));

        maxindex = (M-1)/2+1; % calculate on z = focus 
        index1 = zPSF(maxindex:end) <= 0.5*zPSF(maxindex) ;
        index2 = zPSF(1:maxindex) <= 0.5*zPSF(maxindex) ;
        if ~isempty(index1) && ~isempty(index2)
            zFWHM1(i,1) = Z_exc_new(maxindex + find(index1,1,'first')) ;
            zFWHM2(i,1) = Z_exc_new(find(index2,1,'last')-1);
            if abs(zFWHM1(i,1)) == abs(zFWHM2(i,1))
                zFWHM(i,1) = abs(zFWHM1(i,1)) + abs(zFWHM2(i,1));
            elseif abs(Z_exc_new(maxindex) - zFWHM2(i,1)) > abs(Z_exc_new(maxindex) - zFWHM2(i,1))
                zFWHM(i,1) = abs(Z_exc_new(maxindex) - zFWHM1(i,1))*2;
            else
                zFWHM(i,1) = abs(Z_exc_new(maxindex) - zFWHM2(i,1))*2;
            end
        else
            zFWHM(i,1) = "N/A";
        end

        % plot(Z_exc_new,zPSF)
        % hold on 
        % drawnow
        % yline(0.5*zPSF(maxindex))
        % xlim([-10,10])
        % grid on
    end

    % scatter(Y_exc,zFWHM)
    % hold on
    % xline(Y_exc(yFWHM2));
    % grid on
    % ylabel("zFWHM (\lambda/n)")
    % xlabel("Y (\lambda/n)")
    % xlim([0,50])
end

