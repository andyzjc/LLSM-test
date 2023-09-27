function [IncoherentSumI,absIncoherentSumI,IncoherentIFWHM,half] = IFWHM(newpupil)
    getParameters; %modify image parameter here
    CalculatePhysics;
    xzPSFbound = abs(fftshift(ifft2(ifftshift(newpupil)))).^2;

% first let's integrate the unbound one, to get the energy FWHM 
    IncoherentSumI = cumsum(sum(xzPSFbound,1));
    IncoherentSumI = IncoherentSumI/max(IncoherentSumI,[],'all');
    [~,yzeropoint1,~,~] = intersections(1:length(IncoherentSumI),IncoherentSumI, ...
                                        [(N+1)/2,(N+1)/2],[0,1]);
    [~,yzeropoint2,~,~] = intersections(1:length(IncoherentSumI),IncoherentSumI, ...
                                        [(N+1)/2-1,(N+1)/2-1],[0,1]);
    yzeropoint = (yzeropoint1(1,1)+yzeropoint2(1,1))/2;
    absIncoherentSumI = abs(IncoherentSumI-yzeropoint);
    temp = zeros(1,length(absIncoherentSumI)+1);
    temp(1:(N+1)/2-1) = absIncoherentSumI(1:(N+1)/2-1);
    temp(1,(N+1)/2) = 0;
    temp(1,(N+1)/2+1:end) = absIncoherentSumI((N+1)/2:end);
    absIncoherentSumI = 2* temp(1:end-1);
    % find IFWHM
    [~,minindex] = min(absIncoherentSumI);
    [~,maxindex] = max(absIncoherentSumI);
    half = (absIncoherentSumI(maxindex)-absIncoherentSumI(minindex))/2;
    index1 = (absIncoherentSumI(minindex:end) <= half);
    index2 = (absIncoherentSumI(1:minindex) <= half);
    if ~isempty(index1) && ~isempty(index2)
        IncoherentIFWHM1 = Z_exc(minindex + find(index1,1,'last')-1);
        IncoherentIFWHM2 = Z_exc(find(index2,1,'first'));
        if abs(IncoherentIFWHM1) == abs(IncoherentIFWHM2)
            IncoherentIFWHM = abs(IncoherentIFWHM1) + abs(IncoherentIFWHM2);
        elseif abs(Z_exc(minindex) - IncoherentIFWHM1) > abs(Z_exc(minindex) - IncoherentIFWHM2)
            IncoherentIFWHM = abs(Z_exc(minindex) - IncoherentIFWHM1)*2;
        else
            IncoherentIFWHM = abs(Z_exc(minindex) - IncoherentIFWHM2)*2;
        end
    else
        IncoherentIFWHM = "N/A";
    end
    IncoherentIFWHM = IncoherentIFWHM/0.25; %to pixel



end