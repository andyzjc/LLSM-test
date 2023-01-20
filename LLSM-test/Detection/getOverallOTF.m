function xzOTFOverall = getOverallOTF(xzOTFexc)
    xzOTFdet = ScaleDetectionOTF;
    xzOTFOverall = conv2(xzOTFexc,xzOTFdet,'same');
    
%     plot(KZ_exc,log10(xzOTFOverall((N+1)/2,:)))
%     hold on
%     plot(KZ_exc,log10(xzOTFOverall(:,(N+1)/2)))
%     xlim([-1,1])
%     ylim([-3,0])
%     grid on