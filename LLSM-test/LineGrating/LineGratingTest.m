testPSF = OverallAberratedPSFIncoherent;
zPSF = squeeze(testPSF(:,(N+1)/2,(N+1)/2));

shiftdistance = 0.2; %um
shiftdistance_pixel = round(shiftdistance/deltax)

close all
figure
plot(Z_exc,zPSF)

figure
plot(Z_exc,zPSF+circshift(zPSF,shiftdistance_pixel))
hold on
xline(0)

figure
plot(Z_exc,zPSF+circshift(zPSF,shiftdistance_pixel) + circshift(zPSF,shiftdistance_pixel*2) + circshift(zPSF,shiftdistance_pixel*3) + circshift(zPSF,shiftdistance_pixel+4))
hold on
xline(0)

figure
plot(Z_exc,zPSF + circshift(zPSF,4) )
hold on
xline(0)