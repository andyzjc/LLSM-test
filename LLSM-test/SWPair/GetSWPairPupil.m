function [SWPupil,SWMask] = GetSWPairPupil(ProfileType,...
                                  NA1Ideal,NA2Ideal,...
                                  deltaNA1,deltaNA2,...
                                  WeightRatio)
% generates a pair SW pupil functions 
% SWPupil and SWMask = NxNx2
% SWPupil(:,:,1) = NA1Pupil, SWMask(:,:,1) = SWMask

getParameters;
CalculatePhysics;


