function Normalization_factors = zernfun_normalization(n,m)

% Zernike Normalization factors for generating orthonormal Zernikes.  This
% will make polynomials of unit length RMS  (e.g. um RMS).
%
    Normalization_factors =  sqrt( (2. * n + 2) / (1. + (m == 0)) ); 
end

