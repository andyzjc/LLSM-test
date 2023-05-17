function [ComplexPhase,Phase] = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(0.65./n*k_wave),kz_exc./(0.65./n*k_wave));
    idx = r<=1;
    Zmode = zeros(size(kx_exc));
    Zmode(idx) = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx),'norm');
    Phase = (PhaseAmplitude.* Zmode) .* 2*pi ./ wavelength_exc;
    ComplexPhase = exp( 1i .* Phase );

    