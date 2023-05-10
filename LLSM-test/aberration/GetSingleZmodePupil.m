function Phase_factor = GetSingleZmodePupil(RadialOrder,AngularFrequency,PhaseAmplitude)
    getParameters; %modify image parameter here
    CalculatePhysics;
    
    [theta,r] = cart2pol(kx_exc./(0.65./n*k_wave),kz_exc./(0.65./n*k_wave));
    idx = r<=1;
    Phase = zeros(size(kx_exc));
    Phase(idx) = zernfun(RadialOrder,AngularFrequency,r(idx),theta(idx),'norm');
    Phase_factor = exp(PhaseAmplitude.* 1i .* Phase );