function r = frost1972_acf(t, fc, B, Tc)
%FROST1972_ACF  Normalised autocorrelation of one band-limited source.
%   r(t) = (1 - |t|/Tc) * sinc(B t) * cos(2 pi fc t)  for |t| <= Tc, else 0.
%   sinc(x) = sin(pi x)/(pi x) (MATLAB convention).  Times in tau, frequencies in 1/tau.
r = (1 - abs(t)/Tc) .* sinc(B*t) .* cos(2*pi*fc*t);
r(abs(t) > Tc) = 0;
end
