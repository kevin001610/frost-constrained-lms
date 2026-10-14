function S = frost1972_setup()
%FROST1972_SETUP  Analytic model of the simulation in Frost (1972), Sec. VII / Table I.
%
%   S = FROST1972_SETUP() returns a struct with the array geometry, the source
%   table, the analytic tap-voltage correlation matrix Rxx, the constraint
%   matrix C, the look-direction filter Fvec, the projection P, the shortest
%   constraint-plane vector F, the optimum weight vector Wopt (eq. 16) and the
%   performance quantities of Sec. IV (eigenvalues, mu bounds, misadjustment).
%
%   Model (all times in units of the tap delay tau, frequencies in units of 1/tau):
%     * K = 4 sensors on a line, inter-sensor propagation time tau (so a plane
%       wave from angle theta off broadside reaches sensor k with delay k*tau*sin(theta)).
%     * J = 4 taps per sensor, tap delay tau  ->  KJ = 16 weights.
%     * Each point source is a zero-mean Gaussian process with power P and a
%       band-pass spectrum of centre fc and (full) bandwidth B.  The paper states
%       that all temporal correlations are identically zero for lags > 25 tau; we
%       realise this with the triangular (Bartlett) taper, i.e.
%           r(t) = P * (1 - |t|/25) * sinc(B t) * cos(2 pi fc t),   |t| <= 25,
%       which is a valid (non-negative-PSD) autocorrelation because it is the
%       product of three valid autocorrelations.
%     * Independent white noise of power 0.1 is added on every tap.
%
%   The model is chosen only from the information printed in the paper; its
%   quality is checked in frost1972_verify.m (eigenvalues 0.111 ... 8.355,
%   misadjustment 15.2 ... 17.0 %, optimum output power of Fig. 11).

S.K   = 4;            % sensors
S.J   = 4;            % taps per sensor
S.N   = S.K*S.J;      % number of weights (16)
S.Tc  = 25;           % correlation support in tau (Sec. VII: zero beyond 25 tau)
S.Delta = 60;         % adaptation interval in tau (paper: > 58 tau -> independent snapshots)

% ---- Table I ------------------------------------------------------------
S.src = struct( ...
    'name', {'Look-direction signal','Noise A','Noise B'}, ...
    'P',    {0.1, 1.0, 1.0}, ...          % power
    'theta',{0,   45,  60 }, ...          % degrees from broadside
    'fc',   {0.3, 0.2, 0.4}, ...          % centre frequency (1.0 = 1/tau)
    'B',    {0.1, 0.05,0.07});            % bandwidth
S.Pwhite = 0.1;                            % white noise power per tap

% ---- look-direction filter (Sec. VII) -----------------------------------
S.Fvec = [1; -2; 1.5; 2];                  % script-F, J x 1

% ---- constraint matrix C (eq. 7-8): c_j selects the j-th column of K taps --
% weight index n = (j-1)*K + k   (k = sensor 1..K, j = tap 1..J)
S.C = kron(eye(S.J), ones(S.K,1));         % KJ x J, full rank J
S.P = eye(S.N) - S.C*((S.C'*S.C)\S.C');    % eq. 20b
S.F = S.C*((S.C'*S.C)\S.Fvec);             % eq. 20a (= W(0))

% ---- analytic correlation matrix Rxx (eq. 2a) ---------------------------
n  = (1:S.N)';                             % weight index n = j*K + k + 1  (tap-major, matches C)
kk = mod(n-1, S.K);                        % sensor index 0..K-1
jj = floor((n-1)/S.K);                     % tap index    0..J-1
R = zeros(S.N);
for m = 1:numel(S.src)
    d = sind(S.src(m).theta);              % inter-sensor delay (in tau)
    lag = (jj - jj') + (kk - kk')*d;       % lag between tap a and tap b
    R = R + S.src(m).P * frost1972_acf(lag, S.src(m).fc, S.src(m).B, S.Tc);
end
R = R + S.Pwhite*eye(S.N);
S.Rxx = (R + R')/2;

% ---- optimum weight vector (eq. 16) -------------------------------------
RiC    = S.Rxx \ S.C;
S.Wopt = RiC * ((S.C'*RiC) \ S.Fvec);
S.Popt = S.Fvec' * ((S.C'*RiC) \ S.Fvec);  % optimum output power (Sec. IV)

% ---- Sec. IV performance quantities -------------------------------------
S.lambda  = sort(eig(S.Rxx));              % eigenvalues of Rxx
PRP       = S.P*S.Rxx*S.P;  PRP = (PRP+PRP')/2;
S.sigma_all = sort(eig(PRP));              % J of them are (numerically) zero
S.sigma   = S.sigma_all(S.J+1:end);        % the KJ-J non-zero eigenvalues
S.trPRP   = trace(PRP);
S.mu      = 0.01;                          % Sec. VII
S.mu_max26 = 1/(max(S.sigma) + 0.5*S.trPRP);     % eq. 26
S.mu_max28 = 2/(3*trace(S.Rxx));                 % eq. 28
S.mu_stab  = 1/max(S.sigma);                     % mean-convergence condition
mu = S.mu;
S.M_lo = (mu/2)*S.trPRP / (1 - (mu/2)*(S.trPRP + 2*min(S.sigma)));   % eq. 27 lower
S.M_hi = (mu/2)*S.trPRP / (1 - (mu/2)*(S.trPRP + 2*max(S.sigma)));   % eq. 27 upper
S.tau_slow = 1/(mu*min(S.sigma));          % eq. 25, in adaptations
S.tau_fast = 1/(mu*max(S.sigma));
end
