function ok = frost1972_verify()
%FROST1972_VERIFY  Self-tests for the Frost (1972) reproduction.  Prints PASS/FAIL per test.
%   Run this before trusting any figure produced by frost1972_run.m.
S = frost1972_setup();
ok = true; nt = 0;
    function check(cond, msg, varargin)
        nt = nt + 1;
        if cond, tag = 'PASS'; else, tag = 'FAIL'; ok = false; end
        fprintf('[%s] T%02d %s\n', tag, nt, sprintf(msg, varargin{:}));
    end

fprintf('=== Frost (1972) reproduction: self-verification ===\n');

% ---- T1-T4: constraint algebra of Sec. II/III ----------------------------
check(rank(S.C) == S.J, 'C has full rank J = %d (eq. 8)', S.J);
check(norm(S.P*S.P - S.P, 'fro') < 1e-12 && norm(S.P - S.P','fro') < 1e-12, ...
      'P is a symmetric idempotent projector (P^2 = P)');
check(norm(S.C'*S.P, 'fro') < 1e-12, 'C''P = 0 : P annihilates the constraint directions');
check(norm(S.C'*S.F - S.Fvec) < 1e-12 && norm(S.P*S.F) < 1e-12, ...
      'F = C(C''C)^{-1}Fvec satisfies C''F = Fvec and lies normal to the plane (PF = 0)');

% ---- T5: Rxx is a valid covariance and the PSD model is non-negative ------
check(min(eig(S.Rxx)) > 0, 'Rxx is positive definite (min eig = %.4f)', min(eig(S.Rxx)));
[~, info] = frost1972_waveforms(S, 4, 1, 50);
minS = min(cellfun(@(p) min(p.S), info.psd));
check(minS >= -1e-9, 'model PSD of every source is non-negative (min = %.2e)', minS);

% ---- T6: Wopt from eq. (16) equals the KKT solution of (11) ---------------
KKT = [S.Rxx, S.C; S.C', zeros(S.J)] \ [zeros(S.N,1); S.Fvec];
check(norm(KKT(1:S.N) - S.Wopt) < 1e-9, 'Wopt (eq. 16) = KKT solution of min W''RW s.t. C''W = Fvec');
check(norm(S.C'*S.Wopt - S.Fvec) < 1e-12, 'Wopt satisfies the constraint');
% random feasible perturbations never lower the output power
Z = S.P*randn(S.N, 200);
Pz = sum((S.Wopt + Z) .* (S.Rxx*(S.Wopt + Z)), 1);
check(all(Pz >= S.Popt - 1e-12), 'no feasible perturbation lowers W''RW below Popt = %.4f', S.Popt);

% ---- T7: eigen-structure of P Rxx P (Sec. IV) ----------------------------
check(sum(abs(S.sigma_all) < 1e-9) == S.J, 'P Rxx P has exactly J = %d zero eigenvalues', S.J);
check(min(S.sigma) >= min(S.lambda) - 1e-9 && max(S.sigma) <= max(S.lambda) + 1e-9, ...
      'lambda_min <= sigma_i <= lambda_max  (%.4f <= [%.4f, %.4f] <= %.4f)', ...
      min(S.lambda), min(S.sigma), max(S.sigma), max(S.lambda));

% ---- T8: agreement with the numbers printed in the paper -----------------
check(abs(min(S.lambda) - 0.111) < 0.002 && abs(max(S.lambda) - 8.355) < 0.02, ...
      'eigenvalues of Rxx = [%.3f, %.3f]  (paper: 0.111 ... 8.355)', min(S.lambda), max(S.lambda));
check(abs(S.M_lo - 0.152) < 0.006 && abs(S.M_hi - 0.170) < 0.006, ...
      'misadjustment bounds (eq. 27) at mu = 0.01: %.1f %% ... %.1f %%  (paper: 15.2 ... 17.0 %%)', ...
      100*S.M_lo, 100*S.M_hi);
fprintf('       eq.(26) bound = %.4f, eq.(28) bound = %.4f, 1/sigma_max = %.4f  (paper quotes 0.074 for eq. 26)\n', ...
        S.mu_max26, S.mu_max28, S.mu_stab);

% ---- T9: time-domain waveform generator reproduces the analytic Rxx ------
S2 = S; S2.Delta = 32;                           % > 25 + tap/sensor span -> snapshots exactly independent
[X, ~] = frost1972_waveforms(S2, 20000, 7, 50);
Rhat = X*X'/size(X,2);
relerr = norm(Rhat - S.Rxx, 'fro')/norm(S.Rxx, 'fro');
check(relerr < 0.03, 'sample covariance of generated waveforms matches analytic Rxx (rel. Frobenius err = %.3f)', relerr);
check(abs(trace(Rhat) - trace(S.Rxx))/trace(S.Rxx) < 0.02, 'total tap power %.2f vs analytic %.2f', trace(Rhat), trace(S.Rxx));

% ---- T10: deterministic algorithm (eq. 21) converges to Wopt -------------
W = S.F; for k = 1:20000, W = S.P*(W - S.mu*S.Rxx*W) + S.F; end   % slowest time constant ~ 900 adaptations
check(norm(W - S.Wopt) < 1e-6, 'deterministic constrained gradient descent (eq. 21) converges to Wopt');

% ---- T11: constrained LMS keeps the constraint exactly (double precision) -
r = frost1972_adapt(S, X(:,1:2000), S.mu, 'clms', []);
check(max(r.cerr) < 1e-24, 'CLMS (eq. 22): ||C''W(k) - Fvec||^2 <= %.1e at every iteration', max(r.cerr));

% ---- T12: mean weight vector converges to Wopt with the predicted rate ----
Nrun = 400; Nit = 300;
Wmean = zeros(S.N, Nit+1);
L = chol(S.Rxx, 'lower');
for n = 1:Nrun
    rng(1000+n);
    Xn = L*randn(S.N, Nit);                     % independent Gaussian snapshots
    rr = frost1972_adapt(S, Xn, S.mu, 'clms', []);
    Wmean = Wmean + rr.W/Nrun;
end
dev  = sqrt(sum((Wmean - S.Wopt).^2, 1));
dev0 = norm(S.F - S.Wopt);
PRP = S.P*S.Rxx*S.P; Vpred = zeros(S.N, Nit+1); Vpred(:,1) = S.F - S.Wopt;
for k = 1:Nit, Vpred(:,k+1) = Vpred(:,k) - S.mu*PRP*Vpred(:,k); end      % eq. (24)
errpred = max(sqrt(sum((Wmean - S.Wopt - Vpred).^2, 1)))/dev0;
check(errpred < 0.05, 'ensemble mean E[W(k)] follows eq. (24) V(k+1) = (I - mu P Rxx P) V(k): max dev %.3f of ||V(0)||', errpred);
check(dev(end) < dev0*exp(-Nit*S.mu*min(S.sigma))*1.3, '||E[W]-Wopt|| drops from %.3f to %.3f after %d adaptations (slow-mode prediction %.3f)', dev0, dev(end), Nit, dev0*exp(-Nit*S.mu*min(S.sigma)));
check(all(dev(2:end) <= dev0*(1-S.mu*min(S.sigma)).^(1:Nit)*1.15), ...
      'mean deviation stays under the slow-mode bound (1 - mu sigma_min)^k ||V(0)||');
Vk = Wmean - S.Wopt;                             % mean difference process
check(norm(S.C'*Vk, 'fro') < 1e-9, 'difference vector V(k) stays in the constraint subspace (C''V = 0)');

% ---- T13: steady-state misadjustment lies inside the eq. (27) bounds -----
Nrun = 400; Nit = 1500;  Pss = 0;
for n = 1:Nrun
    rng(5000+n);
    Xn = L*randn(S.N, Nit);
    rr = frost1972_adapt(S, Xn, S.mu, 'clms', []);
    Pss = Pss + mean(rr.Pout(1001:end))/Nrun;    % average steady-state output power
end
Mhat = (Pss - S.Popt)/S.Popt;
check(Mhat > S.M_lo - 0.02 && Mhat < S.M_hi + 0.02, ...
      'measured misadjustment %.1f %% is within [%.1f, %.1f] %% (eq. 27, tolerance 2 %%)', 100*Mhat, 100*S.M_lo, 100*S.M_hi);

% ---- T14: with 6-digit arithmetic GP drifts, CLMS does not --------------
Xq = X(:, 1:2000);
c6 = frost1972_adapt(S, Xq, S.mu, 'clms', 6);
g6 = frost1972_adapt(S, Xq, S.mu, 'gp',   6);
check(mean(c6.cerr(1001:end)) < 5*mean(c6.cerr(2:500)) , ...
      'CLMS 6-digit constraint error does not grow (first 500: %.2e, last 1000: %.2e)', mean(c6.cerr(2:500)), mean(c6.cerr(1001:end)));
check(mean(g6.cerr(1501:end)) > 5*mean(g6.cerr(2:200)), ...
      'GP 6-digit constraint error grows (first 200: %.2e, last 500: %.2e)', mean(g6.cerr(2:200)), mean(g6.cerr(1501:end)));
% GP random walk: variance roughly linear in k  -> fit log-log slope near 1
kk = 50:2000; p = polyfit(log(kk), log(g6.cerr(kk+1)+eps), 1);
fprintf('       GP constraint error growth exponent (log-log slope) = %.2f  (> 0: errors accumulate, cf. Fig. 12)\n', p(1));

if ok, fprintf('=== %d tests: ALL TESTS PASSED ===\n', nt); else, fprintf('=== %d tests: SOME TESTS FAILED ===\n', nt); end
end
