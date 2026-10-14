%FROST1972_RUN  Reproduce the simulation of Frost (1972), Sec. VII, Figs. 9-12.
%   Writes PNG figures and a results summary into ./results (next to this file).
%   Run frost1972_verify first; this script also prints the key numbers.
clear; clc;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end

S = frost1972_setup();

% ---- colours (deck palette: deep blue / teal = constrained LMS / orange = gradient projection)
cBlue = [27 58 107]/255; cTeal = [14 158 142]/255; cOrange = [232 119 34]/255; cGrey = [110 110 110]/255;
set(groot, 'defaultAxesFontName', 'Calibri', 'defaultAxesFontSize', 13, 'defaultTextFontName', 'Calibri');
savef = @(h, name) exportgraphics(h, fullfile(outdir, name), 'Resolution', 200, 'BackgroundColor', 'white');

%% ---- 1. look-direction frequency response (Fig. 9) -------------------------
f = linspace(0, 0.5, 501);
H = exp(-1i*2*pi*f'*(0:S.J-1)) * S.Fvec;                % H(f) = sum_j f_j e^{-i 2 pi f (j-1)}
h = figure('Color', 'w', 'Position', [100 100 800 450]);
plot(f, abs(H).^2, 'Color', cBlue, 'LineWidth', 2.2); grid on;
xlabel('Frequency  (1.0 = 1/\tau)'); ylabel('|H(f)|^2'); xlim([0 0.5]);
title('Look-direction response of  F^T = [1, -2, 1.5, 2]   (cf. Frost Fig. 9)');
savef(h, 'fig09_look_direction_response.png');

%% ---- 2. source spectra (Fig. 10) -------------------------------------------
[~, info] = frost1972_waveforms(S, 2, 1, 50);
h = figure('Color', 'w', 'Position', [100 100 800 450]); hold on; grid on;
cols = {cBlue, cTeal, cOrange};
for m = 1:3
    ff = info.psd{m}.f; SS = info.psd{m}.S; sel = ff >= 0 & ff <= 0.5;
    [ffs, idx] = sort(ff(sel)); SSs = SS(sel); SSs = SSs(idx);
    plot(ffs, SSs, 'Color', cols{m}, 'LineWidth', 2.2, 'DisplayName', ...
        sprintf('%s  (P = %.1f, %d^\\circ)', S.src(m).name, S.src(m).P, S.src(m).theta));
end
xlabel('Frequency  (1.0 = 1/\tau)'); ylabel('Power spectral density'); xlim([0 0.5]);
legend('Location', 'northeast'); title('Power spectra of the simulated sources (cf. Frost Fig. 10 / Table I)');
savef(h, 'fig10_source_spectra.png');

%% ---- 3. single-run adaptation on real waveforms, 6-digit arithmetic (Figs. 11, 12)
Nad = 200;
[X, ~] = frost1972_waveforms(S, Nad, 2026, 100);        % snapshots 60 tau apart
clms = frost1972_adapt(S, X, S.mu, 'clms', 6);
gp   = frost1972_adapt(S, X, S.mu, 'gp',   6);           % exactly the same data
k = 0:Nad;
bounds = S.Popt*[1, 1+S.M_lo, 1+S.M_hi];

for which = 1:2
    if which == 1, r = clms; c = cTeal; nm = 'Constrained LMS (eq. 22)'; fn = 'fig11_clms_output_power.png';
    else,          r = gp;   c = cOrange; nm = 'Gradient projection (eq. 31)'; fn = 'fig12_gp_output_power.png'; end
    h = figure('Color', 'w', 'Position', [100 100 900 620]);
    t = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile([2 1]); hold on; grid on;
    plot(k, r.Pout, '.', 'Color', c, 'MarkerSize', 9, 'DisplayName', nm);
    yline(bounds(1), '-', 'Color', cBlue, 'LineWidth', 1.6, 'DisplayName', sprintf('optimum output power = %.2f', bounds(1)));
    yline(bounds(2), '--', 'Color', cGrey, 'LineWidth', 1.2, 'DisplayName', sprintf('optimum + misadjustment (%.1f%% / %.1f%%)', 100*S.M_lo, 100*S.M_hi));
    yline(bounds(3), '--', 'Color', cGrey, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    ylim([0 7]); xlim([0 Nad]); ylabel('Output power  W^T R_{xx} W'); legend('Location', 'northeast');
    title(sprintf('%s,  \\mu = %.2g,  6-digit arithmetic  (MATLAB reproduction of Frost Fig. %d)', nm, S.mu, 10+which));
    ax2 = nexttile; hold on; grid on;
    plot(k, r.cerr, '.', 'Color', c, 'MarkerSize', 8);
    xlim([0 Nad]); ylabel('||C^T W - F||^2'); xlabel('Number of adaptations');
    ax2.YAxis.Exponent = -11;
    savef(h, fn);
end

% constraint-error comparison on one axis
h = figure('Color', 'w', 'Position', [100 100 900 420]); hold on; grid on;
plot(k, clms.cerr, '.-', 'Color', cTeal, 'MarkerSize', 8, 'LineWidth', 0.8, 'DisplayName', 'Constrained LMS');
plot(k, gp.cerr,   '.-', 'Color', cOrange, 'MarkerSize', 8, 'LineWidth', 0.8, 'DisplayName', 'Gradient projection');
xlabel('Number of adaptations'); ylabel('||C^T W(k) - F||^2'); legend('Location', 'northwest');
ax = gca; ax.YAxis.Exponent = -11; xlim([0 Nad]);
title('Constraint error with 6-digit arithmetic: same data, same \mu');
savef(h, 'fig_constraint_error_compare.png');

% long run: GP random walk vs CLMS over 2000 adaptations (independent Gaussian snapshots)
L = chol(S.Rxx, 'lower'); rng(99); Xl = L*randn(S.N, 2000);
cl = frost1972_adapt(S, Xl, S.mu, 'clms', 6); gl = frost1972_adapt(S, Xl, S.mu, 'gp', 6);
h = figure('Color', 'w', 'Position', [100 100 900 420]); hold on; grid on;
semilogy(0:2000, max(cl.cerr, 1e-13), '-', 'Color', cTeal, 'LineWidth', 1.3, 'DisplayName', 'Constrained LMS');
semilogy(0:2000, max(gl.cerr, 1e-13), '-', 'Color', cOrange, 'LineWidth', 1.3, 'DisplayName', 'Gradient projection');
set(gca, 'YScale', 'log'); xlabel('Number of adaptations'); ylabel('||C^T W(k) - F||^2  (log, floor 10^{-13})');
legend('Location', 'southeast'); title('Error accumulation over 2000 adaptations (6-digit arithmetic)');
savef(h, 'fig_constraint_error_long.png');

%% ---- 4. ensemble: mean weight convergence (Sec. IV) -----------------------
Nrun = 300; Nit = 300; Wm = zeros(S.N, Nit+1); Pm = zeros(1, Nit+1);
for n = 1:Nrun
    rng(2000+n); Xn = L*randn(S.N, Nit);
    rr = frost1972_adapt(S, Xn, S.mu, 'clms', []);
    Wm = Wm + rr.W/Nrun; Pm = Pm + rr.Pout/Nrun;
end
dev = sqrt(sum((Wm - S.Wopt).^2, 1)); d0 = norm(S.F - S.Wopt); kk = 0:Nit;
h = figure('Color', 'w', 'Position', [100 100 900 450]); hold on; grid on;
plot(kk, dev/d0, '-', 'Color', cTeal, 'LineWidth', 2.2, 'DisplayName', sprintf('ensemble mean of %d runs', Nrun));
plot(kk, (1-S.mu*min(S.sigma)).^kk, '--', 'Color', cBlue, 'LineWidth', 1.4, 'DisplayName', '(1-\mu\sigma_{min})^k  slow bound');
plot(kk, (1-S.mu*max(S.sigma)).^kk, ':',  'Color', cBlue, 'LineWidth', 1.4, 'DisplayName', '(1-\mu\sigma_{max})^k  fast bound');
xlabel('Number of adaptations k'); ylabel('||E[W(k)] - W_{opt}|| / ||W(0) - W_{opt}||');
legend('Location', 'northeast'); ylim([0 1.05]);
title('Convergence of the mean weight vector (Sec. IV, eq. 24-25)');
savef(h, 'fig_mean_convergence.png');

%% ---- 5. eigenvalues of Rxx and P Rxx P ------------------------------------
h = figure('Color', 'w', 'Position', [100 100 900 420]);
hold on; grid on;
stem(1:S.N, S.lambda, 'filled', 'Color', cBlue, 'LineWidth', 1.5, 'DisplayName', 'eigenvalues of R_{xx}');
stem((1:S.N)+0.25, S.sigma_all, 'filled', 'Color', cTeal, 'LineWidth', 1.5, 'DisplayName', 'eigenvalues of P R_{xx} P');
xlabel('index (sorted)'); ylabel('eigenvalue'); legend('Location', 'northwest');
title(sprintf('\\lambda(R_{xx}) \\in [%.3f, %.3f]  (paper: 0.111 ... 8.355);  P R_{xx} P has J = %d zero eigenvalues', min(S.lambda), max(S.lambda), S.J));
savef(h, 'fig_eigenvalues.png');

%% ---- 6. summary -------------------------------------------------------------
fid = fopen(fullfile(outdir, 'results_summary.txt'), 'w');
pr = @(varargin) prf(fid, varargin{:});
pr('Frost (1972) reproduction  --  K = %d sensors, J = %d taps, KJ = %d weights\n', S.K, S.J, S.N);
pr('F^T = [%g %g %g %g]\n', S.Fvec);
pr('eigenvalues of Rxx      : %.4f ... %.4f    (paper 0.111 ... 8.355)\n', min(S.lambda), max(S.lambda));
pr('nonzero eig of P Rxx P  : %.4f ... %.4f\n', min(S.sigma), max(S.sigma));
pr('tr(Rxx) = %.3f   tr(P Rxx P) = %.3f\n', trace(S.Rxx), S.trPRP);
pr('mu = %.3f ;  1/sigma_max = %.4f ; bound (26) = %.4f (paper 0.074) ; bound (28) = %.4f\n', S.mu, S.mu_stab, S.mu_max26, S.mu_max28);
pr('misadjustment (27)      : %.1f %% ... %.1f %%   (paper 15.2 ... 17.0 %%)\n', 100*S.M_lo, 100*S.M_hi);
pr('time constants (25)     : fastest %.1f, slowest %.0f adaptations\n', S.tau_fast, S.tau_slow);
pr('optimum output power    : %.4f ; initial W(0)=F output power : %.4f\n', S.Popt, S.F'*S.Rxx*S.F);
pr('single run, last 50 of 200: CLMS mean power %.3f, GP mean power %.3f\n', mean(clms.Pout(end-49:end)), mean(gp.Pout(end-49:end)));
pr('ensemble (%d runs) mean output power at k=300: %.3f ; bounds %.3f ... %.3f\n', Nrun, Pm(end), bounds(2), bounds(3));
pr('6-digit constraint error, mean over 200 adaptations: CLMS %.2e, GP %.2e ; GP at k=2000: %.2e\n', mean(clms.cerr(2:end)), mean(gp.cerr(2:end)), gl.cerr(end));
pr('Wopt^T = [%s]\n', sprintf('%.4f ', S.Wopt));
fclose(fid);
save(fullfile(outdir, 'frost1972_results.mat'), 'S', 'clms', 'gp', 'Wm', 'Pm', 'dev', 'cl', 'gl');
close all;
fprintf('figures written to %s\n', outdir);

function prf(fid, varargin)
fprintf(1, varargin{:}); fprintf(fid, varargin{:});
end
