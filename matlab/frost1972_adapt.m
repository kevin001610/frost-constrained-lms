function out = frost1972_adapt(S, X, mu, algo, ndigits)
%FROST1972_ADAPT  Run the constrained LMS (eq. 22) or gradient-projection (eq. 31) algorithm.
%
%   out = FROST1972_ADAPT(S, X, mu, algo, ndigits)
%     S       : struct from frost1972_setup
%     X       : KJ x Nsnap tap-voltage snapshots
%     mu      : convergence constant
%     algo    : 'clms'  ->  W(k+1) = P[W(k) - mu y(k) X(k)] + F          (22)
%               'gp'    ->  W(k+1) =  W(k) - mu P y(k) X(k)              (31)
%     ndigits : [] for exact double arithmetic, otherwise every stored
%               quantity is rounded to ndigits significant decimal digits
%               (the paper used 6-digit floating point on an HP-2116).
%
%   out.W      : KJ x (Nsnap+1) weight history, W(:,1) = F = W(0)
%   out.Pout   : 1 x (Nsnap+1) true output power W(k)' Rxx W(k)   (as in Fig. 11/12)
%   out.cerr   : 1 x (Nsnap+1) constraint error ||C'W(k) - Fvec||^2
%   out.y      : 1 x Nsnap array output y(k) = W(k)' X(k)

Nsnap = size(X, 2);
q = @(v) v;
if ~isempty(ndigits)
    q = @(v) quant(v, ndigits);
end
W = q(S.F);                             % W(0) = F = C (C'C)^{-1} Fvec
Wh = zeros(S.N, Nsnap+1); Wh(:,1) = W;
y  = zeros(1, Nsnap);
for k = 1:Nsnap
    x = X(:, k);
    y(k) = q(W' * x);                   % array output
    g = q(mu * y(k) * x);               % mu * (unconstrained gradient estimate)
    switch lower(algo)
        case 'clms'
            W = q(q(S.P * q(W - g)) + S.F);   % eq. 22: project, then add F
        case 'gp'
            W = q(W - q(S.P * g));             % eq. 31: add projected gradient only
        otherwise
            error('unknown algorithm %s', algo);
    end
    Wh(:, k+1) = W;
end
out.W    = Wh;
out.y    = y;
out.Pout = sum(Wh .* (S.Rxx * Wh), 1);
E        = S.C' * Wh - S.Fvec;
out.cerr = sum(E.^2, 1);
end

function v = quant(v, nd)
% round to nd significant decimal digits (emulates short-mantissa floating point)
nz = v ~= 0;
e  = floor(log10(abs(v(nz))));
sc = 10.^(nd - 1 - e);
v(nz) = round(v(nz) .* sc) ./ sc;
end
