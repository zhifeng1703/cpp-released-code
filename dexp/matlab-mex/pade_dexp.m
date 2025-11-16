function Y = pade_dexp(S, X, dir, p)
% PADE_DEXP  Differential of exp(S) using Padé approximations (MEX backend).
%
%   Y = pade_dexp(S, X, dir, p)
%
%   Inputs:
%     S   : n×n real SKEW-SYMMETRIC matrix
%     X   : n×n real skew-symmetric matrix (direction)
%     dir : 'fwd' or 'inv'
%     p   : Padé order (default: 13 for 'fwd', 7 for 'inv')
%
%   Output:
%     Y : n×n real matrix, differential result
%
%   Notes:
%     - This version uses compiled MEX gateways:
%           mex_pade_dexpfwd(S, X, p)
%           mex_pade_dexpinv(S, X, p)
%     - No DLL / loadlibrary is used anymore.
%

    % ---------------------------------------------------------------
    % Defaults & argument handling
    % ---------------------------------------------------------------
    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end

    if nargin < 4 || isempty(p)
        if startsWith(lower(dir), 'i')
            p = 7;   % good inverse default
        else
            p = 13;  % good forward default
        end
    end

    % ---------------------------------------------------------------
    % Ensure skew-symmetric inputs (for safety)
    % ---------------------------------------------------------------
    S = 0.5 * (S - S.');
    X = 0.5 * (X - X.');
    

    % ---------------------------------------------------------------
    % Dispatch to MEX gateways
    % ---------------------------------------------------------------
    if startsWith(lower(dir), 'f')      % ---------- FORWARD
        Y = mex_pade_dexpfwd(S, X, p);
    else                                % ---------- INVERSE
        [R, A] = mex_sblas_so2schur(S);
        s = compute_s_min(A, 2.88e-1);
        Y = mex_pade_dexpinv(R, A, X, p, s);
    end
end

function s_min = compute_s_min(theta, delta)
% COMPUTE_S_MIN  Compute minimal s such that rho(D^(1/2^s) - I) <= delta.
%
%   theta : vector of eigen-angles in [0,pi] (only the m positive ones)
%   delta : tolerance (0 < delta <= 2)
%
%   s_min : smallest integer s >= 0 satisfying
%           2*sin(theta_max/2^(s+1)) <= delta.

    % if any(theta < -pi) || any(theta > pi)
    %     error('Angles theta must lie in [-pi, pi].');
    % end

    % -------------------------------------------------------------
    % Input check
    % -------------------------------------------------------------

    
    if ~(delta > 0 && delta <= 2)
        error('delta must lie in (0, 2].');
    end

    theta_max = max(abs(theta));  % automatically handles zero-case

    if theta_max == 0
        s_min = 0;
        return;
    end

    rhs = delta / 2;
    alpha = asin(rhs);  % arcsin(delta/2)

    s_raw = log2(theta_max / alpha) - 1;

    % -------------------------------------------------------------
    % Minimal integer s >= 0
    % -------------------------------------------------------------
    s_min = max(0, ceil(s_raw));
end

