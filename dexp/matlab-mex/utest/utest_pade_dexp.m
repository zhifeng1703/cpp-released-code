function utest_pade_dexp(n, scale)
% UTEST_PADE_DEXP
%   Validate the forward/inverse differential of the matrix exponential
%   using the MATLAB Padé approximation routine pade_dexp.m.
%
%   utest_pade_dexp(n)
%
%   Inputs:
%       n - matrix size (default = 8)
%
%   Tests performed:
%       1) Generate random skew-symmetric S and X.
%       2) Evaluate Dexp_S(X) with Padé orders p = 3 and 13.
%       3) Verify inverse Dexp^{-1}_S(Y) ≈ X using p = 1 and 7.
%       4) Validate forward differential using finite differences.
%       5) Check convergence of FD error as h → 0.
%
%   Requires:
%       pade_dexp.m
%
%   Example:
%       utest_pade_dexp(10);

    if nargin < 1, n = 8; end
    if nargin < 2, scale = 4.0;end
    fprintf('==== Padé Dexp Test (n = %d) ====\n', n);

    % ============================================================
    % 0. Add ./src to path if missing
    % ============================================================
    thisFile = mfilename('fullpath');
    rootDir  = fileparts(fileparts(thisFile));       % .../skewblas
    srcDir   = fullfile(rootDir, 'src');

    if ~contains(path, srcDir)
        addpath(srcDir);
        fprintf('[utest] Added path: %s\n', srcDir);
    end

    % ============================================================
    % 1. Construct skew-symmetric matrix S and direction X
    % ============================================================
    S = 2 * rand(n) - 1;  S = 0.5*(S - S');
    S = (scale / norm(S,2)) * S;

    X = 2 * rand(n) - 1;  X = 0.5*(X - X');

    fprintf('[utest] Generated random skew-symmetric S and X.\n');

    h = 1e-6;
    ExpS  = expm(S);
    ExpSh = expm(S + h*X);

    Y_fd = (ExpSh - ExpS) / h;
    Y_fd_adj = ExpS' * Y_fd;

    % ============================================================
    % 2. Forward Padé Dexp  (p = 3, 13)
    % ============================================================
    fprintf('\n[Forward Direction Tests]\n');
    pf_list = [3, 13];

    Yf = cell(size(pf_list));
    for k = 1:numel(pf_list)
        p = pf_list(k);
        fprintf('  - Evaluating forward Padé order p = %d ...\n', p);
        Yf{k} = mex_pade_dexpfwd(S, X, p);
        err_fd = norm(Yf{k} - Y_fd, 'fro') / max(norm(Y_fd,'fro'), 1e-14);

        fprintf('  relative FD error = %.3e\n', err_fd);
    end

    % ============================================================
    % 3. Inverse Padé Dexp  (p = 1, 7)
    % ============================================================
    fprintf('\n[Inverse Direction Tests]\n');
    pi_list = [1, 7];

    [R, A] = mex_sblas_so2schur(expm(S));
    s = compute_s_min(A, 2.88e-1);

    Yi = cell(size(pi_list));
    for k = 1:numel(pi_list)
        p = pi_list(k);
        fprintf('  - Evaluating inverse Padé order p = %d ...\n', p);
        Yi{k} = mex_pade_dexpinv(R, A, Yf{2}, p, s);    % invert the high-order forward
        rel_err_inv = norm(Yi{k} - X, 'fro') / max(norm(X,'fro'), 1e-14);
        fprintf('      relative inverse error = %.3e\n', rel_err_inv);
    end


    err_fd = norm(Yf{2} - Y_fd_adj, 'fro') / max(norm(Yf{2},'fro'), 1e-14);

    fprintf('  h = %.1e\n', h);
    fprintf('  relative FD error = %.3e\n', err_fd);


    disp('--- Pade Differential Test ---');
    disp('X (input direction):'); disp(X);
    disp('Yf-[3/3] (forward differential):'); disp(ExpS' * Yf{1});
    disp('Yf-[13/13] (forward differential):'); disp(ExpS' * Yf{2});
    disp('Yf_fd (forward finite-diff):'); disp(ExpS' * Y_fd);
    disp('Yi-[1/1] (inverse differential, recover X if norm(S, 2) < \pi):'); disp(Yi{1});
    disp('Yi-[7/7] (inverse differential, recover X if norm(S, 2) < \pi):'); disp(Yi{2});

    % ============================================================
    % 5. Convergence study
    % ============================================================
    % fprintf('\n[FD Convergence Test]\n');
    % hs   = logspace(-1, -7, 7);
    % errs = zeros(size(hs));
    % 
    % for kk = 1:numel(hs)
    %     hk = hs(kk);
    %     Y_temp = ExpS' * (expm(S + hk*X) - ExpS) / hk;
    %     errs(kk) = norm(Yf{2} - Y_temp, 'fro') / max(norm(Yf{2},'fro'), 1e-14);
    % end
    % 
    % figure('Name','Padé Dexp finite-difference convergence','Color','w');
    % loglog(hs, errs, '-o', 'LineWidth', 1.8);
    % grid on; hold on;
    % loglog(hs, hs, '--k', 'DisplayName','O(h)');
    % xlabel('h'); ylabel('Relative Error');
    % title(sprintf('Finite Difference Convergence for Padé Dexp (p = 13), n = %d', n));
    % set(gca, 'XDir','reverse');
    % legend('Error','O(h)', 'Location','best');

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
