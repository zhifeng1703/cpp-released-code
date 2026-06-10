function utest_sblas_dexp(n)
% UTEST_SBLAS_DEXP  
%   Validate the forward/inverse differential of the matrix exponential
%   using the skewblas MEX interface.
%
%   utest_sblas_dexp(n)
%
%   Inputs:
%       n - matrix size (default = 8)
%
%   This test performs:
%       1) Generate a random skew-symmetric S.
%       2) Compute Schur data (R, A) via mex_sblas_ss2schur.
%       3) Compute DK parameters via mex_sblas_dexp_para(A).
%       4) Evaluate forward and inverse differential using mex_sblas_dexp.
%       5) Verify inverse consistency: Dexp^{-1}(Dexp(X)) = X.
%       6) Validate forward differential using finite differences.
%       7) Plot convergence of FD errors as h→0.
%
%   Requires MEX:
%       mex_sblas_ss2schur
%       mex_sblas_dexp_para
%       mex_sblas_dexp
%
%   Example:
%       utest_sblas_dexp(10);

    if nargin < 1, n = 8; end
    fprintf('==== SkewBLAS Dexp Test (n = %d) ====\n', n);

    % ============================================================
    % 0. Environment setup: add mex and build if missing
    % ============================================================
    thisFile = mfilename('fullpath');
    rootDir  = fileparts(fileparts(thisFile));       % .../skewblas
    mexDir   = fullfile(rootDir, 'mex');

    if ~contains(path, mexDir)
        addpath(mexDir);
        fprintf('[utest] Added path: %s\n', mexDir);
    end

    % required mex files
    mex_required = { ...
        'mex_sblas_ss2schur', ...
        'mex_sblas_dexp_para', ...
        'mex_sblas_dexp' ...
    };

    need_build = false;
    for k = 1:numel(mex_required)
        if exist(mex_required{k}, 'file') ~= 3
            fprintf('[utest] Missing MEX: %s — will rebuild.\n', mex_required{k});
            need_build = true;
        end
    end

    if need_build
        buildScript = fullfile(mexDir, 'build_mex.m');
        if ~isfile(buildScript)
            error('build_mex.m not found in: %s', mexDir);
        end
        fprintf('[utest] Building MEX files...\n');
        run(buildScript);
    end


    % ============================================================
    % 1. Construct skew-symmetric matrix S
    % ============================================================
    S = randn(n);  
    S = 0.5 * (S - S');     % enforce skew-symmetric

    % ============================================================
    % 2. Schur decomposition S = R * skew(A) * R'
    % ============================================================
    [R, A] = mex_sblas_ss2schur(S);
    m = length(A);

    % Construct the block-skew matrix D(A)
    DA = zeros(n);
    for i = 1:m
        DA(2*i-1, 2*i) = -A(i);
        DA(2*i,   2*i-1) =  A(i);
    end

    % ============================================================
    % 3. Compute DK parameters (forward & inverse)
    % ============================================================
    [ParaF, ParaI] = mex_sblas_dexp_para(A);

    para = struct('R', R, 'A', A, 'Fwd', ParaF, 'Inv', ParaI);

    % ============================================================
    % 4. Prepare random test direction X
    % ============================================================
    X = randn(n);  
    X = 0.5*(X - X');

    % ============================================================
    % 5. Forward & inverse differential via MEX
    % ============================================================
    Yf = mex_sblas_dexp(S, X, 'fwd', para);
    Yi = mex_sblas_dexp(S, Yf, 'inv', para);

    % ============================================================
    % 6. Finite-difference validation for Dexp_S(X)
    % ============================================================
    h = 1e-6;
    ExpS  = expm(S);
    ExpSh = expm(S + h*X);

    Y_fd = (ExpSh - ExpS) / h;

    % Theoretical adjoint correction:
    Y_fd_adj = ExpS' * Y_fd;

    err_fd = norm(Yf - Y_fd_adj, 'fro') / max(norm(Yf,'fro'), 1e-14);
    

    disp('--- SkewBLAS Differential Test ---');
    disp('X (input direction):'); disp(X);
    disp('Yf (forward differential):'); disp(Yf);
    disp('Yf_fd (forward finite-diff):'); disp(ExpS' * Y_fd);
    disp('Yi (inverse differential, should recover X):'); disp(Yi);
 

    fprintf('\nFinite difference validation for DK Dexp_S(X):\n');
    fprintf('  step size h = %.1e\n', h);
    fprintf('  rel. Frobenius error = %.3e\n', err_fd);

    % ============================================================
    % 7. Convergence test as h→0
    % ============================================================
    hs = logspace(-1, -7, 7);
    errs = zeros(size(hs));
    for kk = 1:numel(hs)
        Y_temp = ExpS' * (expm(S + hs(kk)*X) - ExpS) / hs(kk);
        errs(kk) = norm(Yf - Y_temp, 'fro') / max(norm(Yf,'fro'), 1e-14);
    end

    figure('Name','SBLAS Dexp finite-difference convergence','Color','w');
    loglog(hs, errs, '-o', 'LineWidth', 1.8);
    hold on; grid on;
    loglog(hs, hs, '--k', 'DisplayName', 'O(h) slope');
    xlabel('h'); ylabel('Relative Error');
    title(sprintf('Finite Difference Convergence for SBLAS Dexp, n = %d', n));
    set(gca,'XDir','reverse');
    legend('Error','O(h)', 'Location','best');

end
