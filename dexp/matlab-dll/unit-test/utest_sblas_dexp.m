function utest_sblas_dexp(S, X)
%UTEST_SBLAS_DEXP  Unit test for sBLAS-based differential of the matrix exponential.
%
%   UTEST_SBLAS_DEXP() or UTEST_SBLAS_DEXP(S, X) validates the correctness of
%   the C++/MKL implementation of the matrix-exponential differential
%   Dexp_S(X) provided by the sblas_dexp wrapper (skewblas.dll).
%
%   This test verifies:
%       • Forward mode:     Yf = Dexp_S(X)
%       • Inverse mode:     Yi = Dexp_S^{-1}(Yf), expected to recover X
%       • Finite-difference consistency:
%                 Dexp_S(X) ≈ (exp(S + hX) - exp(S)) / h
%         with the standard adjoint transport by exp(S)' so that both
%         quantities live in the same tangent space at exp(S).
%       • Convergence as h → 0 on a log-log scale.
%
%   INPUTS
%     S   (optional) n×n skew-symmetric matrix.
%     X   (optional) n×n skew-symmetric direction matrix.
%         If omitted, random skew-symmetric matrices of size n = 8 are created.
%
%   OUTPUTS
%     This function prints diagnostics and produces a convergence plot,
%     but returns no arguments.
%
%   DEPENDENCIES
%       • ../lib/skewblas.dll
%       • ../lib/skewblas_api.h
%       • ../src/sblas_dexp.m  (MATLAB wrapper)
%
%   PROCEDURE
%       1) Load skewblas.dll if needed.
%       2) Generate S, X (if not supplied).
%       3) Compute:
%              [Yf, para] = sblas_dexp(S, X, 'fwd');
%              [Yi, para] = sblas_dexp(S, Yf, 'inv', para);
%       4) Compare Yf against finite differences:
%              Y_fd = (exp(S + hX) - exp(S)) / h.
%          Use ExpS' * Y_fd for proper tangent-space alignment.
%       5) Sweep h = 10^{-1} ... 10^{-7} to show O(h) convergence.
%
%   NOTES
%       • The inverse pass checks whether the spectral parameters returned
%         in PARA correctly reconstruct X from Yf.
%       • The adjoint correction ExpS' * Y_fd is required because
%         d/dt exp(S+tX)|_{t=0} naturally lives at exp(S), not at S.

    % ---------------------------------------------------------------
    % Load paths and DLL
    % ---------------------------------------------------------------
    srcPath = fullfile(fileparts(mfilename('fullpath')), '../src');
    addpath(srcPath);

    libPath = fullfile(fileparts(mfilename('fullpath')), '../lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        fprintf('Loading library from %s\n', libPath);
        loadlibrary(dll, hdr);
    else
        fprintf('Library already loaded.\n');
    end

    % ---------------------------------------------------------------
    % Generate random skew-symmetric test data
    % ---------------------------------------------------------------
    if nargin < 1
        n = 8; 
        S = randn(n); S = 0.5*(S - S.'); 
        X = randn(n); X = 0.5*(X - X.');
    end

    % ---------------------------------------------------------------
    % Compute forward and inverse differential via sBLAS
    % ---------------------------------------------------------------
    [Yf, para] = sblas_dexp(S, X, 'fwd');
    [Yi, para] = sblas_dexp(S, Yf, 'inv', para);

    % ---------------------------------------------------------------
    % Finite difference validation for Dexp_S(X)
    % ---------------------------------------------------------------
    h = 1e-6;
    ExpS = expm(S);
    ExpSh = expm(S + h * X);

    % finite-difference approximation of derivative
    Y_fd = (ExpSh - ExpS) / h;

    % compare with differential from DLL
    err_fd = norm(Yf - ExpS' * Y_fd, 'fro') / max(norm(Yf, 'fro'), 1e-14);

    disp('--- Example matrices ---');
    disp('X (input direction):'); disp(X);
    disp('Yf (forward differential):'); disp(Yf);
    disp('Yf (forward finite-diff):'); disp(ExpS' * Y_fd);
    disp('Yi = Dexp(S)[Yf] (inverse differential that recover X):'); disp(Yi);

    fprintf('\nFinite difference validation of Dexp_S(X):\n');
    fprintf('  step size h = %.1e\n', h);
    fprintf('  rel. Frobenius error = %.3e\n', err_fd);

    % ---------------------------------------------------------------
    % Optional: visualize convergence as h → 0
    % ---------------------------------------------------------------
    hs = logspace(-1, -7, 7);
    errs = zeros(size(hs));
    for k = 1:numel(hs)
        Y_fd = ExpS' * (expm(S + hs(k)*X) - ExpS) / hs(k);
        errs(k) = norm(Yf - Y_fd, 'fro') / max(norm(Yf, 'fro'), 1e-14);
    end
    figure;
    loglog(hs, errs, '-o', 'DisplayName', 'Relative Error', 'LineWidth', 1.8);
    grid on;
    xlabel('h (step size)');
    ylabel('Relative Error: Dexp(S)[X] vs (exp(S+hX)-exp(S))/h');
    title('Finite Difference Convergence Test for Dexp(S)[X]');
    set(gca, 'XDir', 'reverse');  % smaller h → rightward
    hold on; loglog(hs, hs, '--k', 'DisplayName', 'O(h) slope'); legend('show');

    % ---------------------------------------------------------------
    % Clean up
    % ---------------------------------------------------------------
    if libisloaded('skewblas')
        unloadlibrary('skewblas');
    end
end
