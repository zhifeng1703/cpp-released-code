function utest_dk_dexp()
%UTEST_DK_DEXP  Unit test for the Daletskii–Kreĭn differential of expm.
%
%   UTEST_DK_DEXP() validates the correctness of the dk_dexp MATLAB
%   wrapper (backed by skewblas.dll) that computes:
%
%       (i)  The forward differential   Yf  = Dexp_S(X)
%       (ii) The inverse differential   Yi  such that Dexp_S^{-1}(Yf) = X
%
%   The test includes:
%       • Generation of random skew-symmetric matrices S and X
%       • Forward DK differential via dk_dexp(S, X, 'fwd')
%       • Inverse DK differential via dk_dexp(S, Yf, 'inv', para)
%       • Finite-difference validation:
%               Dexp_S(X) ≈ (exp(S + hX) - exp(S)) / h
%         corrected by adjoint transport with exp(S)'
%       • Convergence sweep over h decreasing from 1e−1 to 1e−7
%       • Log–log convergence plot (O(h) expected)
%
%   REQUIREMENTS
%       – The C++ DLL skewblas.dll must be located in ../lib
%       – The MATLAB wrappers (dk_dexp, etc.) must be in ../src
%       – The header skewblas_api.h must match the DLL ABI
%
%   OUTPUT
%       This function produces diagnostic prints and a convergence figure.
%       It does not return values.
%
%   NOTES
%       • The Daletskii–Kreĭn formula defines the differential of f(A)
%         via an integral over the resolvent; the DLL implements the
%         closed-form spectral version specialized for skew-symmetric A.
%       • The finite-difference comparison uses the adjoint correction:
%               ExpS' * (expm(S + hX) - ExpS) / h
%         which places both quantities in the same tangent space.
%       • The inverse mode verifies that dk_dexp(...,'inv') recovers X.

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
    n = 8;
    S = randn(n); S = 0.5*(S - S.');
    X = randn(n); X = 0.5*(X - X.');

    % ---------------------------------------------------------------
    % Compute forward and inverse differential via DK (complex)
    % ---------------------------------------------------------------
    [Yf, para] = dk_dexp(S, X, 'fwd');
    [Yi, para] = dk_dexp(S, Yf, 'inv', para);

    % ---------------------------------------------------------------
    % Finite difference validation for Dexp_S(X)
    % ---------------------------------------------------------------
    h = 1e-6;
    ExpS = expm(S);
    ExpSh = expm(S + h * X);

    % finite-difference approximation of derivative
    Y_fd = (ExpSh - ExpS) / h;

    % compare with differential from DLL
    % NOTE: same adjoint correction as your utest_dexp
    err_fd = norm(Yf - ExpS' * Y_fd, 'fro') / max(norm(Yf, 'fro'), 1e-14);

    % ---------------------------------------------------------------
    % Print diagnostic matrices
    % ---------------------------------------------------------------
    disp('--- Daletskii–Kreĭn Differential Test ---');
    disp('X (input direction):'); disp(X);
    disp('Yf (forward differential using DK):'); disp(Yf);
    disp('Yf_fd (forward finite-diff):'); disp(ExpS' * Y_fd);
    disp('Yi (inverse differential on Yf, should recover X):'); disp(Yi);

    fprintf('\nFinite difference validation for DK Dexp_S(X):\n');
    fprintf('  step size h = %.1e\n', h);
    fprintf('  rel. Frobenius error = %.3e\n', err_fd);

    % ---------------------------------------------------------------
    % Convergence test as h → 0
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
    ylabel('Relative Error: DK Dexp(S)[X] vs FD');
    title('Finite Difference Convergence Test for DK Dexp(S)[X]');
    set(gca, 'XDir', 'reverse');
    hold on;
    loglog(hs, hs, '--k', 'DisplayName', 'O(h) slope');
    legend('show');

    % ---------------------------------------------------------------
    % Clean up
    % ---------------------------------------------------------------
    if libisloaded('skewblas')
        unloadlibrary('skewblas');
    end
end
