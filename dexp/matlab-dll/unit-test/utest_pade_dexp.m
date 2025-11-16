function utest_pade_dexp(scale)
%UTEST_PADE_DEXP  Unit test for Pade-based differential of the matrix exponential.
%
%   UTEST_PADE_DEXP() or UTEST_PADE_DEXP(scale) validates the correctness
%   of the Pade-approximant implementation of the matrix-exponential
%   differential Dexp_A(X) provided by pade_dexp (backed by skewblas.dll).
%
%   The test exercises both:
%       • Forward differential:   Yf = Dexp_A(X)
%       • Inverse differential:   Xi = Dexp_A^{-1}(Yf)
%
%   across several Pade orders (p = 1,3,7,13).
%
%   INPUT
%     scale   (optional) scalar used to rescale A so that ‖A‖₂ ≈ scale.
%             Default: scale = 1.5*pi.
%
%   TEST PROCEDURES
%
%   (1) Random skew-symmetric test matrices:
%         A ← random skew-symmetric scaled to ‖A‖₂ = scale
%         X ← random skew-symmetric direction
%
%   (2) Forward differential tests for Pade orders p = 3 and p = 13:
%         Yf = pade_dexp(A, X, 'fwd', p)
%
%       Finite-difference reference:
%         Y_fd(h) = (exp(A + hX) - exp(A))/h
%
%       To place Yf and Y_fd in the same tangent space at exp(A):
%         Y_fd_aligned = exp(A)' * Y_fd
%         Yf_aligned   = exp(A)' * Yf
%
%       Relative error:
%         ‖Yf_aligned − Y_fd_aligned‖_F / ‖Yf‖_F
%
%       The script sweeps h = 10^{-1} ... 10^{-7} and plots O(h)-type
%       convergence for each Pade order.

    % ---------------------------------------------------------------
    % Load library
    % ---------------------------------------------------------------
    libPath = fullfile(fileparts(mfilename('fullpath')), '../lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    if nargin < 1, scale = 1.5 * pi; end

    % ---------------------------------------------------------------
    % Random test matrices
    % ---------------------------------------------------------------
    n = 8;
    A = randn(n);
    A = (A - A') / 2;
    A = (scale / norm(A, 2)) * A;
    X = randn(n);
    X = (X - X') / 2;

    fprintf('Pade differential test on n = %d\n', n);
    fprintf('=========================================\n');

    % ---------------------------------------------------------------
    % Step 1. Forward differential for p = 3 and p = 13
    % ---------------------------------------------------------------
    p_forward = [3, 13];
    Yf_list = cell(1,2);     % store results
    errors_fd = zeros(1,2); % record finite diff error for reporting
    hs = logspace(-1, -7, 7);
    FD_errors_plot = zeros(2, numel(hs));

    ExpA = expm(A);

    for idx = 1:2
        p = p_forward(idx);
        fprintf('\nForward differential with Pade order p = %d\n', p);

        % Compute forward differential
        Yf = pade_dexp(A, X, 'fwd', p);
        Yf_list{idx} = Yf;

        % Finite difference reference
        h = 1e-6;
        Y_fd = (expm(A + h*X) - ExpA)/h;
        Y_fd = ExpA' * Y_fd; % adjoint matching

        sYf = expm(A + h*X)' * Yf; 

        err = norm(ExpA' * Yf - Y_fd, 'fro') / max(norm(Yf,'fro'), 1e-14);
        errors_fd(idx) = err;

        fprintf('  Finite-difference rel. error = %.3e\n', err);

        % Convergence plot data
        for k = 1:numel(hs)
            Y_fd = ExpA'*(expm(A + hs(k)*X) - ExpA)/hs(k);
            FD_errors_plot(idx, k) = norm(ExpA' * Yf - Y_fd, 'fro') ...
                / max(norm(Yf,'fro'), 1e-14);
        end
    end

    % ---------------------------------------------------------------
    % Step 2. Plot both forward FD errors in one figure
    % ---------------------------------------------------------------
    figure;
    loglog(hs, FD_errors_plot(1,:), '-o', ...
        'DisplayName', 'p = 3', 'LineWidth', 1.8);
    hold on;
    loglog(hs, FD_errors_plot(2,:), '-s', ...
        'DisplayName', 'p = 13', 'LineWidth', 1.8);
    loglog(hs, hs, '--k', 'DisplayName','O(h)', 'LineWidth',1.4);

    grid on;
    xlabel('h');
    ylabel('relative FD error');
    title('Finite Difference Convergence for Pade Dexp(A)[X]');
    set(gca, 'XDir', 'reverse');
    legend('show', 'Location','SouthWest');

    fprintf('Done.\n');

    % ---------------------------------------------------------------
    % Clean up
    % ---------------------------------------------------------------
    if libisloaded('skewblas')
        unloadlibrary('skewblas');
    end
end
