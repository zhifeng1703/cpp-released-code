function debug_pade_dexp(scale)
%UTEST_PADE_DEXP  Unit test for Padé-based differential of the matrix exponential.
%
%   UTEST_PADE_DEXP() or UTEST_PADE_DEXP(scale) validates the correctness
%   of the Padé-approximant implementation of the matrix-exponential
%   differential Dexp_S(X) using the compiled MEX routines:
%
%       Y = mex_pade_dexpfwd(S, X, p)
%       Y = mex_pade_dexpinv(S, X, p)
%
%   The test exercises:
%       • Forward differential: Yf = Dexp_S(X)
%       • Convergence of Yf vs finite difference reference
%
%   INPUT
%     scale  (optional) target value for ‖S‖₂. Default: 1.5*pi.
%

    if nargin < 1
        scale = 1.5*pi;
    end

    % ---------------------------------------------------------------
    % Random test matrices
    % ---------------------------------------------------------------
    n = 4;
    S = randn(n);  S = 0.5*(S - S');
    S = (scale / norm(S,2)) * S;

    X = randn(n);  X = 0.5*(X - X');

    ExpS = expm(S);

    fprintf('Padé Dexp unit test (n = %d)\n', n);
    fprintf('============================================\n');

    % ---------------------------------------------------------------
    % Padé forward tests p = 3, 13
    % ---------------------------------------------------------------
    p_list = [3, 13];
    hs = logspace(-1, -7, 7);
    FD_errors = zeros(numel(p_list), numel(hs));

    for ip = 1:numel(p_list)
        p = p_list(ip);
        fprintf('\nForward differential with Padé order p = %d\n', p);

        % --- Forward differential ---
        Yf = mex_pade_dexpfwd(S, X, p);

        % FD reference for h = 1e-6 first
        h = 1e-6;
        Y_fd = (expm(S + h*X) - ExpS) / h;
        Y_fd = ExpS' * Y_fd;      % adjoint-transported FD
        Yf_aligned = ExpS' * Yf;  % same tangent space

        rel_err = norm(Yf_aligned - Y_fd, 'fro') / max(norm(Yf_aligned,'fro'),1e-14);
        fprintf('  FD rel.error (h = 1e-6) = %.3e\n', rel_err);

        % --- Convergence sweep ---
        for ih = 1:numel(hs)
            h = hs(ih);
            Y_fd = (expm(S + h*X) - ExpS) / h;
            Y_fd = ExpS' * Y_fd;
            FD_errors(ip, ih) = norm(Yf_aligned - Y_fd, 'fro') ...
                                / max(norm(Yf_aligned,'fro'),1e-14);
        end
    end

    % ---------------------------------------------------------------
    % Plot convergence curves
    % ---------------------------------------------------------------
    figure; hold on;
    % loglog(hs, FD_errors(1,:), '-o', 'LineWidth',1.8, ...
    %     'DisplayName','p = 3');
    loglog(hs, FD_errors(2,:), '-s', 'LineWidth',1.8, ...
        'DisplayName','p = 13');
    loglog(hs, hs, '--k', 'LineWidth',1.4, ...
        'DisplayName','O(h)');

    set(gca,'XDir','reverse');
    grid on;
    xlabel('h');
    ylabel('Relative FD error');
    title('Finite Difference Convergence of Padé Dexp(S)[X]');
    legend('Location','SouthWest');
end
