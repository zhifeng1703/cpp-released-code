function test_logm(dims, opts)
%TEST_LOGM  Benchmark and accuracy comparison of four logarithm implementations.
%
%   TEST_LOGM(dims)
%
%   This routine performs a timing and reconstruction-error experiment for
%   the following logarithm-map implementations on special orthogonal matrices:
%
%       (1) built-in         – real(logm(Q))
%       (2) mex_sblas_logm   – principal logarithm
%       (3) mex_sblas_diffeo_logm(Q,S)
%       (4) mex_sblas_spdiff_logm(Q,S)
%
%   For each dimension n in dims:
%       • Generate a random skew-symmetric T0 and define
%             S = mex_sblas_logm(expm(T0))
%         so that S lies in the special component.
%       • Generate a random skew-symmetric T and define Q = expm(T).
%       • Evaluate each logm implementation `loop = 20` times.
%       • Average the elapsed times over the 20 runs.
%       • Record the minimum observed reconstruction error
%             ||Q - expm(X)||_F
%         over all repetitions.
%
%   INPUT
%     dims   Vector of matrix sizes to test.
%            Default: 3:100.
%
%     opts   Struct with optional fields:
%              opts.verbose : true/false
%              opts.savefig : true/false
%
%   OUTPUT
%     This driver mainly prints progress and produces figures.
%
%   See also:
%       test_expm, mex_sblas_logm, mex_sblas_diffeo_logm, mex_sblas_spdiff_logm

    if nargin < 1, dims = 3:100; end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'savefig'), opts.savefig = true; end

    fprintf('=============LOGM experiment begins=============\n\n');

    test_logm_data(dims, opts);

    fprintf('=============LOGM experiment ends=============\n\n');

    fprintf(['IMPORTANT NOTE: The compute time of the first execution ' ...
        'could be unreliable even WITH the warm-up run.\n']);
    fprintf(['IMPORTANT NOTE: Make sure to run multiple calls to test_logm for ' ...
        'more reliable timings.\n']);
end


function [elapsed, errors, fig_elapsed, fig_error] = test_logm_data(dims, opts)
    if nargin < 1, dims = 3:100; end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'savefig'), opts.savefig = true; end

    addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'mex'));
    compile_mex();

    nd = length(dims);
    loop = 20;

    elapsed = zeros(nd, 4);
    errors  = +inf(nd, 4);

    idx = 1;
    for n = dims
        % ------------------------------------------------------------
        % Build S in the special component:
        %   S = mex_sblas_logm(expm(T0))
        % ------------------------------------------------------------
        T0 = 2 * rand(n) - 1;
        T0 = 0.5 * (T0 - T0');
        %Qs = expm(T0);
        %S  = mex_sblas_logm(Qs);
        S = (pi / norm(T0, 2)) * T0;
        % ------------------------------------------------------------
        % Build random special orthogonal target Q
        % ------------------------------------------------------------
        T = 2 * rand(n) - 1;
        T = 0.5 * (T - T');
        Q = expm(T);

        elap_err_logm(n, Q, S, struct('verbose', true));

        for k = 1:loop
            [t, e] = elap_err_logm(n, Q, S, struct('verbose', false));

            elapsed(idx, :) = elapsed(idx, :) + t(:)';
            errors(idx, :)  = min(errors(idx, :), e(:)');
        end

        elapsed(idx, :) = elapsed(idx, :) / loop;
        idx = idx + 1;
    end

    [fig_elapsed, fig_error] = logm_fig(dims, elapsed, errors, opts);
end


function [elapsed, errors] = elap_err_logm(n, Q, S, opt)
    if nargin < 1, n = 9; end
    if nargin < 2
        T = rand(n); T = 0.5 * (T - T');
        Q = expm(T);
    end
    if nargin < 3
        T0 = rand(n); T0 = 0.5 * (T0 - T0');
        S = (pi / norm(T0, 2)) * T0;
    end
    if nargin < 4
        opt = struct();
        opt.verbose = true;
    end
    if ~isfield(opt, 'verbose'), opt.verbose = true; end

    elapsed = zeros(4, 1);
    errors  = zeros(4, 1);

    % ------------------------------------------------------------
    % (1) built-in real(logm(Q))
    % ------------------------------------------------------------
    tic;
    X_builtin = logm(Q);
    elapsed(1) = toc;
    X_builtin = real(X_builtin - X_builtin')/2;

    % ------------------------------------------------------------
    % (2) principal logarithm
    % ------------------------------------------------------------
    tic;
    X_principal = mex_sblas_logm(Q);
    elapsed(2) = toc;

    % ------------------------------------------------------------
    % (3) diffeomorphic logarithm at S
    % ------------------------------------------------------------
    tic;
    X_diffeo = mex_sblas_diffeo_logm(Q, S);
    elapsed(3) = toc;

    % ------------------------------------------------------------
    % (4) special diffeomorphic logarithm at S
    % ------------------------------------------------------------
    tic;
    X_spdiff = mex_sblas_spdiff_logm(Q, S);
    elapsed(4) = toc;

    % ------------------------------------------------------------
    % Reconstruction errors: ||Q - expm(X)||_F
    % ------------------------------------------------------------
    errors(1) = norm(Q - expm(X_builtin),   'fro');
    errors(2) = norm(Q - expm(X_principal), 'fro');
    errors(3) = norm(Q - expm(X_diffeo),    'fro');
    errors(4) = norm(Q - expm(X_spdiff),    'fro');

    if opt.verbose
        fprintf('Testing matrix logarithms with size n = %d\n', n);
        fprintf(['‖Q - expm(X_builtin)‖_F = %.3e, ' ...
                 '‖Q - expm(X_principal)‖_F = %.3e, ' ...
                 '‖Q - expm(X_diffeo)‖_F = %.3e, ' ...
                 '‖Q - expm(X_spdiff)‖_F = %.3e\n'], ...
                 errors(1), errors(2), errors(3), errors(4));
        fprintf(['Elapsed time: built-in logm: %.2f us,\t ' ...
                 'principal: %.2f us,\t diffeo: %.2f us,\t spdiff: %.2f us\n'], ...
                 elapsed(1)*1e6, elapsed(2)*1e6, elapsed(3)*1e6, elapsed(4)*1e6);
    end
end

function [fig_elapsed, fig_error] = logm_fig(dims, elapsed, errors, opts)
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'savefig'), opts.savefig = true; end

    labels_time = {'built-in logm', ...
                   'diffeomorphic logm (unknown shift)', ...
                   'diffeomorphic logm (known shift)'};

    labels_error = {'built-in logm', 'Schur-based logm'};

    fig_elapsed = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1 1 5 5], 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);

    plot(dims, elapsed(:,1), '-o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
    plot(dims, elapsed(:,3), '-^', 'LineWidth', 1.2, 'MarkerSize', 4);
    plot(dims, elapsed(:,4), '-d', 'LineWidth', 1.2, 'MarkerSize', 4);
    grid on;
    xlabel('dimension n');
    ylabel('time (s)');
    title('Elapsed time for logarithm computations');
    legend(labels_time, 'Interpreter', 'none', 'Location', 'northwest');

    fig_error = figure('Color', 'w', 'Units', 'inches', ...
        'Position', [1 1 5 5], 'PaperUnits', 'inches', ...
        'PaperPosition', [0 0 5 5], 'PaperSize', [5 5]);

    semilogy(dims, errors(:,1), '-o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
    semilogy(dims, errors(:,2), '-s', 'LineWidth', 1.2, 'MarkerSize', 4);
    grid on;
    xlabel('dimension n');
    ylabel('reconstruction error');
    title('Reconstruction error: ||Q - expm(X)||_F');
    legend(labels_error, 'Interpreter', 'none', 'Location', 'northwest');

    T_elapsed = array2table([dims(:), elapsed], ...
        'VariableNames', {'n', ...
                          'builtin_logm', ...
                          'schur_based_logm', ...
                          'diffeomorphic_logm_unknown_D_component', ...
                          'diffeomorphic_logm_known_D_component'});

    T_error = array2table([dims(:), errors(:,1), errors(:,2)], ...
        'VariableNames', {'n', ...
                          'builtin_logm', ...
                          'schur_based_logm'});

    if opts.savefig
        exportgraphics(fig_elapsed, 'logm_elapsed.pdf', 'ContentType', 'vector');
        exportgraphics(fig_error, 'logm_error.pdf', 'ContentType', 'vector');
        writetable(T_elapsed, 'test_logm_elapsed.csv');
        writetable(T_error, 'test_logm_error.csv');
    end
end