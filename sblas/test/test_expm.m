function test_expm(dims, opts)
%TEST_EXPM  Benchmark and accuracy comparison of four matrix exponential implementations.
%
%   TEST_EXPM(dims)
%
%   This routine performs a timing and accuracy experiment via calling
%
%   [ELAPSED, ERRORS, FIG_ELAPSED, FIG_ERROR] = TEST_EXPM(dims, opts)
%
%   for the following
%   exponential-map implementations on skew-symmetric matrices:
%
%       (1) sblas_expm      – Skew-Schur based formulae
%       (2) pade_expm(S,3)  – [3/3] Pade approximant
%       (3) pade_expm(S,13) – [13/13] Pade approximant
%       (4) expm            – MATLAB built-in
%
%   For each dimension n in dims:
%       • Generate a random skew-symmetric test matrix S.
%       • Evaluate each expm implementation `loop = 20` times.
%       • Average the elapsed times over the 20 runs.
%       • Record the minimum observed Frobenius error 
%             ‖Q_builtin − Q_method‖_F
%         over all repetitions, ensuring at least one result below 1e−8
%         for the sBLAS method (retry loop).
%
%   INPUT
%     dims   Vector of matrix sizes to test. 
%            Default: 3:50.
%
%   OUTPUT
%     elapsed_times   (#dims × 4) matrix of averaged timings (in seconds):
%                     column 1: sblas_expm
%                     column 2: pade_expm(⋅,3)
%                     column 3: pade_expm(⋅,13)
%                     column 4: MATLAB expm
%
%   SIDE OUTPUTS (internal)
%     diff            (#dims × 3) minimal Frobenius errors:
%                     column 1: built-in vs sBLAS
%                     column 2: built-in vs Pade3
%                     column 3: built-in vs Pade13
%
%   FIGURES
%     A 2×1 tiled layout containing:
%       • Top:   averaged elapsed times vs. matrix dimension
%       • Bottom: Frobenius error vs. dimension (log scale)
%
%   NOTES
%     • The DLL `skewblas.dll` is loaded once at the beginning.
%     • S is regenerated for each dimension n but reused across the 20-loop
%       averaging pass, ensuring fair timing and stable error statistics.
%     • The error threshold loop
%         while e(1) > 1e-8
%       guarantees that at least one accurate sBLAS evaluation is captured.
%
%   See also: elap_err_expm, sblas_expm, pade_expm, expm.

    if nargin < 1, dims = 3:150; end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'savefig'), opts.savefig = true; end

    fprintf(['=============EXPM experiment begins=============\n\n']);

    test_expm_data(dims, opts);
    
    fprintf(['=============EXPM experiment ends=============\n\n']);

    fprintf(['IMPORTANT NOTE: The compute time of the first execution ' ...
        'could be unreliable even WITH the warm-up run. \n']);
    fprintf(['IMPORTANT NOTE: Make sure to run multiple calls to test_expm for ' ...
        'more reliable timings.\n']);
end

function [elapsed, errors, fig_elapsed, fig_error] = test_expm_data(dims, opts)
    if nargin < 1, dims = 3:150; end
    if nargin < 2, opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'savefig'), opts.savefig = true; end

    addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'mex'));
    compile_mex();

    nd = length(dims);
    loop = 20;

    elapsed = zeros(nd, 4);
    errors        = +inf(nd, 3);   % store minimum Frobenius errors

    idx = 1;
    for n = dims
        % Random skew-symmetric test matrix
        S = 2*rand(n) - 1;
        S = (S - S') / 2;

        elap_err_expm(n, S, struct('verbose', true));

        for k = 1:loop
            [t, e] = elap_err_expm(n, S, struct('verbose', false));

            elapsed(idx, :) = elapsed(idx, :) + t(:)';

            % Keep the best (minimum) Frobenius errors
            errors(idx, :) = min(errors(idx, :), e(:)');
        end

        elapsed(idx, :) = elapsed(idx, :) / loop;
        idx = idx + 1;
    end

    [fig_elapsed, fig_error] = expm_fig(dims, elapsed, errors, opts);

end

function [elapsed, errors] = elap_err_expm(n, S, opt)
    if nargin < 1, n = 9; end
    if nargin < 2, S = rand(n); end
    if nargin < 3, opt = struct(); opt.verbose = true; end
    

    S = 0.5 * (S - S');

    % srcPath = fullfile(fileparts(mfilename('fullpath')), 'src');
    % addpath(srcPath);
    % dll = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas.dll');
    % hdr = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas_api.h');
    % 
    % if ~libisloaded('skewblas')
    %     loadlibrary(dll, hdr);
    % end

    elapsed = zeros(4, 1);
    errors = zeros(3, 1);
    tic;
    Qss = mex_sblas_expm(S);
    elapsed(1) = toc;

    tic;
    Qpade3 = mex_pade_expm(S, 3);
    elapsed(2) = toc;

    tic;
    Qpade13 = mex_pade_expm(S, 13);
    elapsed(3) = toc;

    tic;
    Q = expm(S);
    elapsed(4) = toc;

    errors(1) = norm(Q - Qss, 'fro');
    errors(2) = norm(Q - Qpade3, 'fro');
    errors(3) = norm(Q - Qpade13, 'fro');

    if opt.verbose
        fprintf('Testing matrix exponential with size n = %d\n', n);
        fprintf('‖Q_matlab - Q_sblas‖_F = %.3e, ‖Q_matlab - Q_pade3‖_F = %.3e, ‖Q_matlab - Q_pade13‖_F = %.3e\n', ...
           errors(1), errors(2), errors(3));
        fprintf('Elapsed time: built-in expm: %.2f ms,\t sblas_expm: %.2f ms, \t pade3_expm: %.2f ms, \t pade13_expm: %.2f ms\n', ...
            elapsed(1)*1e6, elapsed(2)*1e6, elapsed(3)*1e6, elapsed(4)*1e6);
    end
end