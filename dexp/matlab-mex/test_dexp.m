function [elapsed_forward, elapsed_inverse, fig_fwd, fig_inv] = test_dexp(dims, scale, opts)
%TEST_DEXP  Timing experiment for forward and inverse differentials of expm.
%
%   [EFWD, EINV] = TEST_DEXP(dims, scale)
%
%   This routine benchmarks four implementations of the matrix exponential
%   differential Dexp_S(X), both in **forward** mode:
%
%        Y = Dexp_S(X)
%
%   and in **inverse** mode:
%
%        X = (Dexp_S)^{-1}(Y)
%
%   The four compared methods are:
%
%       (1) sBLAS Dexp      – DLL backend (skewblas.dll)
%       (2) Pade_dexp p=3   – forward; p=1 for inverse
%       (3) Pade_dexp p=13  – forward; p=7 for inverse
%       (4) Daletskii–Kreĭn spectral formula (dk_dexp)
%
%   For each dimension n in dims:
%       • A random skew-symmetric S and direction X are generated.
%       • One warm-up call is executed for forward and inverse paths.
%       • Each method is timed over `loop = 50` repetitions.
%       • The timings are averaged over repetitions.
%
%   INPUT
%     dims    Vector of matrix sizes to test. Default: 3:50.
%     scale   Scalar used to scale S by multiplication. Default: 1.
%
%   OUTPUT
%     elapsed_forward   (#dims × 4) matrix of mean execution times (seconds)
%                       columns correspond to:
%                           1. sBLAS Dexp (fwd)
%                           2. Pade p=3 (fwd)
%                           3. Pade p=13 (fwd)
%                           4. Daletskii–Kreĭn (fwd)
%
%     elapsed_inverse   (#dims × 4) matrix of mean execution times (seconds)
%                       columns correspond to:
%                           1. sBLAS Dexp (inv)
%                           2. Pade p=1 (inv)
%                           3. Pade p=7 (inv)
%                           4. Daletskii–Kreĭn (inv)
%
%   PRODUCED FIGURES
%     (1) Forward-mode Dexp timing vs. matrix dimension
%     (2) Inverse-mode Dexp timing vs. matrix dimension
%     x-axes of both figures are linked for synchronized zooming.
%
%   NOTES
%     • S = (S - S')/2 * scale and X = (X - X')/2 are freshly generated for
%       each dimension n and reused during averaging for fair comparison.
%
%     • No accuracy test is performed here; only elapsed time is measured.
%
%     • The DLL 'skewblas.dll' is loaded automatically if needed.
%
%   See also: sblas_dexp, pade_dexp, dk_dexp.

    if nargin < 1, dims = 3:50; end
    if nargin < 2, scale = 1; end
    if nargin < 3, opts = struct('savefig', true, 'verbose', true); end 
    
    if ~isfield(opts, 'savefig'), opts.savefig = true; end
    if ~isfield(opts, 'savefig'), opts.verbose = true; end

    % --- add src path ---

    addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'mex'));
    compile_mex();

    % % --- load DLL (for legacy wrappers that still call loadlibrary) ---
    % libPath = fullfile(fileparts(mfilename('fullpath')), 'lib');
    % dll = fullfile(libPath, 'skewblas.dll');
    % hdr = fullfile(libPath, 'skewblas_api.h');
    % if ~libisloaded('skewblas')
    %     loadlibrary(dll, hdr);
    % end

    loop = 30;

    elapsed_forward = zeros(length(dims), 4);
    elapsed_inverse = zeros(length(dims), 4);

    ind = 1;
    for dim = dims
        fprintf("Timing Dexp for n = %d\n", dim);

        S = 2 * rand(dim) - 1; S = 0.5*(S - S') * scale;
        X = 2 * rand(dim) - 1; X = 0.5*(X - X');

        % warm-up
        elap_forward(dim, S, X, struct('verbose', true));
        elap_inverse(dim, S, X, struct('verbose', true));

        % repetitions
        for k = 1:loop
            t_f = elap_forward(dim, S, X, struct('verbose', false));
            t_i = elap_inverse(dim, S, X, struct('verbose', false));

            elapsed_forward(ind, :) = elapsed_forward(ind,:) + t_f(:)';
            elapsed_inverse(ind, :) = elapsed_inverse(ind,:) + t_i(:)';
        end

        elapsed_forward(ind, :) = elapsed_forward(ind,:) / loop;
        elapsed_inverse(ind, :) = elapsed_inverse(ind,:) / loop;

        ind = ind + 1;
    end

    [fig_fwd, fig_inv] = dexp_fig(dims, elapsed_forward, elapsed_inverse, opts);
end
function elapsed = elap_forward(n, S, X, opt)
    if nargin < 4, opt.verbose = true; end
    elapsed = zeros(4,1);

    % 1. sBLAS forward
    tic;
    sblas_dexp(S, X, 'fwd');
    elapsed(1) = toc;

    % 2. Pade p=3
    tic;
    pade_dexp(S, X, 'fwd', 3);
    elapsed(2) = toc;

    % 3. Pade p=13
    tic;
    pade_dexp(S, X, 'fwd', 13);
    elapsed(3) = toc;

    % 4. Daletskii–Krein spectral
    tic;
    dk_dexp(S, X, 'fwd');
    elapsed(4) = toc;

    if opt.verbose
        fprintf("Forward: sBLAS %.1f μs, p3 %.1f μs, p13 %.1f μs, DK %.1f μs\n", ...
            elapsed(1)*1e6, elapsed(2)*1e6, elapsed(3)*1e6, elapsed(4)*1e6);
    end
end

function elapsed = elap_inverse(n, S, X, opt)
    if nargin < 4, opt.verbose = true; end
    elapsed = zeros(4,1);

    % 1. sBLAS inverse
    tic;
    sblas_dexp(S, X, 'inv');
    elapsed(1) = toc;

    % 2. Padé inverse p=1
    tic;
    pade_dexp(S, X, 'inv', 1);
    elapsed(2) = toc;

    % 3. Padé inverse p=7
    tic;
    pade_dexp(S, X, 'inv', 7);
    elapsed(3) = toc;

    % 4. DK inverse
    tic;
    dk_dexp(S, X, 'inv');
    elapsed(4) = toc;

    if opt.verbose
        fprintf("Inverse: sBLAS %.1f μs, p1 %.1f μs, p7 %.1f μs, DK %.1f μs\n", ...
            elapsed(1)*1e6, elapsed(2)*1e6, elapsed(3)*1e6, elapsed(4)*1e6);
    end
end
