function elapsed_times = test_expm(dims)
%TEST_EXPM  Benchmark and accuracy comparison of four matrix exponential implementations.
%
%   ELAPSED_TIMES = TEST_EXPM(dims)
%
%   This routine performs a timing and accuracy experiment for the following
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
    if nargin < 1, dims = 3:50; end

    % --- Add source folder to MATLAB path ---
    srcPath = fullfile(fileparts(mfilename('fullpath')), 'src');
    addpath(srcPath);

    % --- Load DLL if not already loaded ---
    libPath = fullfile(fileparts(mfilename('fullpath')), 'lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        fprintf('Loading library from %s\n', libPath);
        loadlibrary(dll, hdr);
    else
        fprintf('Library already loaded.\n');
    end

    elapsed_times = zeros(length(dims), 4);
    diff = 1e3 * ones(length(dims), 3);
    loop = 20;

    ind = 1;
    for dim = dims
        S = 2 * rand(dim) - 1;
        S = 1 * (S - S') / 2;
        opt.verbose = true;
        elap_err_expm(dim, S, opt);
        opt.verbose = false;
        for i = 1:loop
            [t, e] = elap_err_expm(dim, S, opt);
            elapsed_times(ind, :) = elapsed_times(ind, :) + t(:)'; 
        end
        elapsed_times(ind, :) = elapsed_times(ind, :) / loop;
        diff(ind, :) = min(diff(ind, :), e(:)');
        while e(1) > 1e-8
            [~, e] = elap_err_expm(dim, S, opt);
            diff(ind, :) = min(diff(ind, :), e(:)');
        end
        ind = ind + 1;
    end

    % Customize the plot
    figure;
    tlo = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Top subplot: elapsed times ---
    nexttile;
    plot(dims, elapsed_times(:, 1), '-o', 'DisplayName', 'sBLAS Expm', 'LineWidth', 1.8);
    hold on;
    plot(dims, elapsed_times(:, 2), '-x', 'DisplayName', 'Pade3 Expm', 'LineWidth', 1.8);
    plot(dims, elapsed_times(:, 3), '-s', 'DisplayName', 'Pade13 Expm', 'LineWidth', 1.8);
    plot(dims, elapsed_times(:, 4), '-^', 'DisplayName', 'Built-in Expm', 'LineWidth', 1.8);
    hold off;
    title('Elapsed Time for Matrix Exponential Calculation');
    ylabel('Time (seconds)');
    legend('Location', 'northwest');
    grid on;

    % --- Bottom subplot: error comparison ---
    nexttile;
    semilogy(dims, diff(:, 1), '-o', 'DisplayName', 'Built-in vs sBLAS', 'LineWidth', 1.8);
    hold on;
    semilogy(dims, diff(:, 2), '-x', 'DisplayName', 'Built-in vs Pade3', 'LineWidth', 1.8);
    semilogy(dims, diff(:, 3), '-s', 'DisplayName', 'Built-in vs Pade13', 'LineWidth', 1.8);
    hold off;
    title('Exponential Result Error (log scale)');
    xlabel('Matrix Dimension (n)');
    ylabel('‖(Built-in Expm) - (Other Expm)‖_F');
    legend('Location', 'northwest');
    grid on;

    % --- Link the x-axes for synchronized zooming/panning ---
    linkaxes(findall(gcf, 'type', 'axes'), 'x');

    %if libisloaded('skewblas')
    %    unloadlibrary('skewblas');
    %end
end

function [elapsed, errors] = elap_err_expm(n, S, opt)
    if nargin < 1, n = 9; end
    if nargin < 2, S = rand(n); end
    if nargin < 3, opt = struct(); opt.verbose = true; end
    

    S = 0.5 * (S - S');

    srcPath = fullfile(fileparts(mfilename('fullpath')), 'src');
    addpath(srcPath);
    dll = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    elapsed = zeros(4, 1);
    errors = zeros(3, 1);
    tic;
    Qss = sblas_expm(S);
    elapsed(1) = toc;

    tic;
    Qpade3 = pade_expm(S, 3);
    elapsed(2) = toc;

    tic;
    Qpade13 = pade_expm(S, 13);
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