function [elapsed, errors] = utest_expm(n, S, opt)
%UTEST_EXPM  Unit test for matrix exponential implementations.
%
%   [elapsed, errors] = UTEST_EXPM(n, S, opt)
%
%   This routine compares four ways to compute expm(S):
%       1) sblas_expm        — C++/MKL backend (skewblas.dll)
%       2) pade_expm(S,3)    — custom [3/3] Pade approximant
%       3) pade_expm(S,13)   — custom [13/13] Pade approximant
%       4) expm(S)           — MATLAB built-in
%
%   INPUTS
%     n     (optional) Integer. Matrix size. Default: n = 9.
%     S     (optional) Real matrix. If supplied, it will be skew-symmetrized
%           internally by S <- (S - S')/2. Default: random n-by-n matrix.
%     opt   (optional) Struct with fields:
%               verbose  – logical flag, display diagnostic output (default: true)
%
%   OUTPUTS
%     elapsed   4×1 vector of wall-clock timings (in seconds):
%               elapsed(1): sblas_expm
%               elapsed(2): pade_expm(S,3)
%               elapsed(3): pade_expm(S,13)
%               elapsed(4): MATLAB expm(S)
%
%     errors    3×1 vector of Frobenius-norm discrepancies:
%               errors(1) = ‖expm(S) − sblas_expm(S)‖_F
%               errors(2) = ‖expm(S) − pade_expm(S,3)‖_F
%               errors(3) = ‖expm(S) − pade_expm(S,13)‖_F
%
%   NOTES
%     - The routine auto-loads 'skewblas.dll' via loadlibrary if it is not yet loaded.
%     - Assumes the wrapper functions sblas_expm and pade_expm are located in ./src.
%     - Timing uses MATLAB's tic/toc wall-clock measurement.
%
%   EXAMPLE
%     [t, e] = utest_expm(20);
%
%   See also: expm, loadlibrary.
%
%   (skewblas MATLAB test wrapper)

    if nargin < 1, n = 9; end
    if nargin < 2, S = rand(n); end
    if nargin < 3, opt = struct(); opt.verbose = true; end
    

    S = 0.5 * (S - S');

    srcPath = fullfile(fileparts(mfilename('fullpath')), '../src');
    addpath(srcPath);
    dll = fullfile(fileparts(mfilename('fullpath')), 'lib', '../skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), 'lib', '../skewblas_api.h');

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