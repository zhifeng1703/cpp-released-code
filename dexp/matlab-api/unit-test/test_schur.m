function [expm_time, expss_time] = test_schur(n, S, opt)
    if nargin < 1, n = 9; end
    if nargin < 2, S = rand(n); end
    if nargin < 3, opt = struct(); opt.verbose = false; end
    

    S = 0.5 * (S - S.');

    srcPath = fullfile(fileparts(mfilename('fullpath')), 'src');
    addpath(srcPath);
    dll = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    
    tic;
    Q = expm(S);
    Qrec = Q; QreclogQ = Q;
    expm_time = toc;

    % ---- Test SkewSymmSchur ----
    tic;
    [R, A] = sblas_ss2schur(S);
    Qrec = sblas_schur2so(R, A);
    expss_time = toc;

    Srec = sblas_schur2ss(R, A);
    [R, A] = sblas_so2schur(Q);
    QreclogQ = sblas_schur2so(R, A);

    if opt.verbose
        fprintf('Testing Schur decompodition with size n = %d\n', n);
        fprintf('‖S - RAR^T‖_F = %.3e, ‖exp(S) - Rexp(A)R^T‖_F = %.3e, ‖Q - Rexp(A)R^T‖_F = %.3e\n', ...
        norm(S - Srec, 'fro'), norm(Q - Qrec, 'fro'), norm(Q - QreclogQ, 'fro'));
        fprintf('Elapsed time: built-in expm: %.2f mili-seconds,\t expschur: %.2f mili-seconds\n', expm_time*1e6, expss_time*1e6);
    end
    


end