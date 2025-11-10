function Q = pade_expm(A, p)
% pade_expm  MATLAB wrapper for the C++ DLL routine ExpmPade.
%
%   Q = pade_expm(A, p)
%
%   Computes the matrix exponential exp(A) using a fixed-order
%   Padé approximation implemented in the C++ library.
%
%   Inputs:
%     A : n×n real matrix (double)
%     p : Padé order (integer, e.g. 7, 9, 13)
%
%   Outputs:
%     Q : n×n real matrix, approximating exp(A)
%
%   Requires:
%     skewblas.dll and skewblas_api.h in ../lib relative to this file.

    % ---------------------------------------------------------------------
    % Load the shared library (only once)
    % ---------------------------------------------------------------------
    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    % ---------------------------------------------------------------------
    % Prepare input/output arguments
    % ---------------------------------------------------------------------
    n = size(A, 1);
    if nargin < 2, p = 7; end  % default Padé order if not provided

    MatQ = libpointer('doublePtr', zeros(n));   % output: exp(A)

    % ---------------------------------------------------------------------
    % Call the C++ API function
    % ---------------------------------------------------------------------
    calllib('skewblas', 'ExpmPade', MatQ, A, n, p);

    % ---------------------------------------------------------------------
    % Retrieve output
    % ---------------------------------------------------------------------
    Q = MatQ.Value;
end
