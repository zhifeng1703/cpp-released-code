function Q = sblas_expm(S)
% sblas_expm  MATLAB wrapper for the C++ DLL routine ExpmSkewSymm.
%
%   [Q, ssf] = sblas_expm(S)
%
%   Computes the exponential of a real skew-symmetric matrix S
%   using the Schur-based algorithm implemented in the C++ library.
%
%   Inputs:
%     S : n×n real skew-symmetric matrix.
%
%   Outputs:
%     Q   : n×n orthogonal matrix representing exp(S).
%     ssf : struct containing the Schur decomposition factors:
%           ssf.R  - n×n orthogonal matrix
%           ssf.A  - n/2 × 1 vector of angles
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
    n = size(S, 1);

    TmpQ = libpointer('doublePtr', zeros(n));          % exp(S)
    %TmpR = libpointer('doublePtr', zeros(n));          % Schur orthogonal
    %TmpA = libpointer('doublePtr', zeros(1, floor(n/2))); % angle vector

    % ---------------------------------------------------------------------
    % Call the C++ API function
    % ---------------------------------------------------------------------
    %calllib('skewblas', 'ExpmSkewSymm', TmpQ, TmpR, TmpA, S, n);
    calllib('skewblas', 'ExpmSkewSymm', TmpQ, S, n);

    % ---------------------------------------------------------------------
    % Retrieve updated values
    % ---------------------------------------------------------------------
    Q = TmpQ.Value;
    %ssf.R = TmpR.Value;
    %ssf.A = TmpA.Value;
end
