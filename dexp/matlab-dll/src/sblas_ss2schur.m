function [R, A] = sblas_ss2schur(S)
% SkewSymmSchur  MATLAB wrapper for the C++ DLL routine SkewSymmSchur.
%
%   [R, A] = SkewSymmSchur(S)
%
%   Inputs:
%     S : n×n real skew-symmetric matrix.
%
%   Outputs:
%     R : n×n orthogonal matrix.
%     A : n/2 × 1 vector of angles.


    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    n = size(S, 1);

    TmpR = libpointer('doublePtr', zeros(n));
    TmpA = libpointer('doublePtr', zeros(1, floor(n/2)));

    calllib('skewblas', 'SkewSymmSchur', TmpR, TmpA, S, n);

    % MATLAB passes by reference, so fetch updated values:
    R = TmpR.Value;
    A = TmpA.Value;
end
