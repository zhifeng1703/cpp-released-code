function [R, A] = sblas_so2schur(Q)
% SpecOrthSchur  MATLAB wrapper for the C++ DLL routine SpecOrthSchur.
%
%   [R, A] = SpecOrthSchur(Q)
%
%   Inputs:
%     Q : n×n real special-orthogonal matrix.
%
%   Outputs:
%     R : n×n orthogonal matrix.
%     A : n/2 × 1 vector of angles.

    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    n = size(Q, 1);

    TmpR = libpointer('doublePtr', zeros(n));
    TmpA = libpointer('doublePtr', zeros(1, floor(n/2)));

    calllib('skewblas', 'SpecOrthSchur', TmpR, TmpA, Q, n);

    % MATLAB passes by reference, so fetch updated values:
    R = TmpR.Value;
    A = TmpA.Value;
end
