function S = sblas_Schur2ss(R, A)
% BuildSkewSymmSchur  Reconstruct the skew-symmetric matrix S = R*A*R'.
%
%   S = BuildSkewSymmSchur(R, A)
%
%   Inputs:
%     R : n-by-n real orthogonal matrix
%     A : vector of size n/2 containing angular parameters
%
%   Output:
%     S : n-by-n real skew-symmetric matrix

    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    n = size(R, 1);

    Sptr = libpointer('doublePtr', zeros(n));
    Rptr = libpointer('doublePtr', R);
    Aptr = libpointer('doublePtr', A);

    calllib('skewblas', 'BuildSkewSymmSchur', Sptr, Rptr, Aptr, n);

    S = Sptr.Value;
end
