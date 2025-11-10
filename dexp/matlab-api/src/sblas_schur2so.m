function Q = sblas_schur2so(R, A)
% BuildSpecOrthSchur  Reconstruct the orthogonal matrix Q = R*exp(A)*R'.
%
%   Q = BuildSpecOrthSchur(R, A)
%
%   Inputs:
%     R : n-by-n real orthogonal matrix
%     A : vector of size n/2 containing principal angles
%
%   Output:
%     Q : n-by-n special orthogonal matrix (rotation)

    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    n = size(R, 1);

    Qptr = libpointer('doublePtr', zeros(n));

    calllib('skewblas', 'BuildSpecOrthSchur', Qptr, R, A, n);

    Q = Qptr.Value;
end
