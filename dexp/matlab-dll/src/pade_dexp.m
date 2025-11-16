function Y = pade_dexp(S, X, dir, p)
% pade_dexp  Differential of exp(S) using Pade approximations.
%
%   Y = pade_dexp(S, X, dir, p)
%
%   Inputs:
%     S   : n×n real SKEW-SYMMETRIC matrix
%     X   : n×n real skew-symmetric matrix (direction)
%     dir : 'fwd' or 'inv'
%     p   : Pade order (default 13 for fwd, 7 for inv)
%
%   Output:
%     Y : differential result
%
%   Note:
%     - FOR INVERSE: Schur decomposition (R,A) is required.
%     - This version does NOT store intermediate results.

    % ---------------------------------------------------------------
    % Load library if needed
    % ---------------------------------------------------------------
    libPath = fullfile(fileparts(mfilename('fullpath')), '..', 'lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    n = size(S, 1);

    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end

    if nargin < 4 || isempty(p)
        if startsWith(lower(dir), 'i')
            p = 7;   % good inverse default
        else
            p = 13;  % good forward default
        end
    end

    % ---------------------------------------------------------------
    % Prepare input/output pointers
    % ---------------------------------------------------------------
    PtrY = libpointer('doublePtr', zeros(n));

    if startsWith(lower(dir), 'f')      % ------------------ FORWARD
        calllib('skewblas', 'DexpPadeForward', ...
            PtrY, X, S, n, p);

    else                                % ------------------ INVERSE
        % 1. Compute Schur decomposition of S
        R  = zeros(n); 
        A  = zeros(floor(n/2), 1);

        PtrR = libpointer('doublePtr', R(:));
        PtrA = libpointer('doublePtr', A(:));

        calllib('skewblas', 'SkewSymmSchur', PtrR, PtrA, S, n);

        % 2. Compute matrix 1-norm
        normS = norm(S, 1);

        % 3. Call inverse differential with new API
        calllib('skewblas', 'DexpPadeInverse', ...
            PtrY, X, PtrR, PtrA, normS, n, p);
    end

    % ---------------------------------------------------------------
    % Retrieve result
    % ---------------------------------------------------------------
    Y = PtrY.Value;
end
