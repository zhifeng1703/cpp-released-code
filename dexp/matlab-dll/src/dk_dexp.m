function [Y, para] = dk_dexp(S, X, dir, para)
% dk_dexp  Daletskii–Kreĭn differential of exp(S), skew-symmetric case
%
%   [Y, para] = dk_dexp(S, X, dir, para)
%
%   Implements Dexp_S(X) using the complex DK formulation.

    % ---------------------------------------------------------------
    % Load DLL if needed
    % ---------------------------------------------------------------
    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');
    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end
    
    n = size(S, 1);

    % ===============================================================
    % Build parameters if missing
    % ===============================================================
    if nargin < 4 || isempty(para) || ~isfield(para, 'EigVec')
        % -----------------------------------------------------------
        % 1. Compute real Schur decomposition  (required pipeline)
        % -----------------------------------------------------------
        TmpR = libpointer('doublePtr', zeros(n));
        TmpA = libpointer('doublePtr', zeros(1, floor(n/2)));
        calllib('skewblas', 'SkewSymmSchur', TmpR, TmpA, S, n);

        para.R = TmpR.Value;
        para.A = TmpA.Value;

        % -----------------------------------------------------------
        % 2. Compute full complex eigenvectors and eigenvalues of S
        % -----------------------------------------------------------
        [para.EigVec, para.EigVal] = schur_to_eig(para.R, para.A);

        % -----------------------------------------------------------
        % 3. Allocate interleaved complex arrays for parameter buffers
        % -----------------------------------------------------------
        CmpxParaForward = libpointer('doublePtr', zeros(2*n*n,1));
        CmpxParaInverse = libpointer('doublePtr', zeros(2*n*n,1));
        CmpxEigVal      = libpointer('doublePtr', zeros(2*n,1));

        % pack eigenvalues into interleaved buffer
        ev = para.EigVal(:);
        m = length(ev);
        ev_int = zeros(2*n,1); 
        ev_int(1:2:(2*m)) = real(ev);
        ev_int(2:2:(2*m)) = imag(ev);
        CmpxEigVal.Value = ev_int;

        % -----------------------------------------------------------
        % 4. Generate DK parameters
        % -----------------------------------------------------------
        calllib('skewblas', 'DexpDalKreinParaSkewSymm', ...
            CmpxParaForward, CmpxParaInverse, CmpxEigVal, n);

        % unpack returned parameters
        pf = CmpxParaForward.Value;
        pi = CmpxParaInverse.Value;

        para.ParaForward = complex(pf(1:2:end), pf(2:2:end));
        para.ParaInverse = complex(pi(1:2:end), pi(2:2:end));
    end
    
    % ===============================================================
    % Prepare call for DexpDalKreinForward/Inverse
    % ===============================================================
    MatY = libpointer('doublePtr', zeros(n*n,1));   % n*n output buffer
    MatX = libpointer('doublePtr', X(:));           % n*n input buffer (flattened)
    
    % pack EigenVec (n×n complex) into interleaved double vector
    V = para.EigVec(:);
    Vint = zeros(2*n*n,1);
    Vint(1:2:end) = real(V);
    Vint(2:2:end) = imag(V);
    PtrEigVec = libpointer('doublePtr', Vint);
    
    % pack forward / inverse parameters (n×n complex)
    pf = para.ParaForward(:);
    pi = para.ParaInverse(:);
    
    pfint = zeros(2*n*n,1);  pfint(1:2:end)=real(pf);  pfint(2:2:end)=imag(pf);
    piint = zeros(2*n*n,1);  piint(1:2:end)=real(pi);  piint(2:2:end)=imag(pi);
    
    PtrForward = libpointer('doublePtr', pfint);
    PtrInverse = libpointer('doublePtr', piint);
    
    % ===============================================================
    % Forward / Inverse differential call
    % ===============================================================
    if startsWith(lower(dir), 'i')  % inverse
        calllib('skewblas', 'DexpDalKreinInverse', ...
            MatY, MatX, PtrEigVec, PtrInverse, n);
    else                            % forward
        calllib('skewblas', 'DexpDalKreinForward', ...
            MatY, MatX, PtrEigVec, PtrForward, n);
    end
    
    % ===============================================================
    % Retrieve output
    % ===============================================================
    Y = reshape(MatY.Value, n, n);

end

function [EigVec, EigVal] = schur_to_eig(R, A)
% Convert skew-symmetric Schur form (R, A) into full complex eigenvectors
%
% Inputs:
%   R : n×n real orthogonal matrix
%   A : (n/2)×1 vector of Schur angles
%
% Outputs:
%   EigVec : n×n complex eigenvector matrix
%   EigVal : n×1 complex eigenvalue vector

    n = size(R,1);
    m = length(A);   % number of 2×2 blocks

    % Canonical block eigenvectors
    Vtilde = zeros(n, n);  % real part first
    Vtilde_i = zeros(n, n); % imaginary part

    col = 1;

    for k = 1:m
        i1 = 2*k-1;
        i2 = 2*k;

        % +i*a_k eigenvector in canonical basis
        Vtilde(i1, col) = 1/sqrt(2);
        Vtilde_i(i2, col) = 1/sqrt(2);    % imaginary = [0;1]/sqrt(2)
        col = col + 1;

        % -i*a_k eigenvector
        Vtilde(i1, col) = 1/sqrt(2);
        Vtilde_i(i2, col) = -1/sqrt(2);
        col = col + 1;
    end

    % Combine real/imag parts
    Vtilde_c = complex(Vtilde, Vtilde_i);

    % Rotate to original basis
    EigVec = R * Vtilde_c;

    % Eigenvalues ± i a_k
    EigVal = reshape([1i*A(:), -1i*A(:)] .', [], 1);
end
