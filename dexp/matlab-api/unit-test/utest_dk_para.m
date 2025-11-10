function utest_dk_para()
    % ===============================================================
    % Load DLL
    % ===============================================================
    libPath = fullfile(fileparts(mfilename('fullpath')), 'lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        fprintf('Loading library...\n');
        loadlibrary(dll, hdr);
    end

    % ===============================================================
    % Generate small skew-symmetric matrix
    % ===============================================================
    n = 8;
    S = randn(n); 
    S = 0.5*(S - S.');

    fprintf('Testing DK parameter computation on n = %d\n', n);

    % ===============================================================
    % Step 1. Real Schur decomposition via DLL
    % ===============================================================
    TmpR = libpointer('doublePtr', zeros(n));
    TmpA = libpointer('doublePtr', zeros(1, floor(n/2)));

    calllib('skewblas', 'SkewSymmSchur', TmpR, TmpA, S, n);

    R = TmpR.Value;
    A = TmpA.Value(:);

    fprintf('Schur decomposition completed.\n');
    disp('R (orthogonal):'); disp(R);
    disp('A (Schur angles):'); disp(A.');

    % ===============================================================
    % Step 2. Convert Schur → Complex Eigenpairs
    % ===============================================================
    [EigVec, EigVal] = schur_to_eig(R, A);

    fprintf('Complex eigen-decomposition completed.\n');
    disp('EigVal (first few):');
    disp(EigVal(1:6).');

    % pack eigenvalues into interleaved double array
    ev = EigVal(:);
    ev_int = zeros(2*n,1); 
    ev_int(1:2:end) = real(ev);
    ev_int(2:2:end) = imag(ev);

    PtrEigVal = libpointer('doublePtr', ev_int);

    % ===============================================================
    % Step 3. Allocate complex buffers for DK parameters
    % ===============================================================
    CmpxParaForward = libpointer('doublePtr', zeros(2*n*n,1));
    CmpxParaInverse = libpointer('doublePtr', zeros(2*n*n,1));

    % ===============================================================
    % Step 4. Call C++ to compute DK parameters
    % ===============================================================
    calllib('skewblas', 'DexpDalKreinParaSkewSymm', ...
        CmpxParaForward, CmpxParaInverse, PtrEigVal, n);

    fprintf('DK parameter generation completed.\n');

    % ===============================================================
    % Step 5. Unpack returned DK parameter matrices
    % ===============================================================
    pf = CmpxParaForward.Value;
    pi = CmpxParaInverse.Value;

    ParaF = complex(pf(1:2:end), pf(2:2:end));
    ParaI = complex(pi(1:2:end), pi(2:2:end));

    ParaF = reshape(ParaF, n, n);
    ParaI = reshape(ParaI, n, n);

    % ===============================================================
    % Show result
    % ===============================================================
    disp('Forward DK parameter matrix (ParaForward):');
    disp(ParaF);

    disp('Inverse DK parameter matrix (ParaInverse):');
    disp(ParaI);

    % ===============================================================
    % Quick sanity checks
    % ===============================================================
    fprintf('\nSanity checks:\n');
    fprintf('  ‖ParaForward‖_F = %.3e\n', norm(ParaF, 'fro'));
    fprintf('  ‖ParaInverse‖_F = %.3e\n', norm(ParaI, 'fro'));

    % Some minimal structure checks (not exact identities)
    % Usually ParaForward and ParaInverse should have block symmetry
    fprintf('  isfinite: %d\n', all(isfinite([ParaF(:); ParaI(:)])));

    % ===============================================================
    % Done
    % ===============================================================
    if libisloaded('skewblas')
        unloadlibrary('skewblas');
    end
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
