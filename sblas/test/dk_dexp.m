function [Y, para] = dk_dexp(S, X, dir, para)
% DK_DEXP  Daletskii–Kreĭn differential of expm(S) (skew-symmetric case).
%
%   [Y, para] = dk_dexp(S, X, dir)
%   [Y, para] = dk_dexp(S, X, dir, para)
%
%   Computes either:
%       Y = Dexp_S(X)          (dir = 'fwd')
%       Y = Dexp_S^{-1}(X)     (dir = 'inv')
%
%   Uses:
%       • Real skew-Schur decomposition of S
%       • Complex eigenvectors/eigenvalues
%       • DK forward & inverse parameters from C++ backend
%
%   para (optional) may contain:
%       para.R              real Schur vectors
%       para.A              Schur block parameters
%       para.EigVec         complex n×n eigenvectors
%       para.EigVal         complex n×1 eigenvalues
%       para.PF             DK forward parameters (interleaved real)
%       para.PI             DK inverse parameters (interleaved real)
%
%   All fields are computed if missing.
%

    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end
    if nargin < 4
        para = struct();
    end

    % ===============================================================
    % Ensure skew-symmetric inputs
    % ===============================================================
    S = 0.5*(S - S.');
    X = 0.5*(X - X.');
    n = size(S,1);

    % ===============================================================
    % 1. Real Schur decomposition if missing
    % ===============================================================
    if ~isfield(para,'R') || ~isfield(para,'A')
        [R, A] = mex_sblas_ss2schur(S);
        para.R = R;
        para.A = A;
    end

    % ===============================================================
    % 2. Build eigenvectors/eigenvalues if missing
    % ===============================================================
    if ~isfield(para,'EigVec') || ~isfield(para,'EigVal')
        [EigVec, EigVal] = schur_to_eig(para.R, para.A);
        para.EigVec = EigVec;
        para.EigVal = EigVal;
    end

    % ===============================================================
    % 3. DK parameters (PF_int, PI_int) if missing
    % ===============================================================
    if ~isfield(para,'PF_int') || ~isfield(para,'PI_int')
        [PF, PI] = mex_dk_dexp_para(para.EigVal, 'skew');
        para.PF = PF;
        para.PI = PI;
    end

    % ===============================================================
    % 4. Apply DK differential (forward or inverse)
    % ===============================================================
    if startsWith(lower(dir),'f')
        Y = mex_dk_dexp(para.EigVec, para.PF, para.PI, X, 'fwd');
    else
        Y = mex_dk_dexp(para.EigVec, para.PF, para.PI, X, 'inv');
    end
end



function [EigVec, EigVal] = schur_to_eig(R, A)
% SCHUR_TO_EIG  Convert real Schur form of a skew-symmetric matrix
%               into complex eigenpairs, packed in C++-style layout.
%
%   [EigVec, EigVal] = schur_to_eig(R, A)
%
%   Inputs:
%       R : n×n real orthogonal matrix (from skew-symmetric real Schur)
%       A : m×1 real vector of Schur angles (one per 2×2 block)
%
%   Outputs (C++-style interleaved real/imag):
%       EigVec : 2*n*n × 1 real vector
%                [Re v_11; Im v_11; Re v_21; Im v_21; ...] in column-major
%       EigVal : 2*n   × 1 real vector
%                [Re λ_1; Im λ_1; Re λ_2; Im λ_2; ...]
%
%   For each 2×2 skew-symmetric block J(a_k) = [0 -a_k; a_k 0], we use
%       v_- = (e1 + i e2)/sqrt(2),   λ_- = -i a_k
%       v_+ = (e1 - i e2)/sqrt(2),   λ_+ =  i a_k
%   where e1 = R(:,2k-1), e2 = R(:,2k).

    n = size(R, 1);
    m = numel(A);

    Vc = zeros(n, n);      % complex eigenvectors
    Dc = zeros(n, 1);      % complex eigenvalues

    col = 1;
    for k = 1:m
        ak = A(k);

        e1 = R(:, 2*k - 1);
        e2 = R(:, 2*k);

        % Orthonormal eigenvectors for the 2×2 block
        v_minus = (e1 + 1i*e2) / sqrt(2);   % λ = -i*ak
        v_plus  = (e1 - 1i*e2) / sqrt(2);   % λ =  i*ak

        Vc(:, col)     = v_minus;
        Dc(col)        = -1i * ak;

        Vc(:, col + 1) = v_plus;
        Dc(col + 1)    =  1i * ak;

        col = col + 2;
    end

    % Handle odd dimension: one zero eigenvalue with eigenvector R(:,end)
    if 2*m < n
        Vc(:, end) = R(:, end);
        Dc(end)    = 0;
    end

    % Pack into interleaved real/imag vectors (C++ CMPX_TYPE layout)
    EigVec = zeros(2*n*n, 1);
    EigVal = zeros(2*n,   1);

    vlin = Vc(:);   % n*n complex, column-major
    dlin = Dc(:);   % n×1 complex

    EigVec(1:2:end) = real(vlin);
    EigVec(2:2:end) = imag(vlin);

    EigVal(1:2:end) = real(dlin);
    EigVal(2:2:end) = imag(dlin);
end
