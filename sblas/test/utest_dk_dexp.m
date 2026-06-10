function utest_dk_dexp(n)
% UTEST_DK_DEXP  Validate Daletskii–Kreĭn differential of matrix exponential.
%
%   utest_dk_dexp(n)
%
%   Inputs:
%       n - matrix size (default 8)
%
%   This test:
%       1) Constructs a random skew-symmetric matrix S = R D R'
%          using Schur blocks [0, -a_i; a_i, 0].
%       2) Computes DK forward/inverse differentials via MEX interface.
%       3) Validates the forward differential using finite differences.
%       4) Plots the convergence of FD errors as h → 0.
%
%   Requires:
%       mex_sblas_dkpara.mex*
%       mex_dk_dexp.mex*
%
%   Example:
%       utest_dk_dexp(6);

    % ============================================================
    % 0. Add skewblas/mex to path & build MEX if needed
    % ============================================================
    thisFile = mfilename('fullpath');
    rootDir  = fileparts(fileparts(thisFile));      % .../skewblas
    mexDir   = fullfile(rootDir, 'mex');

    if ~isfolder(mexDir)
        error('Cannot find skewblas/mex folder at: %s', mexDir);
    end

    % Add path if not already added
    if ~contains(path, mexDir)
        addpath(mexDir);
        fprintf('[utest] Added path: %s\n', mexDir);
    end

    % Required mex files for this test
    mex_required = { ...
        'mex_sblas_ss2schur', ...
        'mex_dk_dexp_para', ...
        'mex_dk_dexp' ...
    };

    % Check if any required MEX is missing
    need_build = false;
    for k = 1:numel(mex_required)
        if exist(mex_required{k}, 'file') ~= 3   % 3 = MEX-file
            need_build = true;
            fprintf('[utest] Missing MEX: %s.mex* — rebuilding...\n', mex_required{k});
        end
    end

    % Build if needed
    if need_build
        buildScript = fullfile(mexDir, 'build_mex.m');
        if ~isfile(buildScript)
            error('build_mex.m not found in: %s', mexDir);
        end
        run(buildScript);
    end

    if nargin < 1, n = 8; end
    fprintf('==== Daletskii–Kreĭn Differential Test (n = %d) ====\n', n);

    % ============================================================
    % 1. Construct skew-symmetric test matrix S = R D R'
    % ============================================================
    S = randn(n);
    S = 0.5 * (S - S');
    [R, A] = mex_sblas_ss2schur(S);
    DA = zeros(n);
    for i = 1:floor(n/2)
        DA(2*i, 2*i-1) = A(i);
        DA(2*i-1, 2*i) = -A(i);
    end

    % ============================================================
    % 2. Convert to eigenpairs (interleaved complex form)
    % ============================================================
    [EigVec, EigVal] = schur_to_eig(R, A);
    % [CmpxEigVec, CmpxEigVal] = unpack_cpp_complex(EigVec, EigVal, n);


    % ============================================================
    % 3. Compute Daletskii–Kreĭn parameters
    % ============================================================
    [PF, PI] = mex_dk_dexp_para(EigVal, 'skew');

    % ============================================================
    % 4. Prepare random test direction
    % ============================================================
    X = randn(n);
    X = 0.5 * (X - X');

    % ============================================================
    % 5. Forward and inverse DK differential (via MEX)
    % ============================================================
    Yf = mex_dk_dexp(EigVec, PF, PI, X, 'fwd');
    Yi = mex_dk_dexp(EigVec, PF, PI, Yf, 'inv');

    fprintf('Relative inverse consistency error: %.3e\n', ...
        norm(Yi - X, 'fro') / max(norm(X,'fro'),1e-14));

    % ============================================================
    % 6. Finite-difference validation for Dexp_S(X)
    % ============================================================
    h = 1e-6;
    ExpS  = expm(S);
    ExpSh = expm(S + h*X);
    Y_fd  = (ExpSh - ExpS) / h;

    % adjoint correction (ExpS' * Y_fd)
    err_fd = norm(Yf - ExpS' * Y_fd, 'fro') / max(norm(Yf, 'fro'), 1e-14);

    % ------------------------------------------------------------
    % Print diagnostics
    % ------------------------------------------------------------
    disp('--- Daletskii–Kreĭn Differential Test ---');
    disp('X (input direction):'); disp(X);
    disp('Yf (forward differential using DK):'); disp(Yf);
    disp('Yf_fd (forward finite-diff):'); disp(ExpS' * Y_fd);
    disp('Yi (inverse differential, should recover X):'); disp(Yi);

    fprintf('\nFinite difference validation for DK Dexp_S(X):\n');
    fprintf('  step size h = %.1e\n', h);
    fprintf('  rel. Frobenius error = %.3e\n', err_fd);

    % ============================================================
    % 7. Convergence test as h → 0
    % ============================================================
    hs = logspace(-1, -7, 7);
    errs = zeros(size(hs));
    for k = 1:numel(hs)
        Y_fd = ExpS' * (expm(S + hs(k)*X) - ExpS) / hs(k);
        errs(k) = norm(Yf - Y_fd, 'fro') / max(norm(Yf, 'fro'), 1e-14);
    end

    figure('Name','DK Dexp finite-diff convergence','Color','w');
    loglog(hs, errs, '-o', 'DisplayName', 'Relative Error', 'LineWidth', 1.8);
    grid on; hold on;
    xlabel('h (step size)');
    ylabel('Relative Error: DK Dexp(S)[X] vs FD');
    title(sprintf('Finite Difference Convergence Test for DK Dexp(S)[X], n=%d', n));
    set(gca, 'XDir', 'reverse');
    loglog(hs, hs, '--k', 'DisplayName', 'O(h) slope');
    legend('show', 'Location','best');
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


function [EigVec, EigVal] = unpack_cpp_complex(EigVec_int, EigVal_int, n)
    % ------------------------------
    % Unpack eigenvalues
    % ------------------------------
    reV = EigVal_int(1:2:end);   % odd entries
    imV = EigVal_int(2:2:end);   % even entries
    EigVal = complex(reV, imV);

    % ------------------------------
    % Unpack eigenvectors
    % ------------------------------
    reM = EigVec_int(1:2:end);
    imM = EigVec_int(2:2:end);

    V = complex(reM, imM);
    EigVec = reshape(V, n, n);   % column-major reshape
end
