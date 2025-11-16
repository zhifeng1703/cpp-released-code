function [Y, para] = sblas_dexp(S, X, dir, para)
% sblas_dexp  Differential of exp(S) using the sBLAS MEX backend.
%
%   [Y, para] = sblas_dexp(S, X, dir, para)
%
%   Inputs:
%     S    : n×n skew-symmetric matrix
%     X    : n×n skew-symmetric direction
%     dir  : 'fwd' (default) or 'inv'
%     para : optional cached structure containing:
%              .R            – Schur orthogonal matrix
%              .A            – Schur angles
%              .ParaFwd      – forward diff parameters
%              .ParaInv      – inverse diff parameters
%
%   Output:
%     Y    : result of Dexp_S(X) or Dexp_S^{-1}(X)
%     para : cached fields for reuse

    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end

    n = size(S,1);
    a = floor(n/2);
    parvec_size = a*(a-1)*8 + a*4;   % C++ parameter size (real)

    % ---------------------------------------------------------------
    % Compute Schur decomposition and DK parameters if missing
    % ---------------------------------------------------------------
    need_schur = (nargin < 4) || isempty(para) || ...
                 ~isfield(para,'R') || ~isfield(para,'A');

    if need_schur
        % --- Schur ---
        [R, A] = mex_sblas_ss2schur(S);
        para.R = R;
        para.A = A;

        % --- Parameter vectors ---
        % C++ expects REAL arrays, not complex
        [PF, PI] = mex_sblas_dexp_para(A);

        para.Fwd = PF;
        para.Inv = PI;
    else
        % Validate parameter lengths
        if ~isfield(para,'Fwd') || numel(para.ParaFwd) ~= parvec_size
            error('Invalid or missing para.ParaForward.');
        end
        if ~isfield(para,'Inv') || numel(para.ParaInv) ~= parvec_size
            error('Invalid or missing para.ParaInverse.');
        end
    end

    % ---------------------------------------------------------------
    % Call the actual MEX differential
    % ---------------------------------------------------------------
    if startsWith(lower(dir),'i')
        Y = mex_sblas_dexp(S, X, 'inv', para);
    else
        Y = mex_sblas_dexp(S, X, 'fwd', para);
    end
end
