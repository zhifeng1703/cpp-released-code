function [Y, para] = sblas_dexp(S, X, dir, para)
% sblas_dexp  MATLAB wrapper for the differential of exp(S)
%              on skew-symmetric matrices.
%
%   [Y, para] = sblas_dexp(S, X, dir, para)
%
%   Computes the forward or inverse differential of the matrix exponential
%   for a real skew-symmetric generator S, using the C++ sBLAS library.
%
%   Inputs:
%     S    : n×n real skew-symmetric matrix (generator)
%     X    : n×n real skew-symmetric perturbation (direction)
%     dir  : (optional) 'fwd' (default) or 'inv'
%     para : (optional) structure containing precomputed fields:
%              .R            - Schur orthogonal matrix
%              .A            - Schur angles
%              .ParaForward  - forward differential parameters
%              .ParaInverse  - inverse differential parameters
%
%   Outputs:
%     Y    : n×n real skew-symmetric matrix result
%     para : structure (newly created or reused)
%
%   Requires:
%     skewblas.dll and skewblas_api.h in ../lib

    % ---------------------------------------------------------------------
    % Load library if needed
    % ---------------------------------------------------------------------
    dll = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas.dll');
    hdr = fullfile(fileparts(mfilename('fullpath')), '..', 'lib', 'skewblas_api.h');
    if ~libisloaded('skewblas')
        loadlibrary(dll, hdr);
    end

    % ---------------------------------------------------------------------
    % Default arguments
    % ---------------------------------------------------------------------
    if nargin < 3 || isempty(dir)
        dir = 'fwd';
    end
    n = size(S, 1);
    a = floor(n/2);
    parvec_size = a*(a-1)*8 + a*4;

    % ---------------------------------------------------------------------
    % Compute Schur decomposition and differential parameters if needed
    % ---------------------------------------------------------------------
    if nargin < 4 || isempty(para) || ~isfield(para, 'R') || ~isfield(para, 'A')
        % --- Step 1: Schur decomposition ---
        TmpR = libpointer('doublePtr', zeros(n));
        TmpA = libpointer('doublePtr', zeros(1, a));
        calllib('skewblas', 'SkewSymmSchur', TmpR, TmpA, S, n);

        para.R = TmpR.Value;
        para.A = TmpA.Value;

        % --- Step 2: Differential parameters ---
        ParaForward = libpointer('doublePtr', zeros(parvec_size, 1));
        ParaInverse = libpointer('doublePtr', zeros(parvec_size, 1));
        calllib('skewblas', 'DexpSkewSymmPara', ...
                ParaForward, ParaInverse, para.A, n);

        para.ParaForward = ParaForward.Value;
        para.ParaInverse = ParaInverse.Value;
    else
        % Ensure the fields exist with correct size
        if ~isfield(para, 'ParaForward') || numel(para.ParaForward) ~= parvec_size
            para.ParaForward = zeros(parvec_size, 1);
        end
        if ~isfield(para, 'ParaInverse') || numel(para.ParaInverse) ~= parvec_size
            para.ParaInverse = zeros(parvec_size, 1);
        end
    end

    % ---------------------------------------------------------------------
    % Prepare matrices for differential call
    % ---------------------------------------------------------------------
    TmpY = libpointer('doublePtr', zeros(n));

    % ---------------------------------------------------------------------
    % Call the corresponding DLL routine
    % ---------------------------------------------------------------------
    if startsWith(lower(dir), 'i')
        calllib('skewblas', 'DexpSkewSymmInverse', ...
            TmpY, X, para.R, para.ParaInverse, n);
    else
        calllib('skewblas', 'DexpSkewSymmForward', ...
            TmpY, X, para.R, para.ParaForward, n);
    end

    % ---------------------------------------------------------------------
    % Retrieve output and updated parameter arrays
    % ---------------------------------------------------------------------
    Y = reshape(TmpY.Value, n, n);
end
