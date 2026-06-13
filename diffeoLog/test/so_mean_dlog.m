function [MatM, iter, objval, elapsed_ns] = so_mean_dlog(VecMatQ, VecDelta, MatE, opts)
%SO_MEAN_DLOG Diffeomorphic-log Karcher mean on SO(n).
% Inputs:
%   VecMatQ   : cell array of SO(n) matrices Q_i
%   VecDelta  : cell array of initial skew matrices Delta_i
%   opts.MaxIter
%   opts.AbsTol
%
% Outputs:
%   MatM       : final skew-symmetric logarithm of the mean
%   iter       : iteration count
%   objval     : norm of final gradient residual
%   elapsed_ns : elapsed wall-clock time in nanoseconds

    if nargin < 3
        opts = struct();
    end
    if ~isfield(opts, "MaxIter"), opts.MaxIter = 1000; end
    if ~isfield(opts, "AbsTol"),  opts.AbsTol  = 1e-6; end
    if ~isfield(opts, "Step"),  opts.Step  = 1; end

    D = numel(VecMatQ);
    n = size(VecMatQ{1}, 1);

    plogE = real(logm(MatE));
    MatA = skew_part(plogE);
    VecMatDelta = VecDelta;

    t_start = tic;

    MatX = zeros(n, n);
    for d = 1:D
        MatX = MatX + VecMatDelta{d};
    end
    MatX = -(1.0 / D) * MatX;

    objval = norm(MatX, "fro");

    iter = 1;
    % step = opts.Step;

    while iter <= opts.MaxIter && objval > opts.AbsTol

        % Julia code uses stepsize = 1 in so_mean_plog.
        % stepsize = norm(MatX, 2);
        % if stepsize > opts.Step
        %     stepsize = opts.Step / stepsize;
        % else 
        %     stepsize = 1;
        % end

        stepsize = opts.Step;
        MatE = MatE * expm(stepsize * MatX);

        % Diffeomorphic logarithm of the current mean.
        MatA = mex_sblas_spdiff_logm(MatE, MatA);

        % Update diffeomorphic-log tangents:
        for d = 1:D
            MateY = VecMatQ{d}' * MatE;
            VecMatDelta{d} = mex_sblas_spdiff_logm(MateY, VecMatDelta{d});
        end

        MatX(:,:) = 0.0;
        for d = 1:D
            MatX = MatX + VecMatDelta{d};
        end
        MatX = -(1.0 / D) * MatX;

        objval = norm(MatX, "fro");
        iter = iter + 1;
    end

    elapsed_ns = ceil(toc(t_start) * 1e9);

    if objval >= opts.AbsTol
        warning("SO(n) diffeomorphic-log mean did not converge. Final residual = %.3e.", objval);
    end

    MatM = MatA;
end

function S = skew_part(A)
    S = 0.5 * (A - A');
end
