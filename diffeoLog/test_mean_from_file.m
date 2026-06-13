function [M_plog, M_dlog, stats, Prob] = test_mean_from_file(path, opts)
% Inputs:
%   path : folder containing saved CSV data.
%
%   opts.MaxIter
%   opts.AbsTol
%   opts.Step
%   opts.InitCnt    optional, default 3
%
% Outputs:
%   M_plog : T-by-n-by-n array of principal-log means.
%   M_dlog : T-by-n-by-n array of diffeomorphic-log means.
%   stats  : table with iteration counts, residuals, and elapsed times.
%   Prob   : loaded T-by-D-by-n-by-n problem array.

    addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'mex'));
    addpath(fullfile(fileparts(mfilename('fullpath')), 'test'));
    compile_mex();

    if nargin < 1 || isempty(path)
        path = "./test/figure";
    end


    if nargin < 2
        opts = struct();
    end

    if ~isfield(opts, "MaxIter"), opts.MaxIter = 1000; end
    if ~isfield(opts, "AbsTol"),  opts.AbsTol  = 1e-9;  end

    [Prob, info] = read_moving_mean_problems(path);

    T = info.T;
    D = info.D;
    n = info.n;

    M_plog = zeros(T, n, n);
    M_dlog = zeros(T, n, n);

    iters_plog = zeros(T, 1);
    objvals_plog = zeros(T, 1);
    elapsed_plog_ns = zeros(T, 1);

    iters_dlog = zeros(T, 1);
    objvals_dlog = zeros(T, 1);
    elapsed_dlog_ns = zeros(T, 1);

    for t = 1:T

        VecMatQ = cell(D, 1);
        VecDelta0 = cell(D, 1);

        Omega = zeros(n, n);

        for d = 1:D
            Dtd = squeeze(Prob(t,d,:,:));
            Dtd = skew_part(Dtd);

            VecDelta0{d} = Dtd;
            VecMatQ{d} = expm(Dtd);

            Omega = Omega + Dtd;
        end

        Omega = skew_part(Omega / D);

        % Bad initial condition:
        %   E_bad = exp(-mean_i D_i).
        MatE_bad = expm(-Omega);

        % ------------------------------------------------------------
        % Principal-log initial Delta_i at E_bad
        % ------------------------------------------------------------

        VecDelta_plog = cell(D, 1);

        for d = 1:D
            Yd = VecMatQ{d}' * MatE_bad;
            VecDelta_plog{d} = skew_part(real(logm(Yd)));
        end

        % ------------------------------------------------------------
        % Diffeomorphic-log initial Delta_i at E_bad
        % ------------------------------------------------------------

        VecDelta_dlog = init_dlog_delta(Omega, VecDelta0, 1);

        % ------------------------------------------------------------
        % Run the two solvers from the same bad initial mean
        % ------------------------------------------------------------

        [MatM_plog, iter_p, objval_p, elapsed_p] = ...
            so_mean_plog(VecMatQ, VecDelta_plog, MatE_bad, opts);

        [MatM_dlog, iter_d, objval_d, elapsed_d] = ...
            so_mean_dlog(VecMatQ, VecDelta_dlog, MatE_bad, opts);

        M_plog(t,:,:) = MatM_plog;
        M_dlog(t,:,:) = MatM_dlog;

        iters_plog(t) = iter_p;
        objvals_plog(t) = objval_p;
        elapsed_plog_ns(t) = elapsed_p;

        iters_dlog(t) = iter_d;
        objvals_dlog(t) = objval_d;
        elapsed_dlog_ns(t) = elapsed_d;

        if mod(t, 100) == 0 || t == 1 || t == T
            fprintf("Problem %d / %d | plog: iter = %d, res = %.3e | dlog: iter = %d, res = %.3e\n", ...
                t, T, iter_p, objval_p, iter_d, objval_d);
        end
    end

    stats = table( ...
        (1:T)', ...
        iters_plog, objvals_plog, elapsed_plog_ns, ...
        iters_dlog, objvals_dlog, elapsed_dlog_ns, ...
        'VariableNames', { ...
            'problem_index', ...
            'plog_iterations', 'plog_residual', 'plog_elapsed_ns', ...
            'dlog_iterations', 'dlog_residual', 'dlog_elapsed_ns' ...
        });
    plot_mean(M_plog, M_dlog, Prob);
end


function S = skew_part(A)
    S = 0.5 * (A - A');
end

function VecDelta = init_dlog_delta(Omega, VecDelta0, cnt)

    if nargin < 3 || isempty(cnt)
        cnt = 3;
    end

    D = numel(VecDelta0);

    Omega = skew_part(Omega);

    VecMatE = cell(cnt, 1);
    for j = 1:cnt
        VecMatE{j} = expm((-j / cnt) * Omega);
    end

    VecDelta = cell(D, 1);

    for d = 1:D
        Dd = skew_part(VecDelta0{d});
        Qd = expm(Dd);
        nlogX = -Dd;

        for j = 1:cnt
            Y = Qd' * VecMatE{j};
            nlogX = mex_sblas_spdiff_logm(Y, nlogX);
            nlogX = skew_part(nlogX);
        end

        VecDelta{d} = nlogX;
    end
end