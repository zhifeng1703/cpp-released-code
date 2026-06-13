function [Prob, info] = read_moving_mean_problems(path)
%READ_MOVING_MEAN_PROBLEMS Read Karcher mean problems saved from Julia.
%
%   [Prob, info] = read_moving_mean_problems(path)
%
% Input:
%   path : folder containing
%          karcher_mean_problems.csv
%          karcher_mean_info.txt
%
% Output:
%   Prob : T-by-D-by-n-by-n array.
%          squeeze(Prob(t,d,:,:)) is the d-th skew-symmetric matrix
%          in the t-th problem.
%
%   info : struct with fields T, D, n.

    if nargin < 1 || isempty(path)
        path = "./";
    end

    info_file = fullfile(path, "karcher_mean_info.txt");
    data_file = fullfile(path, "karcher_mean_problems.csv");

    meta = readmatrix(info_file);
    T = meta(1);
    D = meta(2);
    n = meta(3);

    ProbFlat = readmatrix(data_file);

    expected_size = [T, D*n*n];
    if ~isequal(size(ProbFlat), expected_size)
        error("ProbFlat has size %s, expected [%d, %d].", ...
            mat2str(size(ProbFlat)), T, D*n*n);
    end

    Prob = zeros(T, D, n, n);

    for t = 1:T
        col = 1;
        for d = 1:D
            Prob(t,d,:,:) = reshape(ProbFlat(t, col:col+n*n-1), n, n);
            col = col + n*n;
        end
    end

    info = struct();
    info.T = T;
    info.D = D;
    info.n = n;
end