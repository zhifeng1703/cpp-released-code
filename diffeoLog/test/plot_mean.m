function [plt, PlogAngles, DlogAngles, DataAngles] = plot_mean(Mp, Md, Prob)
%PLOT_MEAN Plot principal-log mean, diffeomorphic-log mean, and data angles.
% Inputs:
%   Mp   : T x 5 x 5 array.
%          squeeze(Mp(t,:,:)) is the principal-log mean, already
%          represented as a 5 x 5 skew-symmetric matrix.
%
%   Md   : T x 5 x 5 array.
%          squeeze(Md(t,:,:)) is the diffeomorphic-log mean, already
%          represented as a 5 x 5 skew-symmetric matrix.
%
%   Prob : T x D x 5 x 5 array.
%          squeeze(Prob(t,d,:,:)) is the d-th skew-symmetric data
%          logarithm in the t-th problem.
%
% Outputs:
%   plt         : figure handle
%   PlogAngles  : T x 2 array of canonical principal-log mean angles
%   DlogAngles  : T x 2 array of canonical diffeomorphic-log mean angles
%   DataAngles  : T x D x 2 array of canonical data angle coordinates

    %% Basic checks for Mp and Md

    if ndims(Mp) ~= 3
        error('Expected Mp to be a 3D array of size T x 5 x 5.');
    end

    if ndims(Md) ~= 3
        error('Expected Md to be a 3D array of size T x 5 x 5.');
    end

    [T, n1, n2] = size(Mp);
    [Td, dn1, dn2] = size(Md);

    if n1 ~= 5 || n2 ~= 5
        error('Expected Mp to have size T x 5 x 5, but got %s.', mat2str(size(Mp)));
    end

    if Td ~= T || dn1 ~= 5 || dn2 ~= 5
        error('Expected Md to have size T x 5 x 5 matching Mp, but got %s.', mat2str(size(Md)));
    end

    %% Basic checks for Prob

    if ndims(Prob) ~= 4
        error('Expected Prob to be a 4D array of size T x D x 5 x 5.');
    end

    [Tp, D, pn1, pn2] = size(Prob);

    if Tp ~= T
        error('Mean arrays and Prob have inconsistent time dimensions: means have T=%d, Prob has T=%d.', T, Tp);
    end

    if pn1 ~= 5 || pn2 ~= 5
        error('Expected Prob to have size T x D x 5 x 5, but got %s.', mat2str(size(Prob)));
    end

    n = 5;
    m = floor(n / 2);

    if m < 2
        error('Need at least two angles to make the 2D angle plot.');
    end

    %% Compute canonical angle coordinates

    PlogAngles = zeros(T, 2);
    DlogAngles = zeros(T, 2);
    DataAngles = zeros(T, D, 2);

    for t = 1:T
        Mt = squeeze(Mp(t,:,:));
        Mt = 0.5 * (Mt - Mt');
        [~, alpha] = mex_sblas_ss2schur(Mt);
        beta = canonical_align_angles(alpha, n);
        PlogAngles(t, :) = beta(1:2).';

        Mt = squeeze(Md(t,:,:));
        Mt = 0.5 * (Mt - Mt');
        [~, alpha] = mex_sblas_ss2schur(Mt);
        beta = canonical_align_angles(alpha, n);
        DlogAngles(t, :) = beta(1:2).';

        for d = 1:D
            Dt = squeeze(Prob(t,d,:,:));
            Dt = 0.5 * (Dt - Dt');

            [~, alpha] = mex_sblas_ss2schur(Dt);
            beta = canonical_align_angles(alpha, n);

            DataAngles(t, d, :) = beta(1:2);
        end
    end

    %% Scatter plot

    plt = figure('Color', 'w');
    hold on;
    box on;
    axis equal;

    h_data = gobjects(D, 1);

    for d = 1:D
        A = squeeze(DataAngles(:, d, :));

        h_data(d) = scatter(A(:,1), A(:,2), ...
            1, 'g', 'filled', ...
            'MarkerFaceAlpha', 0.35, ...
            'MarkerEdgeAlpha', 0.35);
    end

    h_plog = scatter(PlogAngles(:,1), PlogAngles(:,2), ...
        8, 'r', 'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.75);

    h_dlog = scatter(DlogAngles(:,1), DlogAngles(:,2), ...
        1, 'b', 'filled', ...
        'MarkerFaceAlpha', 0.75, ...
        'MarkerEdgeAlpha', 0.75);

    xlabel('\alpha_1', 'Interpreter', 'tex');
    ylabel('\alpha_2', 'Interpreter', 'tex');
    title('Karcher mean experiment');

    % xlim([1.3, 4.5]);
    % ylim([0.0, 3.2]);

    legend([h_data(1), h_plog, h_dlog], ...
        {'Data Points', ...
         'Mean by Principal Logarithm', ...
         'Mean by Diffeomorphic Logarithm'}, ...
        'Location', 'best');

    hold off;
end


function beta = canonical_align_angles(alpha, n)
%CANONICAL_ALIGN_ANGLES Canonically order and sign-align Schur angles.
%
% Input:
%   alpha : raw Schur angles alpha_i
%   n     : matrix dimension
%
% Output:
%   beta  : ordered signed angles

    alpha = alpha(:);
    m = floor(n / 2);

    if length(alpha) < m
        error('alpha has length %d, but expected at least %d.', ...
            length(alpha), m);
    end

    alpha = alpha(1:m);

    theta = mod(alpha - pi, 2*pi) - pi;

    tol = 10 * eps;
    theta(abs(theta + pi) < tol) = pi;

    [~, p] = sort(abs(theta), 'descend');

    theta_p = theta(p);
    delta = ones(m, 1);

    for i = 1:(m-1)
        if theta_p(i) < 0
            theta_p(i)   = -theta_p(i);
            theta_p(i+1) = -theta_p(i+1);

            delta(p(i))   = -delta(p(i));
            delta(p(i+1)) = -delta(p(i+1));
        end
    end

    if mod(n, 2) == 1 && theta_p(m) < 0
        theta_p(m) = -theta_p(m);
        delta(p(m)) = -delta(p(m));
    end

    beta = delta(p) .* alpha(p);
end