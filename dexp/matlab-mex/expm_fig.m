function [fig_elapsed, fig_error] = expm_fig(dims, elapsed_times, errors, opts)
%PLOT_EXPM  Draw and (optionally) save figures for EXPM benchmark.
%
%   plot_expm(dims, elapsed_times, errors)
%   plot_expm(dims, elapsed_times, errors, opts)
%
%   opts.savefig   = true/false (default: false)
%   opts.outdir    = target folder (default: ./figures)

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'savefig'), opts.savefig = true; end
    if ~isfield(opts,'verbose'),  opts.verbose = true; end

    if opts.verbose
        vis = 'on';
    else
        vis = 'off';
    end

    % ---- Output directory ----
    figpath = fullfile(fileparts(mfilename('fullpath')), 'figures');
    if opts.savefig && ~exist(figpath, 'dir')
        mkdir(figpath);
    end

    %% ------------------------------------------------------------
    % FIGURE 1 — Elapsed Times
    % ------------------------------------------------------------
    fig_elapsed = figure('Visible', vis);
    plot(dims, elapsed_times(:,1), '-o', 'LineWidth', 1.8, ...
         'DisplayName','sBLAS Expm'); hold on;
    plot(dims, elapsed_times(:,2), '-x', 'LineWidth', 1.8, ...
         'DisplayName','Pade3 Expm');
    plot(dims, elapsed_times(:,3), '-s', 'LineWidth', 1.8, ...
         'DisplayName','Pade13 Expm');
    plot(dims, elapsed_times(:,4), '-^', 'LineWidth', 1.8, ...
         'DisplayName','Built-in Expm');
    hold off;

    title('Elapsed Time for Matrix Exponential');
    xlabel('Matrix Dimension n');
    ylabel('Time (sec)');
    legend('Location','northwest');
    grid on;

    if opts.savefig
        exportgraphics(fig_elapsed, fullfile(figpath, 'expm_elapsed_matlab.pdf'), ...
            'ContentType','vector','BackgroundColor','none');
    end


    %% ------------------------------------------------------------
    % FIGURE 2 — Error (log-scale)
    % ------------------------------------------------------------
    fig_error = figure('Visible',vis);
    semilogy(dims, errors(:,1), '-o', 'LineWidth',1.8, ...
             'DisplayName','Builtin vs sBLAS'); hold on;
    semilogy(dims, errors(:,2), '-x', 'LineWidth',1.8, ...
             'DisplayName','Builtin vs Pade3');
    semilogy(dims, errors(:,3), '-s', 'LineWidth',1.8, ...
             'DisplayName','Builtin vs Pade13');
    hold off;

    title('Matrix Exponential Error (Frobenius Norm)');
    xlabel('Matrix Dimension n');
    ylabel('‖Q_builtin - Q_method‖_F');
    legend('Location','northwest');
    grid on;

    if opts.savefig
        exportgraphics(fig_error, fullfile(figpath, 'expm_error_matlab.pdf'), ...
            'ContentType','vector','BackgroundColor','none');
    end
end
