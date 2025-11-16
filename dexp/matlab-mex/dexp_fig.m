
function [fig_fwd, fig_inv] = dexp_fig(dims, elapsed_fwd, elapsed_inv, opts)
% SAVE_DEXP_FIG  Redraw forward & inverse Dexp timing figures and save as tight PDFs.
%
%   save_dexp_fig(elapsed_fwd, elapsed_inv)
%
%   Inputs:
%       elapsed_fwd : (#dims x 4) timing matrix for forward Dexp
%       elapsed_inv : (#dims x 4) timing matrix for inverse Dexp
%
%   This function:
%       - infers the dimension vector
%       - redraws both figures with your exact styling
%       - saves them as:
%             dexp_fwd_matlab.pdf
%             dexp_inv_matlab.pdf
%         using tight, vector PDF export.
% ---- Default options ----
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end
    if ~isfield(opts, 'savefig'), opts.savefig = false; end

    % ---- Figure visibility ----
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

    margin = 0.40;   % 10 mm = 0.40 in approx


    %% ------------------------------------------------------------
    % FIGURE 1 — Forward timing
    % ------------------------------------------------------------
    fig_fwd = figure('Visible', vis);

    plot(dims, elapsed_fwd(:,1), '-o', 'LineWidth',1.8, 'DisplayName','sBLAS Dexp (fwd)'); hold on;
    plot(dims, elapsed_fwd(:,2), '-x', 'LineWidth',1.8, 'DisplayName','Padé p=3 (fwd)');
    plot(dims, elapsed_fwd(:,3), '-s', 'LineWidth',1.8, 'DisplayName','Padé p=13 (fwd)');
    plot(dims, elapsed_fwd(:,4), '-^', 'LineWidth',1.8, 'DisplayName','Daletskii–Krein (fwd)');
    hold off;

    title('Forward Dexp Timing');
    xlabel('Matrix Dimension n');
    ylabel('Time (seconds)');
    legend('Location','northwest');
    grid on;


    if opts.savefig
        fig_fwd.Units = 'inches';
        pos = fig_fwd.Position;
        fig_fwd.Position = [pos(1)-margin, pos(2)-margin, pos(3)+2*margin, pos(4)+2*margin];

        exportgraphics(fig_fwd, fullfile(figpath, 'dexp_fwd_matlab.pdf'), ...
            'ContentType','vector', 'BackgroundColor','none');
    end


    %% ------------------------------------------------------------
    % FIGURE 2 — Inverse timing
    % ------------------------------------------------------------
    fig_inv = figure('Visible', vis);

    plot(dims, elapsed_inv(:,1), '-o', 'LineWidth',1.8, 'DisplayName','sBLAS Dexp (inv)'); hold on;
    plot(dims, elapsed_inv(:,2), '-x', 'LineWidth',1.8, 'DisplayName','Padé p=1 (inv)');
    plot(dims, elapsed_inv(:,3), '-s', 'LineWidth',1.8, 'DisplayName','Padé p=7 (inv)');
    plot(dims, elapsed_inv(:,4), '-^', 'LineWidth',1.8, 'DisplayName','Daletskii–Krein (inv)');
    hold off;

    title('Inverse Dexp Timing');
    xlabel('Matrix Dimension n');
    ylabel('Time (seconds)');
    legend('Location','northwest');
    grid on;


    if opts.savefig
        fig_inv.Units = 'inches';
        pos = fig_inv.Position;
        fig_inv.Position = [pos(1)-margin, pos(2)-margin, pos(3)+2*margin, pos(4)+2*margin];

        exportgraphics(fig_inv, fullfile(figpath, 'dexp_inv_matlab.pdf'), ...
            'ContentType','vector', 'BackgroundColor','none');
    end


    if opts.savefig
        fprintf('Saved figures to:\n  %s\n  %s\n', ...
            fullfile(figpath,'dexp_fwd_matlab.pdf'), ...
            fullfile(figpath,'dexp_inv_matlab.pdf'));
    end
end
