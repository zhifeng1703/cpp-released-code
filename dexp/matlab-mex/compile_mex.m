function compile_mex(mex_dir)
% COMPILE_MEX  Check and compile all required MEX files.
%
%   compile_mex()
%   compile_mex(mex_dir)
%
%   INPUT
%       mex_dir (optional) – folder containing build_mex.m
%                            default = './mex/'
%
%   The function:
%       • checks whether required MEX files exist in mex_dir
%       • if any are missing, calls mex_dir/build_mex.m
%       • all compiled binaries stay inside mex_dir
%
%   Example:
%       compile_mex;                        % default ./mex
%       compile_mex('C:/projects/skewblas/mex');

    if nargin < 1 || isempty(mex_dir)
        % default: ./mex/ relative to caller
        caller_path = fileparts(mfilename('fullpath'));
        mex_dir = fullfile(caller_path, 'mex');
    end

    if ~isfolder(mex_dir)
        error('compile_mex:NoFolder', ...
            'MEX directory does not exist: %s', mex_dir);
    end

    addpath(mex_dir);

    fprintf('\n=== Checking MEX folder: %s ===\n', mex_dir);

    required_mex = {
        'mex_sblas_expm'
        'mex_sblas_dexp'
        'mex_sblas_dexp_para'
        'mex_dk_dexp'
        'mex_dk_dexp_para'
        'mex_pade_dexpfwd'
        'mex_pade_dexpinv'
        'mex_sblas_ss2schur'
        'mex_sblas_so2schur'
    };

    need_rebuild = false;

    for k = 1:numel(required_mex)
        mex_file = fullfile(mex_dir, [required_mex{k} '.' mexext]);
        if exist(mex_file, 'file') ~= 3
            fprintf('   Missing: %s\n', mex_file);
            need_rebuild = true;
        end
    end

    if ~need_rebuild
        fprintf('All required MEX files compiled.\n');
        return;
    end

    % Build script must be inside this folder
    build_script = fullfile(mex_dir, 'build_mex.m');
    if exist(build_script, 'file') ~= 2
        error('compile_mex:MissingBuildScript', ...
            'Cannot find build_mex.m at: %s', build_script);
    end

    fprintf('Some MEX files are missing — compiling now...\n');
    old_dir = pwd;
    cd(mex_dir);

    try
        run('build_mex.m');
    catch ME
        cd(old_dir);
        rethrow(ME);
    end

    cd(old_dir);
    fprintf('=== MEX build complete. ===\n\n');
end
