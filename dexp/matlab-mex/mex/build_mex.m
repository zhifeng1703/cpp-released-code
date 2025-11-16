function build_mex
    fprintf('\n=== Building skewblas MEX functions ===\n\n');

    root_mex = fileparts(mfilename('fullpath'));
    root_src = fullfile(root_mex, '..', 'src');

    % --------------------------
    % Collect all .cpp source files in src/
    % --------------------------
    S = dir(fullfile(root_src, '*.cpp'));
    src_cpp = fullfile(root_src, {S.name});   % cell array of strings

    % --------------------------
    % Common MEX options
    % --------------------------
    mex_opts = {
        '-v'
        '-R2017b'
        '-DMATLAB_MEX_BUILD'
        ['-I' root_src]
        '-lmwblas'
        '-lmwlapack'
        '-outdir' 
        root_mex
    };

    % --------------------------
    % MEX gateways to build
    % --------------------------
    gateways = {
        'mex_sblas_expm.cpp'        
        'mex_sblas_ss2schur.cpp'      
        'mex_sblas_so2schur.cpp'
        'mex_sblas_dexp_para.cpp'
        'mex_sblas_dexp.cpp'
        'mex_dk_dexp_para.cpp'
        'mex_dk_dexp.cpp'
        'mex_pade_expm.cpp'
        'mex_pade_dexpfwd.cpp'
        'mex_pade_dexpinv.cpp'
    };

    for k = 1:numel(gateways)
        gw = fullfile(root_mex, gateways{k});
        fprintf('Compiling %s ...\n', gateways{k});

        mex(mex_opts{:}, gw, src_cpp{:});
    end

    fprintf('\n=== skewblas MEX build COMPLETE ===\n');
end
