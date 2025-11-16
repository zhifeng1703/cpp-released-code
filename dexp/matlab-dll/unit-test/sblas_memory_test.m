function sblas_memory_test(dims)
    if nargin < 1, dims = 3:100; end
        
    srcPath = fullfile(fileparts(mfilename('fullpath')), 'src');
    addpath(srcPath);

    % --- Load DLL if not already loaded ---
    libPath = fullfile(fileparts(mfilename('fullpath')), 'lib');
    dll = fullfile(libPath, 'skewblas.dll');
    hdr = fullfile(libPath, 'skewblas_api.h');

    if ~libisloaded('skewblas')
        fprintf('Loading library from %s\n', libPath);
        loadlibrary(dll, hdr);
    else
        fprintf('Library already loaded.\n');
    end

    i = 1;
    for d = dims
        S = rand(d);
        S = 100 * (S - S') / 2;
        Q = sblas_expm(S);
        m = memory;
        fprintf('%3d: MATLAB used %.1f MB\n', i, m.MemUsedMATLAB/1e6);
        i = i+1;
    end
end