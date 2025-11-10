function elapsed_times = time_schur(dims)
    if nargin < 1, dims = 3:50; end

    % --- Add source folder to MATLAB path ---
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

    elapsed_times = zeros(length(dims), 2);
    loop = 20;

    ind = 1;
    for dim = dims
        S = 200 * rand(dim);
        opt.verbose = true;
        test_schur(dim, S, opt);

        opt.verbose = false;
        for i = 1:loop
            [t1, t2] = test_schur(dim, S, opt);
            elapsed_times(ind, 1) = elapsed_times(ind, 1) + t1;
            elapsed_times(ind, 2) = elapsed_times(ind, 2) + t2;
        end
        elapsed_times(ind, :) = elapsed_times(ind, :) / loop;
        ind = ind + 1;
    end

    %if libisloaded('skewblas')
    %    unloadlibrary('skewblas');
    %end
end