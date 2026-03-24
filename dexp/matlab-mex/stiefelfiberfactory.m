function M = stiefelfiberfactory(m, opts)
%STIEFELFIBERFACTORY  Fiber manifold for lifted Stiefel geodesic problem
%
% Point:
%   x.Q ∈ SO(2m)
%   x.S ∈ skew(2m)
%
% Tangent:
%   eta.Z ∈ skew(m)
%
% Retraction:
%   Q+ = [Q1, Q2*exp(tZ)]
%   S+ = mex_sblas_spdiff_logm(Q+, S)
%
% Uses mex_pade_expm cache for dyadic reuse.

    if nargin < 2
        opts = struct();
    end
    if ~isfield(opts,'pade_order')
        opts.pade_order = 13;
    end

    n = 2*m;

    M.name = @() sprintf('Stiefel fiber manifold SO(%d)', n);
    M.dim  = @() m*(m-1)/2;

    M.makepoint   = @makepoint;
    M.retr        = @retr;

    M.proj        = @(x,U) make_tangent(skew(getZ(U)));
    M.tangent     = M.proj;
    M.egrad2rgrad = M.proj;

    M.inner = @(x,eta,zeta) sum(sum(eta.Z .* zeta.Z));
    M.norm  = @(x,eta) norm(eta.Z,'fro');

    M.zerovec = @(x) make_tangent(zeros(m));
    M.lincomb = @lincomb;
    M.transp  = @(x1,x2,eta) eta;

    M.randvec = @(x) make_tangent(skew(randn(m)));

    % ============================================================
    % MAKEPOINT
    % ============================================================

    function x = makepoint(Q, S)

        S = skew(S);

        % Schur + angles
        [R, a] = mex_sblas_ss2schur(S);

        % exponential
        Q = mex_pade_expm(S, opts.pade_order);

        % dexp parameters
        [ParaF, ParaI] = mex_sblas_dexp_para(a);

        x.Q  = Q;
        x.S  = S;

        x.Q1 = Q(:,1:m);
        x.Q2 = Q(:,m+1:end);

        x.R = R;
        x.a = a;

        x.ParaF = ParaF;
        x.ParaI = ParaI;

        x.para = struct('R',R,'A',a,'Fwd',ParaF,'Inv',ParaI);
    end

    % ============================================================
    % RETRACTION
    % ============================================================

    function y = retr(x, eta, t)

        if nargin < 3
            t = 1;
        end

        Z = skew(eta.Z);

        E = cached_pade_expm(Z, t, opts.pade_order);

        Qnew = [x.Q1, x.Q2 * E];

        % branch-following logarithm
        [Snew, Rnew, anew] = mex_sblas_spdiff_logm(Qnew, x.S);

        [ParaF, ParaI] = mex_sblas_dexp_para(anew);

        y.Q  = Qnew;
        y.S  = skew(Snew);

        y.Q1 = Qnew(:,1:m);
        y.Q2 = Qnew(:,m+1:end);

        y.R = Rnew;
        y.a = anew;

        y.ParaF = ParaF;
        y.ParaI = ParaI;

        y.para = struct('R',Rnew,'A',anew,'Fwd',ParaF,'Inv',ParaI);
    end

    % ============================================================
    % LINCOMB
    % ============================================================

    function eta = lincomb(~, varargin)
        switch numel(varargin)
            case 2
                eta = make_tangent(varargin{1} * varargin{2}.Z);
            case 4
                eta = make_tangent(varargin{1} * varargin{2}.Z + ...
                                   varargin{3} * varargin{4}.Z);
            otherwise
                error('Invalid lincomb usage.');
        end
    end

end

% ================================================================
% Cached Padé exponential
% ================================================================

function E = cached_pade_expm(Z, t, p)

    persistent cacheZ cacheP cacheDyadic

    if isempty(cacheZ) || cacheP ~= p || ~isequal(cacheZ,Z)
        cacheZ = Z;
        cacheP = p;
        [~, info] = mex_pade_expm(Z, p);
        cacheDyadic = info.dyadic;
    end

    k = dyadic_level(t);

    if ~isnan(k) && (k+1) <= numel(cacheDyadic)
        E = cacheDyadic{k+1};
    else
        E = mex_pade_expm(t*Z, p);
    end
end

% ================================================================
% Helpers
% ================================================================

function eta = make_tangent(Z)
    eta = struct('Z', skew(Z));
end

function Z = getZ(U)
    if isstruct(U)
        Z = U.Z;
    else
        Z = U;
    end
end

function X = skew(X)
    X = 0.5*(X - X');
end

function k = dyadic_level(t)
    if t <= 0
        k = NaN; return;
    end
    val = -log2(t);
    kr  = round(val);
    if abs(val-kr) < 1e-12
        k = kr;
    else
        k = NaN;
    end
end