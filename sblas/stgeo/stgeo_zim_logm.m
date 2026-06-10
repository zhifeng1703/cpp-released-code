function conv_hist = stgeo_zim_logm(S_init, opts)

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'MaxIter'), opts.MaxIter = 200; end
if ~isfield(opts, 'AbsTol'), opts.AbsTol = 1e-8; end

n = size(S_init,1);
p = round(n/2);
MaxIter = opts.MaxIter;
AbsTol = opts.AbsTol;

S = S_init;
Q = expm(S);
conv_hist = cell(MaxIter, 1);

for k = 1:MaxIter
    conv_hist{k} = S;
    
    C = S(p+1:2*p, p+1:2*p);
    err = norm(C, 'fro');
    if err < AbsTol
        conv_hist = conv_hist(1:k);
        break;
    end

    Q(:,p+1:n) = Q(:,p+1:n) * expm(-C);
    S = real(logm(Q));
    S = 0.5 * (S - S');
end
end