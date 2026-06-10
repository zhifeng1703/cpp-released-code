function stgeo_utest(X, opts)


if nargin < 1
    X = stgeo_problem(); 
end
if nargin < 2, opts = struct(); end
if ~isfield(opts, 'MaxIter'), opts.MaxIter = 200; end
if ~isfield(opts, 'AbsTol'), opts.AbsTol = 1e-8; end

if iscell(X), X = X{1}; end

% Preprocess: Q -> U -> QR -> H -> S_init
n = size(X, 1);
p = round(n / 2);
Q = expm(X);
[H, R] = qr(Q(:, 1:p));
H(:, 1:p) = H(:, 1:p) * round(R(1:p, 1:p));
if det(H) < 0
    H(:, p+1) = -H(:, p+1);
end
S_init = real(logm(H));
S_init = 0.5 *(S_init - S_init');

Q_init = expm(S_init);
if norm(Q_init(:, 1:p) - Q(:, 1:p),"fro")/(sqrt(2)* p) > 1e-13
    disp("Error! Bad initialization.")
    return;
end

% Run algorithms
conv_histories = cell(2, 1);

conv_histories{1} = stgeo_zim_logm(S_init, opts);
conv_histories{2} = stgeo_zim_slog(S_init, opts);

% Quick plot

markers = {'o', 's', 'd', '^', 'v', '<', '>', 'p', 'h', '*'};
figure;
for i = 1:2
    % Extract norms (simplified - assumes conv_hist is cell of S matrices)
    if iscell(conv_histories{i})
        p_dim = size(conv_histories{i}{1}, 1)/2;
        norms = zeros(length(conv_histories{i}), 1);
        for j = 1:length(conv_histories{i})
            C = conv_histories{i}{j}(p_dim+1:2*p_dim, p_dim+1:2*p_dim);
            norms(j) = norm(C, 'fro');
        end
    else
        norms = conv_histories{i}.normC;
    end
    semilogy(1:length(norms), norms, 'LineWidth', 2, 'Marker', markers{mod(i-1, length(markers))+1}); hold on;
end
xlabel('Iteration'); ylabel('||C||_F');
legend(["logm", "slog"], 'Location', 'best');
grid on; hold off;
end