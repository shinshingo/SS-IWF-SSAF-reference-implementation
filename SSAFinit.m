function S = SSAFinit(w0, mu, N, L)
%SSAFINIT Initialize an SSAF state with a cosine-modulated filter bank.
%   S = SSAFINIT(W0, MU, N, L) initializes the fullband coefficient vector
%   W0, step size MU, N subbands, and an L-tap analysis filter bank.

validateattributes(w0, {'numeric'}, {'vector', 'real', 'finite'});
validateattributes(mu, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(N, {'numeric'}, {'scalar', 'integer', '>=', 2});
validateattributes(L, {'numeric'}, {'scalar', 'integer', '>=', N});

% Reuse the filter bank when N and L are unchanged across algorithms.
persistent cached_N cached_L cached_H cached_F
if isempty(cached_H) || cached_N ~= N || cached_L ~= L
    p = design_prototype(N, L);
    l = (0:L-1).';
    H = zeros(L, N);
    % Construct each analysis filter using the manuscript equation.
    for i = 0:N-1
        H(:, i + 1) = 2 * p .* cos((pi / N) * (i + 0.5) * ...
            (l - (L - 1) / 2) + (-1)^i * pi / 4);
    end
    cached_N = N;
    cached_L = L;
    cached_H = H;
    cached_F = flipud(H);
end

S.coeffs = w0(:);
S.step = mu;
S.analysis = cached_H;
S.synthesis = cached_F;
S.iter = 0;
S.AdaptStart = L + numel(w0);
end

function p = design_prototype(N, L)
% Select the passband edge by minimizing the power-complementarity error.
stop_edge = 1 / N;
objective = @(pass_edge) power_complementarity_error( ...
    firpm(L - 1, [0, pass_edge, stop_edge, 1], [1, 1, 0, 0], [5, 1]), N);
options = optimset('Display', 'off', 'TolX', 1e-7);
pass_edge = fminbnd(objective, 1 / (8 * N), 0.9 / N, options);
p = firpm(L - 1, [0, pass_edge, stop_edge, 1], ...
    [1, 1, 0, 0], [5, 1]).';
end

function error_value = power_complementarity_error(p, N)
omega = linspace(0, pi / N, 512);
response = freqz(p, 1, omega);
paired_response = freqz(p, 1, pi / N - omega);
error_value = max(abs(abs(response).^2 + abs(paired_response).^2 - 1));
end
