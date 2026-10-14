function [en, S] = SS_M_SSAF_SSAFadapt(un, dn, S, sigma_v2)
%SS_M_SSAF_SSAFADAPT SS-M-SSAF using the conventional SSAF normalization.
%   SIGMA_V2 is the fullband background-noise variance. UN and DN are the
%   fullband input and desired signals. EN is the reconstructed fullband
%   error, and S contains the final state and optional diagnostic histories.

% Decompose the fullband signals and form the decimated subband regressors.
[un, dn, w, H, F, m, N, L, sample_count] = prepare_inputs(un, dn, S);
[subband_input, subband_desired] = decompose_signals(un, dn, H, N);
padded_input = [zeros(m - 1, N); subband_input];

% Robust-scale and two-threshold selection parameters.
epsilon = 1e-5;
N_a = 21;
kappa_e = 2.576;
lambda_e = 1 - N / (0.5 * m);
C = 1.483 * (1 + 5 / (N_a - 1));
sigma_v_i2 = sigma_v2 / N;
sigma_e = zeros(N, 1);
error_window = zeros(N, N_a);
window_head = 1;
mu_i = repmat(S.step, N, 1);

en = zeros(1, sample_count);
selected_count = zeros(1, sample_count);
synthesis_state = zeros(L, 1);
[track_msd, unknown_system, unknown_norm, msd] = initialize_msd(S, sample_count);

for n = N:N:sample_count
    if isfield(S, 'unknownsys1') && isfield(S, 'reset_p') && n >= S.reset_p
        unknown_system = S.unknownsys1;
        unknown_norm = unknown_system.' * unknown_system;
    end
    rows = (m - 1 + n):-1:n;
    U = padded_input(rows, :);
    subband_error = subband_desired(n, :).' - U.' * w;
    energy = sum(U.'.^2, 2);
    % SSAF uses one normalization term shared by all subbands.
    normalization = repmat(sqrt(sum(energy) + epsilon), N, 1);

    % Estimate the upper threshold from the smoothed median error scale.
    error_window(:, window_head) = abs(subband_error);
    window_head = mod(window_head, N_a) + 1;
    sigma_e = lambda_e * sigma_e + (1 - lambda_e) * ...
        C * median(error_window, 2);
    xi_i = kappa_e * sigma_e;
    % The MSD-motivated lower threshold excludes low-information subbands.
    zeta_i = 0.5 * mu_i .* energy ./ normalization + ...
        sqrt(2 * sigma_v_i2 / pi);
    selected = zeta_i < xi_i & abs(subband_error) > zeta_i & ...
        abs(subband_error) < xi_i;

    if track_msd
        difference = unknown_system - w;
        msd(n - N + 1:n) = (difference.' * difference) / unknown_norm;
    end
    if n >= S.AdaptStart
        if any(selected)
            % Update the fullband coefficient vector using selected subbands.
            w = w + U(:, selected) * ...
                ((mu_i(selected) ./ normalization(selected)) .* ...
                sign(subband_error(selected)));
        end
        S.iter = S.iter + 1;
    end
    selected_count(n - N + 1:n) = nnz(selected);
    % Reconstruct a fullband error block from the subband errors.
    [en(n - N + 1:n), synthesis_state] = ...
        synthesize_error(subband_error, synthesis_state, F, N);
end

S.coeffs = w;
S.IS = selected_count;
if track_msd
    S.MSD = msd;
end
end

function [un, dn, w, H, F, m, N, L, sample_count] = prepare_inputs(un, dn, S)
un = un(:).';
dn = dn(:).';
if numel(un) ~= numel(dn)
    error('un and dn must have the same length.');
end
w = S.coeffs(:);
H = S.analysis;
F = S.synthesis;
[L, N] = size(H);
m = numel(w);
sample_count = numel(un);
if mod(m, N) ~= 0
    error('The adaptive-filter length must be a multiple of the subband count.');
end
end

function [u_sub, d_sub] = decompose_signals(un, dn, H, N)
sample_count = numel(un);
u_sub = zeros(sample_count, N);
d_sub = zeros(sample_count, N);
for band = 1:N
    u_sub(:, band) = filter(H(:, band), 1, un).';
    d_sub(:, band) = filter(H(:, band), 1, dn).';
end
end

function [enabled, unknown_system, unknown_norm, msd] = initialize_msd(S, count)
enabled = isfield(S, 'unknownsys');
unknown_system = zeros(numel(S.coeffs), 1);
unknown_norm = 1;
msd = zeros(1, count);
if enabled
    unknown_system = S.unknownsys(:);
    unknown_norm = unknown_system.' * unknown_system;
end
end

function [output_block, state] = synthesize_error(error, state, F, N)
state = F * error + state;
output_block = state(1:N).';
state = [state(N + 1:end); zeros(N, 1)];
end
