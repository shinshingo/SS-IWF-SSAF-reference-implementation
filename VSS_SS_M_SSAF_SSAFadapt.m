function [en, S] = VSS_SS_M_SSAF_SSAFadapt(un, dn, S, sigma_d2, sigma_v2, reset_type, kappa)
%VSS_SS_M_SSAF_SSAFADAPT VSS-SS-M-SSAF with SSAF normalization.
%   SIGMA_D2 and SIGMA_V2 are the noise-free desired-signal and fullband
%   background-noise variances. RESET_TYPE = 1 enables reset detection, and
%   KAPPA controls the smoothing factors. S is initialized by SSAFinit.

if nargin < 6
    reset_type = 1;
end
if nargin < 7
    kappa = 2;
end

% Decompose the fullband signals and form the decimated subband regressors.
[un, dn, w, H, F, m, N, L, sample_count] = prepare_inputs(un, dn, S);
[subband_input, subband_desired] = decompose_signals(un, dn, H, N);
padded_input = [zeros(m - 1, N); subband_input];

% Robust-scale, VSS, and two-threshold selection parameters.
epsilon = 1e-5;
N_a = 21;
kappa_e = 2.576;
lambda_e = 1 - N / (kappa * m);
lambda_mu = lambda_e;
C = 1.483 * (1 + 5 / (N_a - 1));
sigma_v_i2 = sigma_v2 / N;
sigma_e = zeros(N, 1);
error_window = zeros(N, N_a);
window_head = 1;
sigma_u2 = var(un);
% Initialize every subband at the same bounded maximum step size.
mu_max = min(1, sqrt(sigma_d2 / (sigma_u2 * m)));
mu_i = repmat(mu_max, N, 1);
beta_i = ones(N, 1);
lambda_beta = 0.8;

% Reset-control window and decision parameters from the manuscript.
V_T = m;
V_D = round(0.75 * V_T);
normalized_error = zeros(V_T, 1);
control_value = 0;
previous_control = 0;
eta = 0.005;

en = zeros(1, sample_count);
selected_count = zeros(1, sample_count);
step_history = zeros(1, sample_count);
control_history = zeros(1, sample_count);
reset_history = false(1, sample_count);
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

    % Estimate the M-estimate upper threshold and identify reliable errors.
    error_window(:, window_head) = abs(subband_error);
    window_head = mod(window_head, N_a) + 1;
    sigma_e = lambda_e * sigma_e + (1 - lambda_e) * ...
        C * median(error_window, 2);
    xi_i = kappa_e * sigma_e;
    reliable = abs(subband_error) < xi_i;

    if track_msd
        difference = unknown_system - w;
        msd(n - N + 1:n) = (difference.' * difference) / unknown_norm;
    end
    if n >= S.AdaptStart
        % Reset the VSS state when the normalized-error control rises sharply.
        reset_now = reset_type == 1 && ...
            (control_value - previous_control) / (mean(mu_i) + epsilon) > eta;
        if reset_now
            mu_i(:) = mu_max;
            beta_i(:) = 1;
            reset_history(n) = true;
        else
            beta_i(reliable) = lambda_beta * beta_i(reliable) + ...
                (1 - lambda_beta) * abs(subband_error(reliable));
            % Update each subband step size using its reliable error level.
            psi_i = normalization ./ (energy + epsilon) .* ...
                (beta_i - sqrt(2 * sigma_v_i2 / pi));
            active = psi_i > 0;
            mu_i(active) = lambda_mu * mu_i(active) + ...
                (1 - lambda_mu) * min(psi_i(active), mu_i(active));
        end
        previous_control = control_value;

        % Apply the MSD-motivated lower threshold and M-estimate upper threshold.
        zeta_i = 0.5 * mu_i .* energy ./ normalization + ...
            sqrt(2 * sigma_v_i2 / pi);
        selected = zeta_i < xi_i & abs(subband_error) > zeta_i & reliable;
        if any(selected)
            % Update the fullband coefficient vector using selected subbands.
            w = w + U(:, selected) * ...
                ((mu_i(selected) ./ normalization(selected)) .* ...
                sign(subband_error(selected)));
        end
        S.iter = S.iter + 1;
    else
        selected = false(N, 1);
    end
    selected_count(n - N + 1:n) = nnz(selected);
    % Reconstruct a fullband error block from the subband errors.
    [en(n - N + 1:n), synthesis_state] = ...
        synthesize_error(subband_error, synthesis_state, F, N);

    if reset_type == 1
        % Retain the smallest normalized errors to reduce outlier influence.
        normalized_error = [abs(en(n - N + 1:n)).' / ...
            sqrt(sum(energy) + epsilon); normalized_error(1:end - N)];
        if mod(n, V_T) == 0
            ordered_error = sort(normalized_error);
            retained = ordered_error(1:V_T - V_D);
            % Update the reset control once per fullband window.
            control_value = mean(retained.^2);
        end
    end
    step_history(n - N + 1:n) = mean(mu_i);
    control_history(n - N + 1:n) = control_value;
end

S.coeffs = w;
S.IS = selected_count;
S.mu = step_history;
S.ctrl = control_history;
S.reset_event = reset_history;
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
