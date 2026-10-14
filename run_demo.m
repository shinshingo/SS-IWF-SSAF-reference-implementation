function results = run_demo()
%RUN_DEMO Run the six proposed algorithms with an abrupt system change.
%   RESULTS = RUN_DEMO generates a colored input and a desired response
%   contaminated by background Gaussian and Bernoulli-Gaussian impulsive
%   noise. The unknown system changes from w to -w halfway through the run.

% Simulation settings used by all six algorithms.
rng(1);
sample_count = 4000000;
m = 512;
N = 8;
L = 64;
change_sample = sample_count / 2;

% Generate the colored input and the two unknown-system responses.
input_signal = filter(1, [1, -0.9], randn(1, sample_count));
first_system = randn(m, 1);
first_system = first_system / norm(first_system);
second_system = -first_system;
first_output = filter(first_system, 1, input_signal);
second_output = filter(second_system, 1, input_signal);

% Normalize the input and noise-free responses by the desired-signal scale.
normalization_scale = std(first_output);
input_signal = input_signal / normalization_scale;
first_output = first_output / normalization_scale;
second_output = second_output / normalization_scale;
noise_free_output = [first_output(1:change_sample), ...
    second_output(change_sample + 1:end)];

% Add 20-dB background Gaussian noise and BG impulsive noise.
sigma_d2 = var(noise_free_output);
sigma_vg2 = sigma_d2 * 10^(-20 / 10);
background_noise = sqrt(sigma_vg2) * randn(1, sample_count);
P_r = 0.01;
sigma_eta2 = 1000 * sigma_d2;
impulse_mask = rand(1, sample_count) < P_r;
eta_noise = impulse_mask .* sqrt(sigma_eta2) .* randn(1, sample_count);
desired_signal = noise_free_output + background_noise + eta_noise;

% Initialize the six algorithms with their respective step sizes.
S1 = initialize_state(0.00118);
S2 = initialize_state(0.00018);
S3 = initialize_state(0.0016);
S4 = initialize_state(0.00024);
S5 = initialize_state(0.0025);
S6 = initialize_state(0.0025);

% Run the fixed-step-size and band-dependent VSS algorithms.
[~, results.m_ssaf] = M_SSAF_SSAFadapt(input_signal, desired_signal, S1);
[~, results.m_iwf_ssaf] = M_SSAF_IWF_SSAFadapt( ...
    input_signal, desired_signal, S2);
[~, results.ss_m_ssaf] = SS_M_SSAF_SSAFadapt( ...
    input_signal, desired_signal, S3, sigma_vg2);
[~, results.ss_m_iwf_ssaf] = SS_M_SSAF_IWF_SSAFadapt( ...
    input_signal, desired_signal, S4, sigma_vg2);
[~, results.vss_ss_m_ssaf] = VSS_SS_M_SSAF_SSAFadapt( ...
    input_signal, desired_signal, S5, sigma_d2, sigma_vg2, 1, 2);
[~, results.vss_ss_m_iwf_ssaf] = VSS_SS_M_SSAF_IWF_SSAFadapt( ...
    input_signal, desired_signal, S6, sigma_d2, sigma_vg2, 1, 2);
results.system_change_sample = change_sample;
results.P_r = P_r;
results.sigma_eta2 = sigma_eta2;

% Plot the one-run NMSD curves and mark the system-change sample.
figure;
layout = tiledlayout(2, 1);
nexttile;
plot(10 * log10(results.m_ssaf.MSD), 'LineWidth', 1.2);
hold on;
plot(10 * log10(results.m_iwf_ssaf.MSD), 'LineWidth', 1.2);
plot(10 * log10(results.ss_m_ssaf.MSD), 'LineWidth', 1.2);
plot(10 * log10(results.ss_m_iwf_ssaf.MSD), 'LineWidth', 1.2);
grid on;
ylabel('NMSD (dB)');
legend('M-SSAF (SSAF)', 'M-SSAF (IWF-SSAF)', ...
    'SS-M-SSAF (SSAF)', 'SS-M-SSAF (IWF-SSAF)', 'Location', 'best');

nexttile;
marker_indices = round(linspace(1, sample_count, 40));
plot(10 * log10(results.vss_ss_m_ssaf.MSD), '-o', ...
    'LineWidth', 1.2, 'MarkerIndices', marker_indices(1:2:end), 'MarkerSize', 5);
hold on;
plot(10 * log10(results.vss_ss_m_iwf_ssaf.MSD), '-s', ...
    'LineWidth', 1.2, 'MarkerIndices', marker_indices(2:2:end), 'MarkerSize', 5);
grid on;
xlabel('Fullband input sample');
ylabel('NMSD (dB)');
legend('VSS-SS-M-SSAF (SSAF)', 'VSS-SS-M-SSAF (IWF-SSAF)', ...
    'Location', 'best');
title(layout, 'Tracking under BG impulsive noise after a system change');

    function state = initialize_state(mu)
        % Attach the unknown systems only to evaluate the NMSD trajectory.
        state = SSAFinit(zeros(m, 1), mu, N, L);
        state.unknownsys = first_system;
        state.unknownsys1 = second_system;
        state.reset_p = change_sample;
    end
end
