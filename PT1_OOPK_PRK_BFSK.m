% Digital Communications Final Project
% Part I: Performance of OOK, PRK, and BFSK (Peak Eb/N0)

clear; clc; close all;

%% 1. Simulation Parameters
N_bits = 1e6; 
EbN0_dB = 0:3:60; 

BER_OOK_manual = []; BER_PRK_manual = []; BER_BFSK_manual = [];

%% 2. Generate Random Binary Data
data = randi([0 1], 1, N_bits);

%% Manual Implementation Loop (Peak Eb/N0)
for ebno = EbN0_dB
    EbN0_lin = 10^(ebno/10);
    N0 = 1 / EbN0_lin; 
    noise_variance = N0 / 2; % Baseband noise variance per dimension
    
    % --- PRK (BPSK) ---
    s_prk = 2*data - 1; 
    noise_prk = sqrt(noise_variance) * randn(1, N_bits); 
    detected_prk = (s_prk + noise_prk) > 0;
    BER_PRK_manual = [BER_PRK_manual, sum(data ~= detected_prk) / N_bits];
    
    % --- BFSK (Orthogonal) ---
    s_bfsk = zeros(1, N_bits);
    s_bfsk(data == 0) = 1; s_bfsk(data == 1) = 1i; 
    noise_bfsk = sqrt(noise_variance) * (randn(1, N_bits) + 1i*randn(1, N_bits));
    detected_bfsk = imag(s_bfsk + noise_bfsk) > real(s_bfsk + noise_bfsk);
    BER_BFSK_manual = [BER_BFSK_manual, sum(data ~= detected_bfsk) / N_bits];

    % --- OOK (On-Off Keying) ---
    s_ook = data; 
    noise_ook = sqrt(noise_variance) * randn(1, N_bits);
    detected_ook = (s_ook + noise_ook) > 0.5;
    BER_OOK_manual = [BER_OOK_manual, sum(data ~= detected_ook) / N_bits];
end

%% Theoretical Formulas 
% PRK Coherent Theory: Q(sqrt(2*Eb/N0))
BER_PRK_theory = qfunc(sqrt(2 * 10.^(EbN0_dB/10)));

% BFSK Coherent Theory: Q(sqrt(Eb/N0))
BER_BFSK_theory = qfunc(sqrt(10.^(EbN0_dB/10)));

% OOK Peak Theory: Q(sqrt(0.5 * Eb/N0))  <-- REPLACE WITH THIS LINE
BER_OOK_theory = qfunc(sqrt(0.5 * 10.^(EbN0_dB/10)));


%% Built-in Verification 
BER_PRK_builtin = zeros(1, length(EbN0_dB));

for i = 1:length(EbN0_dB)
    disp(['Processing Built-in PRK Eb/N0: ', num2str(EbN0_dB(i)), ' dB...']); 
    s_prk_built = pskmod(data, 2); 
    r_prk_built = awgn(s_prk_built, EbN0_dB(i), 'measured');
    det_prk_built = pskdemod(r_prk_built, 2);
    [~, BER_PRK_builtin(i)] = biterr(data, det_prk_built);
end

%% Plot Figure 1: Manual Only
figure('Name', 'Figure 1: Manual Simulation', 'Position', [100, 100, 800, 600]);
semilogy(EbN0_dB, BER_OOK_manual, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(EbN0_dB, BER_PRK_manual, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(EbN0_dB, BER_BFSK_manual, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Part I - BER vs Peak E_b/N_0: OOK, PRK, BFSK');
xlabel('Peak E_b/N_0 (dB)'); ylabel('Bit Error Rate (BER)');
legend('OOK (manual)', 'PRK (manual)', 'BFSK (manual)');
ylim([1e-6 1]);

%% Plot Figure 2: Manual vs Theory vs Built-in (Combined)
figure('Name', 'Figure 2: Verification', 'Position', [150, 150, 800, 600]);

% 1. Manual Simulations (Solid lines)
semilogy(EbN0_dB, BER_OOK_manual, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_dB, BER_PRK_manual, 'r-', 'LineWidth', 2);
semilogy(EbN0_dB, BER_BFSK_manual, 'g-', 'LineWidth', 2);

% 2. Theoretical Curves (Dashed lines)
semilogy(EbN0_dB, BER_OOK_theory, 'b--', 'LineWidth', 2);
semilogy(EbN0_dB, BER_PRK_theory, 'r--', 'LineWidth', 2);
semilogy(EbN0_dB, BER_BFSK_theory, 'g--', 'LineWidth', 2);

% 3. Built-in PRK Verification (Dotted line with markers)
semilogy(EbN0_dB, BER_PRK_builtin, 'k:x', 'MarkerSize', 10, 'LineWidth', 2);

grid on;
title('Part I - Verification (Peak Energy Base)');
xlabel('Peak E_b/N_0 (dB)'); ylabel('Bit Error Rate (BER)');
legend('OOK (manual)', 'PRK (manual)', 'BFSK (manual)', ...
       'OOK (theory)', 'PRK (theory)', 'BFSK (theory)', ...
       'PRK (pskmod sim)');
ylim([1e-6 1]);