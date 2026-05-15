% Digital Communications Final Project
% Part I: Performance of OOK, PRK, and BFSK

clear; clc; close all;

%% 1. Simulation Parameters
N_bits = 1e6; % 1 million bits
EbN0_dB = 0:3:60; % Eb/N0 range 

% Initialize empty BER vectors
BER_OOK = []; BER_PRK = []; BER_BFSK = [];

%% 2. Generate Random Binary Data
data = randi([0 1], 1, N_bits);

%% Manual Implementation Loop (Baseband Eb/N0)
for ebno = EbN0_dB
    % Calculate linear Eb/N0 and Noise Spectral Density (N0)
    % Assuming Eb = 1 for all schemes to normalize the curves
    EbN0_lin = 10^(ebno/10);
    N0 = 1 / EbN0_lin; 
    noise_variance = N0 / 2; % Baseband noise is N0/2 per dimension
    
    %% --- PRK (BPSK) ---
    % Distance = 2. 
    s_prk = 2*data - 1; % -1 and +1
    noise_prk = sqrt(noise_variance) * randn(1, N_bits); 
    r_prk = s_prk + noise_prk;
    detected_prk = r_prk > 0;
    BER_PRK = [BER_PRK, sum(data ~= detected_prk) / N_bits];
    
    %% --- BFSK (Orthogonal) ---
    % Distance = sqrt(2). 
    s_bfsk = zeros(1, N_bits);
    s_bfsk(data == 0) = 1;
    s_bfsk(data == 1) = 1i; 
    % Complex noise (N0/2 on real, N0/2 on imag)
    noise_bfsk = sqrt(noise_variance) * (randn(1, N_bits) + 1i*randn(1, N_bits));
    r_bfsk = s_bfsk + noise_bfsk;
    detected_bfsk = imag(r_bfsk) > real(r_bfsk);
    BER_BFSK = [BER_BFSK, sum(data ~= detected_bfsk) / N_bits];

    %% --- OOK (On-Off Keying) ---
    % To match Eb/N0, if average Eb = 1, peak energy must be 2.
    % So amplitude for '1' is sqrt(2). Distance = sqrt(2).
    s_ook = sqrt(2) * data; 
    noise_ook = sqrt(noise_variance) * randn(1, N_bits);
    r_ook = s_ook + noise_ook;
    detected_ook = r_ook > (sqrt(2)/2); % Threshold is halfway
    BER_OOK = [BER_OOK, sum(data ~= detected_ook) / N_bits];
end

%% Plot Manual
figure('Name', 'Manual Implementation Performance', 'Position', [100, 100, 800, 600]);
semilogy(EbN0_dB, BER_OOK, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(EbN0_dB, BER_BFSK, 'b--s', 'LineWidth', 2, 'MarkerSize', 8); % Dashed to show overlap
semilogy(EbN0_dB, BER_PRK, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Probability of Error vs E_b/N_0 (Manual Implementation)');
xlabel('E_b/N_0 (dB)'); ylabel('Bit Error Rate (BER)');
legend('OOK', 'BFSK', 'PRK (BPSK)');
ylim([1e-5 1]);

%% 9. Built-in MATLAB Functions Evaluation
BER_OOK_builtin = zeros(1, length(EbN0_dB));
BER_PRK_builtin = zeros(1, length(EbN0_dB));
BER_BFSK_builtin = zeros(1, length(EbN0_dB));

M = 2; freq_sep = 1; nsamp = 8; Fs = 8; 

for i = 1:length(EbN0_dB)
    ebno = EbN0_dB(i);
    disp(['Processing Built-in Eb/N0: ', num2str(ebno), ' dB...']); 
    
    % --- PRK Built-in ---
    s_prk_built = pskmod(data, 2); 
    % awgn defaults to mapping SNR to Eb/N0 internally if we specify power
    r_prk_built = awgn(s_prk_built, ebno, 'measured');
    det_prk_built = pskdemod(r_prk_built, 2);
    [~, BER_PRK_builtin(i)] = biterr(data, det_prk_built);
    
    % --- BFSK Built-in ---
    s_bfsk_built = fskmod(data, M, freq_sep, nsamp, Fs); 
    % Adjust for oversampling so Eb/N0 remains accurate
    r_bfsk_built = awgn(s_bfsk_built, ebno - 10*log10(nsamp), 'measured');
    det_bfsk_built = fskdemod(r_bfsk_built, M, freq_sep, nsamp, Fs);
    [~, BER_BFSK_builtin(i)] = biterr(data, det_bfsk_built);

    % --- OOK Built-in ---
    % Standard PAM gives +1/-1. We shift to 0 and 1.
    s_pam = pammod(data, 2);
    s_ook_built = (s_pam + 1) / 2; 
    % OOK average power is 3dB lower than peak. To plot against Eb/N0, 
    % we must trick awgn by artificially adding 3dB to the noise target.
    r_ook_built = awgn(s_ook_built, ebno + 3, 'measured');
    det_ook_built = pamdemod(2 * r_ook_built - 1, 2);
    [~, BER_OOK_builtin(i)] = biterr(data, det_ook_built);
end

%% Plot Built-in
figure('Name', 'Built-in Functions Performance', 'Position', [150, 150, 800, 600]);
semilogy(EbN0_dB, BER_OOK_builtin, 'r-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(EbN0_dB, BER_BFSK_builtin, 'b--s', 'LineWidth', 2, 'MarkerSize', 8);
semilogy(EbN0_dB, BER_PRK_builtin, 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
title('Probability of Error vs E_b/N_0 (MATLAB Built-in Verification)');
xlabel('E_b/N_0 (dB)'); ylabel('Bit Error Rate (BER)');
legend('OOK (Built-in)', 'BFSK (Built-in)', 'PRK (Built-in)');
ylim([1e-5 1]);