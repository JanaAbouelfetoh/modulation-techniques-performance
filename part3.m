%% Part III: Compare the Performance of Rectangular BPSK, QPSK, and 4QAM
clc; clear; close all;

%%  Simulation Parameters

N_bits      = 1e6;                    % Number of bits per SNR value
SNR_dB      = 0 : 3 : 60;            % SNR range in dB (3 dB steps)
SNR_linear  = 10 .^ (SNR_dB / 10);   % Convert SNR to linear scale

% Pre-allocate BER vectors
BER_BPSK = zeros(1, length(SNR_dB));
BER_QPSK = zeros(1, length(SNR_dB));
BER_4QAM = zeros(1, length(SNR_dB));
BER_4ASK = zeros(1, length(SNR_dB));   % Re-simulated here for the plot

%% =========================================================================
%  SECTION 1: BPSK
%  1 bit per symbol. Constellation: {-1 (bit 0), +1 (bit 1)}
%  Minimum distance d_min = 2  (between -1 and +1)
% =========================================================================
fprintf('Simulating BPSK...\n');

% Generate random bits
bits_bpsk = randi([0 1], 1, N_bits);

% Modulate: map 0 -> -1, 1 -> +1
tx_bpsk = 2 * bits_bpsk - 1;       % BPSK symbols (real-valued)

% Signal power (should be 1 for BPSK with +/-1)
P_bpsk = mean(abs(tx_bpsk).^2);

for k = 1 : length(SNR_dB)

    % Noise variance: sigma^2 = P / SNR_linear
    % For BPSK the SNR here is Eb/N0 (1 bit per symbol)
    sigma2 = P_bpsk / SNR_linear(k);
    noise  = sqrt(sigma2 /2) * randn(1, N_bits);   % Real AWGN

    % Received signal
    rx_bpsk = tx_bpsk + noise;

    % Decision: threshold at 0
    detected_bpsk = (rx_bpsk >= 0);   % 1 if >= 0, 0 otherwise

    % Count errors and compute BER
    num_errors    = sum(detected_bpsk ~= bits_bpsk);
    BER_BPSK(k)   = num_errors / N_bits;
end

%% =========================================================================
%  SECTION 2: QPSK (Rectangular / Gray-coded)
% =========================================================================
fprintf('Simulating QPSK...\n');

% Generate random bits (must be even number for 2 bits/symbol)
bits_qpsk = randi([0 1], 1, N_bits);

% Group into pairs; each pair selects one QPSK symbol
% Pair [b1 b2]: b1 -> I component, b2 -> Q component
% Mapping: 0 -> -1, 1 -> +1
num_sym_qpsk = N_bits / 2;
b_mat   = reshape(bits_qpsk, 2, num_sym_qpsk);   % 2 x num_sym matrix

% I (in-phase) and Q (quadrature) components
I_tx = 2 * b_mat(1, :) - 1;    % first bit of each pair
Q_tx = 2 * b_mat(2, :) - 1;    % second bit of each pair

% Signal power per symbol: E[I^2 + Q^2] = 1 + 1 = 2
P_qpsk = mean(I_tx.^2 + Q_tx.^2);

for k = 1 : length(SNR_dB)

    % Eb/N0 given; Es = 2*Eb => sigma^2 = Es / (2 * SNR_linear) per dimension
    % sigma^2 per dimension = (Es/2) / SNR_linear = Eb / SNR_linear = P_bpsk / SNR_linear
    % Use the same noise variance as BPSK for fair per-bit comparison
    sigma2_per_dim = (P_qpsk / 2) / SNR_linear(k);   % noise var per I or Q channel
    noise_I = sqrt(sigma2_per_dim) * randn(1, num_sym_qpsk);
    noise_Q = sqrt(sigma2_per_dim) * randn(1, num_sym_qpsk);

    % Received I and Q
    rx_I = I_tx + noise_I;
    rx_Q = Q_tx + noise_Q;

    % Decision: threshold at 0 on each axis independently
    det_b1 = (rx_I >= 0);   % detected first bit
    det_b2 = (rx_Q >= 0);   % detected second bit

    % Reconstruct detected bit stream
    det_bits_qpsk = reshape([det_b1; det_b2], 1, N_bits);

    % BER (per bit)
    num_errors  = sum(det_bits_qpsk ~= bits_qpsk);
    BER_QPSK(k) = num_errors / N_bits;
end

%% =========================================================================
%  SECTION 3: 4QAM (Square QAM, equivalent to QPSK for M=4)
%  4QAM and QPSK have identical BER when the same Gray coding and d_min
% =========================================================================
fprintf('Simulating 4QAM...\n');

bits_4qam = randi([0 1], 1, N_bits);

num_sym_4qam = N_bits / 2;
b_mat_qam    = reshape(bits_4qam, 2, num_sym_4qam);

% Gray code mapping for I axis:  0 -> -1, 1 -> +1
% Gray code mapping for Q axis:  0 -> -1, 1 -> +1
I_tx_qam = 2 * b_mat_qam(1, :) - 1;
Q_tx_qam = 2 * b_mat_qam(2, :) - 1;

P_4qam = mean(I_tx_qam.^2 + Q_tx_qam.^2);   % = 2

for k = 1 : length(SNR_dB)

    sigma2_per_dim_qam = (P_4qam / 2) / SNR_linear(k);
    noise_I_qam = sqrt(sigma2_per_dim_qam) * randn(1, num_sym_4qam);
    noise_Q_qam = sqrt(sigma2_per_dim_qam) * randn(1, num_sym_4qam);

    rx_I_qam = I_tx_qam + noise_I_qam;
    rx_Q_qam = Q_tx_qam + noise_Q_qam;

    % Nearest-neighbour decision (threshold at 0 for 2-level per axis)
    det_b1_qam = (rx_I_qam >= 0);
    det_b2_qam = (rx_Q_qam >= 0);

    det_bits_4qam = reshape([det_b1_qam; det_b2_qam], 1, N_bits);

    num_errors   = sum(det_bits_4qam ~= bits_4qam);
    BER_4QAM(k)  = num_errors / N_bits;
end

%% =========================================================================
%  SECTION 4: 4-ASK
% =========================================================================
fprintf('Simulating 4-ASK...\n');

bits_4ask = randi([0 1], 1, N_bits);

num_sym_4ask = N_bits / 2;
b_mat_ask    = reshape(bits_4ask, 2, num_sym_4ask);

% Convert each 2-bit group to a symbol index (Gray -> natural binary -> symbol)
% Gray code: 00->0, 01->1, 11->2, 10->3
% Symbol levels: index 0->-3, 1->-1, 2->+1, 3->+3
gray_idx  = b_mat_ask(1, :) * 2 + b_mat_ask(2, :);   % gray code value 0..3
nat_idx   = bitxor(gray_idx, bitshift(gray_idx, -1));  % gray to natural binary
tx_4ask   = 2 * nat_idx - 3;                           % map 0,1,2,3 -> -3,-1,+1,+3

% Signal power: E[X^2] for uniform {-3,-1,+1,+3} = (9+1+1+9)/4 = 5
P_4ask = mean(tx_4ask .^ 2);

for k = 1 : length(SNR_dB)

    % Noise for the 4-ASK real channel
    % Use per-bit SNR: Es = 2*Eb => sigma^2 = Es / (2*SNR_linear) = Eb/SNR_linear
    Eb_4ask   = P_4ask / 2;                          % energy per bit
    sigma2_ask = Eb_4ask / SNR_linear(k);
    noise_ask  = sqrt(sigma2_ask) * randn(1, num_sym_4ask);

    rx_4ask = tx_4ask + noise_ask;

    % ML decision: find nearest symbol in {-3,-1,+1,+3}
    levels    = [-3, -1, 1, 3];
    [~, idx]  = min(abs(rx_4ask' - levels), [], 2);  % index of nearest level
    det_nat   = idx - 1;                              % natural binary 0..3

    % Convert natural binary back to Gray code bits
    gray_det  = bitxor(det_nat', bitshift(det_nat', -1));
    det_b1    = floor(gray_det / 2);
    det_b2    = mod(gray_det, 2);

    det_bits_4ask = reshape([det_b1; det_b2], 1, N_bits);

    num_errors   = sum(det_bits_4ask ~= bits_4ask);
    BER_4ASK(k)  = num_errors / N_bits;
end

%% =========================================================================
%  PLOT: BER vs Eb/N0 for all four modulation schemes
% =========================================================================
fprintf('Plotting results...\n');

figure('Name', 'Part III: BER Comparison', 'NumberTitle', 'off', ...
       'Position', [100 100 900 600]);

semilogy(SNR_dB, BER_BPSK, 'b-o',  'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'BPSK');   hold on;
semilogy(SNR_dB, BER_QPSK, 'r-s',  'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'QPSK');
semilogy(SNR_dB, BER_4QAM, 'g-^',  'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', '4QAM');
semilogy(SNR_dB, BER_4ASK, 'm-d',  'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', '4ASK (Part II)');

% Reference floor line
yline(1e-5, 'k--', '10^{-5} floor', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1);

grid on;
xlabel('E_b / N_0 (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)',  'FontSize', 13, 'FontWeight', 'bold');
title('Part III: BER Performance of BPSK, QPSK, 4QAM, and 4ASK', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
xlim([0 60]);
ylim([1e-7 1]);

fprintf('\nDone. Figure displayed.\n');