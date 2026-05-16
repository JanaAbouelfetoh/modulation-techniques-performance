%% Part III: Compare the Performance of Rectangular BPSK, QPSK, and 4QAM
clc; clear; close all;

%%  Simulation Parameters
N_bits      = 1e6;                   % Number of bits per SNR value
SNR_dB      = 0 : 3 : 60;            % Average Eb/N0 range in dB (3 dB steps)
SNR_linear  = 10 .^ (SNR_dB / 10);   % Convert SNR to linear scale

% Pre-allocate BER vectors
BER_BPSK = zeros(1, length(SNR_dB));
BER_QPSK = zeros(1, length(SNR_dB));
BER_4QAM = zeros(1, length(SNR_dB));
BER_4ASK = zeros(1, length(SNR_dB));   

%% =========================================================================
%  SECTION 1: BPSK
% =========================================================================
fprintf('Simulating BPSK...\n');
bits_bpsk = randi([0 1], 1, N_bits);
tx_bpsk = 2 * bits_bpsk - 1;       
P_bpsk = mean(abs(tx_bpsk).^2);    % Eb = 1

for k = 1 : length(SNR_dB)
    Eb = P_bpsk; 
    N0 = Eb / SNR_linear(k);
    noise = sqrt(N0 / 2) * randn(1, N_bits);   % Corrected: N0/2 variance
    
    rx_bpsk = tx_bpsk + noise;
    detected_bpsk = (rx_bpsk >= 0);   
    num_errors    = sum(detected_bpsk ~= bits_bpsk);
    BER_BPSK(k)   = num_errors / N_bits;
end

%% =========================================================================
%  SECTION 2: QPSK (Rectangular / Gray-coded)
% =========================================================================
fprintf('Simulating QPSK...\n');
bits_qpsk = randi([0 1], 1, N_bits);
num_sym_qpsk = N_bits / 2;
b_mat   = reshape(bits_qpsk, 2, num_sym_qpsk);   

I_tx = 2 * b_mat(1, :) - 1;    
Q_tx = 2 * b_mat(2, :) - 1;    
P_qpsk = mean(I_tx.^2 + Q_tx.^2); % Es = 2

for k = 1 : length(SNR_dB)
    Eb = P_qpsk / 2; % Eb = 1
    N0 = Eb / SNR_linear(k);
    
    % Corrected: Divided by 2 for baseband per-dimension noise variance
    noise_I = sqrt(N0 / 2) * randn(1, num_sym_qpsk);
    noise_Q = sqrt(N0 / 2) * randn(1, num_sym_qpsk);
    
    rx_I = I_tx + noise_I;
    rx_Q = Q_tx + noise_Q;
    
    det_b1 = (rx_I >= 0);   
    det_b2 = (rx_Q >= 0);   
    det_bits_qpsk = reshape([det_b1; det_b2], 1, N_bits);
    
    num_errors  = sum(det_bits_qpsk ~= bits_qpsk);
    BER_QPSK(k) = num_errors / N_bits;
end

%% =========================================================================
%  SECTION 3: 4QAM (Square QAM, equivalent to QPSK for M=4)
% =========================================================================
fprintf('Simulating 4QAM...\n');
bits_4qam = randi([0 1], 1, N_bits);
num_sym_4qam = N_bits / 2;
b_mat_qam    = reshape(bits_4qam, 2, num_sym_4qam);

I_tx_qam = 2 * b_mat_qam(1, :) - 1;
Q_tx_qam = 2 * b_mat_qam(2, :) - 1;
P_4qam = mean(I_tx_qam.^2 + Q_tx_qam.^2);   % Es = 2

for k = 1 : length(SNR_dB)
    Eb = P_4qam / 2; % Eb = 1
    N0 = Eb / SNR_linear(k);
    
    % Corrected: Divided by 2 for baseband per-dimension noise variance
    noise_I_qam = sqrt(N0 / 2) * randn(1, num_sym_4qam);
    noise_Q_qam = sqrt(N0 / 2) * randn(1, num_sym_4qam);
    
    rx_I_qam = I_tx_qam + noise_I_qam;
    rx_Q_qam = Q_tx_qam + noise_Q_qam;
    
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

gray_idx  = b_mat_ask(1, :) * 2 + b_mat_ask(2, :);   
nat_idx   = bitxor(gray_idx, bitshift(gray_idx, -1));  
tx_4ask   = 2 * nat_idx - 3;                           

P_4ask = mean(tx_4ask .^ 2); % Es = 5

for k = 1 : length(SNR_dB)
    Eb_4ask = P_4ask / 2; % Eb = 2.5
    N0 = Eb_4ask / SNR_linear(k);
    
    % Divided by 2 for baseband per-dimension noise variance
    noise_ask  = sqrt(N0 / 2) * randn(1, num_sym_4ask);
    
    rx_4ask = tx_4ask + noise_ask;
    
    levels    = [-3, -1, 1, 3];
    [~, idx]  = min(abs(rx_4ask' - levels), [], 2);  
    det_nat   = idx - 1;                              
    gray_det  = bitxor(det_nat', bitshift(det_nat', -1));
    det_b1    = floor(gray_det / 2);
    det_b2    = mod(gray_det, 2);
    det_bits_4ask = reshape([det_b1; det_b2], 1, N_bits);
    
    num_errors   = sum(det_bits_4ask ~= bits_4ask);
    BER_4ASK(k)  = num_errors / N_bits;
end

%% =========================================================================
%  PLOT: BER vs Eb/N0
% =========================================================================
fprintf('Plotting results...\n');
figure('Name', 'Part III: BER Comparison', 'NumberTitle', 'off', 'Position', [100 100 900 600]);

semilogy(SNR_dB, BER_BPSK, 'b-o',  'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BPSK');   hold on;
% Used dashed lines for QPSK and 4QAM so you can see them perfectly overlapping BPSK
semilogy(SNR_dB, BER_QPSK, 'r--s',  'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'QPSK');
semilogy(SNR_dB, BER_4QAM, 'g:^',  'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '4QAM');
semilogy(SNR_dB, BER_4ASK, 'm-d',  'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', '4ASK (Part II)');

yline(1e-5, 'k--', '10^{-5} floor', 'LabelHorizontalAlignment', 'left', 'LineWidth', 1);
grid on;
xlabel('Average E_b / N_0 (dB)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)',  'FontSize', 13, 'FontWeight', 'bold');
title('Part III: BER Performance of BPSK, QPSK, 4QAM, and 4ASK', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
xlim([0 20]); 
ylim([1e-6 1]);