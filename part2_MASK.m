%  Part II: Compare the Performance of M-ASK for M = 2, 4, 8


clc; clear; close all;

%% =========================================================
%  SIMULATION PARAMETERS
%% =========================================================
num_bits      = 1e6;               % Number of bits per SNR value
SNR_dB        = 0:3:60;           % SNR range in dB (0 to 60, step 3)
M_values      = [2, 4, 8];        % ASK orders to simulate
num_M         = length(M_values);

% Pre-allocate BER storage: rows = M-ASK order, cols = SNR points
BER_sim = zeros(num_M, length(SNR_dB));
BER_theory = zeros(num_M, length(SNR_dB));

%% =========================================================
%  IMPORTANT: Minimum Distance Constraint
%  For M-ASK with equally spaced levels and minimum distance d_min,
%  symbol levels are chosen as:  -(M-1), -(M-3), ..., (M-3), (M-1)
%  scaled so that d_min = 2 for all M (i.e., adjacent levels differ by 2).
%  This ensures a FAIR comparison across M values.
%% =========================================================

for m_idx = 1:num_M
    M      = M_values(m_idx);          % Current ASK order (2, 4, or 8)
    k      = log2(M);                  % Bits per symbol
    
    % --- Step 3: Generate M-ASK Constellation ---
    % Levels: -(M-1), -(M-3), ..., (M-3), (M-1)
    % Minimum distance = 2 (same for all M)
    levels = -(M-1) : 2 : (M-1);      % e.g., M=2 -> [-1,1]; M=4 -> [-3,-1,1,3]
    
    % Average symbol power for these levels
    % E_s = mean(levels.^2)
    E_s_avg = mean(levels .^ 2);
    
    % --- Generate random bit stream ---
    bits = randi([0 1], 1, num_bits);
    
    % --- Group bits into symbols ---
    % Zero-pad if necessary so total bits is divisible by k
    num_pad = mod(-num_bits, k);           % padding needed
    bits_padded = [bits, zeros(1, num_pad)];
    num_symbols = length(bits_padded) / k;
    
    % Reshape into (k x num_symbols) and convert each column from binary to decimal
    % No bi2de (requires Communications Toolbox) — use matrix multiply instead:
    %   weights [2^(k-1), 2^(k-2), ..., 1] dot each row of bit_matrix gives decimal value
    bit_matrix   = reshape(bits_padded, k, num_symbols);  % k rows, num_symbols cols
    weights      = 2.^(k-1 : -1 : 0);                    % MSB-first weights, 1 x k
    symbol_index = weights * bit_matrix;                  % 1 x num_symbols, decimal 0..M-1
    
    % Map symbol indices to ASK amplitude levels
    tx_symbols = levels(symbol_index + 1);   % MATLAB indexing: +1
    
    fprintf('Simulating %d-ASK (M=%d, k=%d bits/symbol, E_s=%.4f)...\n', ...
             M, M, k, E_s_avg);
    
    %% -------------------------------------------------------
    %  Loop over each SNR value
    %% -------------------------------------------------------
    BER = [];   % Reset BER vector for this M
    
    for snr_idx = 1:length(SNR_dB)
        
        snr_linear = 10^(SNR_dB(snr_idx)/10);   % Convert SNR from dB to linear
        
        % --- Step 4: Add AWGN Noise ---
        % SNR (per symbol) = E_s / N0
        % => N0 = E_s / SNR
        % Noise variance (one-sided) sigma^2 = N0/2
        N0        = E_s_avg / snr_linear;
        sigma     = sqrt(N0 / 2);               % Std deviation of noise
        noise     = sigma * randn(1, num_symbols);  % Real AWGN (ASK is 1D)
        
        rx_signal = tx_symbols + noise;         % Received signal = sent + noise
        
        % --- Step 5: Decision / Detection ---
        % Minimum distance detector: find the closest ASK level for each sample
        % Expand for vectorised nearest-neighbour search
        % rx_signal: 1 x num_symbols
        % levels:    1 x M
        % Compute distance matrix: num_symbols x M
        dist_matrix   = abs(rx_signal' - levels);        % num_symbols x M
        [~, detected_idx] = min(dist_matrix, [], 2);     % index of closest level
        detected_idx  = detected_idx' - 1;               % back to 0-based indices
        
        % --- Decode detected indices back to bits ---
        % No de2bi (requires Communications Toolbox) — use bitshift/bitand instead:
        %   extract each bit by shifting right and masking with 1
        detected_bits_matrix = zeros(num_symbols, k);
        for bit_pos = 1:k
            detected_bits_matrix(:, bit_pos) = bitand(floor(detected_idx / 2^(k - bit_pos)), 1);
        end
        detected_bits = reshape(detected_bits_matrix', 1, []);       % 1 x (num_symbols*k)
        
        % Trim to original bit length (remove padding)
        detected_bits = detected_bits(1:num_bits);
        
        % --- Step 6: Count bit errors ---
        num_errors = sum(bits ~= detected_bits);
        
        % --- Step 7: Calculate and store BER ---
        prob_error = num_errors / num_bits;
        BER = [BER, prob_error];                          %#ok<AGROW>
    end
    
    BER_sim(m_idx, :) = BER;
    
    %% -------------------------------------------------------
    %  Theoretical BER for M-ASK (Gray-coded, AWGN)
    %  BER ≈ (2*(M-1)/M) / log2(M)  *  Q( sqrt(6*log2(M)*SNR_b / (M^2-1)) )
    %  where SNR_b = Eb/N0 (per bit SNR)
    %  Note: SNR_dB here refers to Es/N0; Eb/N0 = (Es/N0) / log2(M)
    %% -------------------------------------------------------
    for snr_idx = 1:length(SNR_dB)
        snr_linear = 10^(SNR_dB(snr_idx)/10);   % Es/N0 linear
        Eb_N0      = snr_linear / k;              % Eb/N0 linear
        
        % Theoretical BER for M-ASK with Gray coding
        arg = sqrt(6 * k * Eb_N0 / (M^2 - 1));
        BER_theory(m_idx, snr_idx) = (2*(M-1)/M) / k * qfunc(arg);
    end
    
end

%% =========================================================
%  Step 8: Plot BER Curves
%% =========================================================
figure('Name','Part II: M-ASK BER Performance','NumberTitle','off', ...
       'Color','white','Position',[100 100 900 600]);

colors  = {'b', 'r', 'g'};
markers = {'o', 's', '^'};
M_labels = {'2-ASK (BASK/OOK)','4-ASK','8-ASK'};

hold on; grid on;

for m_idx = 1:num_M
    % Simulated BER
    semilogy(SNR_dB, BER_sim(m_idx,:), ...
             [colors{m_idx}, markers{m_idx}, '-'], ...
             'LineWidth', 1.5, 'MarkerSize', 6, ...
             'DisplayName', [M_labels{m_idx}, ' (Simulated)']);
    
    % Theoretical BER
    semilogy(SNR_dB, BER_theory(m_idx,:), ...
             [colors{m_idx}, '--'], ...
             'LineWidth', 1.2, ...
             'DisplayName', [M_labels{m_idx}, ' (Theoretical)']);
end

xlabel('SNR  (E_s/N_0)  [dB]', 'FontSize', 12);
ylabel('Bit Error Rate (BER)', 'FontSize', 12);
title('Part II: BER Performance of M-ASK (M = 2, 4, 8)', 'FontSize', 14);
legend('Location','southwest','FontSize',10);
xlim([0 60]);
ylim([1e-6 1]);
hold off;

fprintf('\nSimulation complete. Figure displayed.\n');

%% =========================================================
%  HELPER FUNCTION: Q-function
%% =========================================================
function y = qfunc(x)
    % Q(x) = 0.5 * erfc(x / sqrt(2))
    y = 0.5 * erfc(x / sqrt(2));
end
