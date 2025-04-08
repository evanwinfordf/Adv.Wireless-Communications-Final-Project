clc; clear; close all;

% Parameters
SNR_dB = -4:2:10;
SNR_linear = 10.^(SNR_dB/10);
N = 1e4;            % Reduced for performance (was 1e5)
Es = 1;             % Symbol energy
Nt = 8; Nr = 4;     % MIMO dimensions (8x4)

ber_mimo = zeros(size(SNR_dB));
snr_eff = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    sigma2 = Es / SNR_linear(i);  % Noise variance
    errors = 0;
    snr_tmp = 0;

    % --- Channel Matrix ---
    H = (randn(Nr, Nt) + 1j*randn(Nr, Nt)) / sqrt(2);  % Rayleigh flat fading

    % --- Beamforming ---
    [~, ~, V] = svd(H);       % SVD of H
    w_tx = V(:,1);            % MRT beamformer
    h_eff = H * w_tx;         % Effective channel vector
    w_rx = h_eff / norm(h_eff);  % MRC receiver

    h_total = w_rx' * H * w_tx;  % Scalar channel gain after beamforming

    % --- Generate BPSK symbols ---
    bits = randi([0 1], 1, N);
    symbols = 2 * bits - 1;  % BPSK: 0 -> -1, 1 -> +1

    % --- Transmit through channel ---
    noise = sqrt(sigma2/2) * (randn(Nr, N) + 1j*randn(Nr, N));
    y = H * (w_tx * symbols) + noise;  % Matrix mult with broadcasting
    y_combined = w_rx' * y;  % Combine via receiver beamforming

    % --- Detection ---
    bits_rx = real(y_combined) > 0;
    errors = sum(bits_rx ~= bits);

    % --- Effective SNR ---
    snr_tmp = abs(h_total)^2 * Es / sigma2;

    % --- Store results ---
    ber_mimo(i) = errors / N;
    snr_eff(i) = snr_tmp;
end

% === Plot BER ===
figure;
semilogy(SNR_dB, ber_mimo, 'b-o', 'LineWidth', 2); grid on;
xlabel('SNR (dB)'); ylabel('BER');
title('BER Performance of 8x4 MIMO Beamforming');
legend('8x4 MIMO Beamforming', 'Location', 'southwest');

% === Plot Effective SNR ===
figure;
plot(SNR_dB, 10*log10(snr_eff), 'r-*', 'LineWidth', 2); grid on;
xlabel('Input SNR (dB)');
ylabel('Effective SNR after Beamforming (dB)');
title('Effective SNR of 8x4 MIMO System');
legend('Effective SNR');