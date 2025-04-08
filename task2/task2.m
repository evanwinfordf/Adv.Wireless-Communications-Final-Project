clc; clear; close all;

SNR_dB = -4:2:14;
SNR_linear = 10.^(SNR_dB/10);
N = 1e5;  % Number of bits
Es = 1;  % BPSK symbol energy

ber_awgn_sim = zeros(size(SNR_dB));
ber_rayleigh_sim = zeros(size(SNR_dB));

for i = 1:length(SNR_dB)
    sigma2 = Es ./ SNR_linear(i);
    
    % --- Generate BPSK symbols ---
    bits = randi([0 1], 1, N);
    symbols = 2*bits - 1;  % Map 0 -> -1, 1 -> +1

    % === AWGN Channel ===
    noise_awgn = sqrt(sigma2/2) * (randn(1, N));
    y_awgn = symbols + noise_awgn;
    bits_rx_awgn = real(y_awgn) > 0;
    ber_awgn_sim(i) = sum(bits_rx_awgn ~= bits) / N;

    % === Rayleigh Fading Channel (coherent detection) ===
    h = (randn(1, N) + 1j*randn(1, N)) / sqrt(2);  % Rayleigh
    noise_rayleigh = sqrt(sigma2/2) * (randn(1, N) + 1j*randn(1, N));
    y_rayleigh = h .* symbols + noise_rayleigh;

    % Coherent detection (assume h is known)
    y_eq = y_rayleigh ./ h;
    bits_rx_rayleigh = real(y_eq) > 0;
    ber_rayleigh_sim(i) = sum(bits_rx_rayleigh ~= bits) / N;
end

% === Theoretical BER ===
ber_awgn_theory = qfunc(sqrt(2*SNR_linear));
ber_rayleigh_theory = 0.5 * (1 - sqrt(SNR_linear ./ (2 + SNR_linear)));

% === Plotting ===
semilogy(SNR_dB, ber_awgn_theory, 'b-', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_awgn_sim, 'bo', 'LineWidth', 1.5);
semilogy(SNR_dB, ber_rayleigh_theory, 'r-', 'LineWidth', 2);
semilogy(SNR_dB, ber_rayleigh_sim, 'ro', 'LineWidth', 1.5);
legend('AWGN Theory', 'AWGN Sim', 'Rayleigh Theory', 'Rayleigh Sim', 'Location', 'southwest');
xlabel('SNR (dB)'); ylabel('BER'); title('BPSK BER under AWGN and Rayleigh Channels');
grid on;