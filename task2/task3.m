clc; clear; close all;

SNR_dB = -4:2:10;
SNR_linear = 10.^(SNR_dB/10);
N = 1e5;       % number of symbols
Es = 1;        % symbol energy
M = [1, 2, 4, 16]; % Number of Rx antennas

ber_mrc = zeros(length(M), length(SNR_dB));
ber_sc = zeros(length(M), length(SNR_dB));

for m = 1:length(M)
    Rx = M(m);
    for i = 1:length(SNR_dB)
        sigma2 = Es / SNR_linear(i);
        bits = randi([0 1], 1, N);
        symbols = 2*bits - 1;  % BPSK: 0 -> -1, 1 -> +1

        % Channel and noise for all Rx antennas
        h = (randn(Rx, N) + 1j*randn(Rx, N)) / sqrt(2);  % Rayleigh
        noise = sqrt(sigma2/2) * (randn(Rx, N) + 1j*randn(Rx, N));
        y = h .* repmat(symbols, Rx, 1) + noise;

        %% --- MRC combining ---
        y_mrc = sum(conj(h) .* y, 1);  % weighted sum
        bits_rx_mrc = real(y_mrc) > 0;
        ber_mrc(m,i) = sum(bits_rx_mrc ~= bits) / N;

        %% --- SC combining ---
        [~, idx] = max(abs(h), [], 1);  % best branch
        linear_idx = sub2ind(size(h), idx, 1:N);  % convert to linear indices
        y_sc = y(linear_idx);  % received signal from best antenna
        h_sc = h(linear_idx);  % channel gain from best antenna
        y_eq = y_sc ./ h_sc;   % equalize
        bits_rx_sc = real(y_eq) > 0;
        ber_sc(m,i) = sum(bits_rx_sc ~= bits) / N;
    end
end

%% === Plotting ===
figure;
markers = {'o', 's', 'd', '^'};
for m = 1:length(M)
    semilogy(SNR_dB, ber_mrc(m,:), ['-b' markers{m}], 'LineWidth', 1.5); hold on;
    semilogy(SNR_dB, ber_sc(m,:), ['--r' markers{m}], 'LineWidth', 1.5);
end
legend('MRC 1Rx','SC 1Rx','MRC 2Rx','SC 2Rx','MRC 4Rx','SC 4Rx','MRC 16Rx','SC 16Rx', 'Location', 'southwest');
xlabel('SNR (dB)'); ylabel('BER'); title('BER vs SNR for MRC and SC Beamforming');
grid on;