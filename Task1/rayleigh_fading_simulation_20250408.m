
clc; clear; close all;

% --- Parameters ---
fc = 3.5e9;                 % Carrier frequency (Hz)
c = 3e8;                    % Speed of light (m/s)
lambda = c / fc;           % Wavelength
Pt_dBm = 43;               % Transmit power in dBm
Pt = 10^((Pt_dBm - 30)/10); % Transmit power in watts
noise_dBm = -100;
sigma_noise = 10^((noise_dBm - 30)/10);
hb = 25;                   % Base station height

% User distances
d1 = 100;  % Near user (m)
d2 = 900;  % Far user (m)

% Shadowing
shadow_std = 6;

% Compute Path Loss (3GPP NLOS)
fc_GHz = fc / 1e9;
PL1 = 13.54 + 39.08*log10(d1) + 20*log10(fc_GHz) - 0.6*(hb - 1.5);
PL2 = 13.54 + 39.08*log10(d2) + 20*log10(fc_GHz) - 0.6*(hb - 1.5);
shadowing1 = normrnd(0, shadow_std);
shadowing2 = normrnd(0, shadow_std);
PL_total_dB1 = PL1 + shadowing1;
PL_total_dB2 = PL2 + shadowing2;
g1 = 10^(-PL_total_dB1 / 10);
g2 = 10^(-PL_total_dB2 / 10);

% Simulation setup
fs = 1000;
t = 0:1/fs:2;
N = 20;  % multipath
v = [3, 30, 100] / 3.6;  % speeds in m/s

for vi = 1:length(v)
    speed = v(vi);
    fd = speed / lambda;
    phi = 2*pi*rand(1,N);
    theta = 2*pi*(1:N)/N;
    beta = 1/sqrt(N);

    Zt = zeros(size(t));
    for n = 1:N
        Zt = Zt + beta * exp(1j*(2*pi*fd*cos(theta(n))*t + phi(n)));
    end

    h1 = sqrt(g1) * Zt;
    h2 = sqrt(g2) * Zt;
    Pr1 = Pt * abs(h1).^2;
    Pr2 = Pt * abs(h2).^2;
    SNR1 = Pr1 / sigma_noise;
    SNR2 = Pr2 / sigma_noise;

    % Plot h(t)
    figure;
    subplot(2,1,1); plot(t, 20*log10(abs(h1)));
    title(['|h(t)| in dB for User 1 @ ', num2str(speed*3.6), ' km/h']);
    ylabel('Magnitude (dB)'); xlabel('Time (s)');
    subplot(2,1,2); plot(t, 20*log10(abs(h2)));
    title(['|h(t)| in dB for User 2 @ ', num2str(speed*3.6), ' km/h']);
    ylabel('Magnitude (dB)'); xlabel('Time (s)');

    % Plot PDF of SNR
    figure;
    histogram(10*log10(SNR1), 'Normalization', 'pdf', 'FaceAlpha', 0.6); hold on;
    histogram(10*log10(SNR2), 'Normalization', 'pdf', 'FaceAlpha', 0.6);
    title(['Empirical SNR PDF @ ', num2str(speed*3.6), ' km/h']);
    legend('User 1', 'User 2');
    xlabel('SNR (dB)'); ylabel('PDF'); grid on;

    % Optional: Overlay theoretical PDF
    figure;
    snr_vals = linspace(0, max([SNR1 SNR2]), 1000);
    bar_gamma1 = mean(SNR1);
    bar_gamma2 = mean(SNR2);
    pdf1 = (1/bar_gamma1) * exp(-snr_vals/bar_gamma1);
    pdf2 = (1/bar_gamma2) * exp(-snr_vals/bar_gamma2);
    plot(snr_vals, pdf1, 'b', 'LineWidth', 2); hold on;
    plot(snr_vals, pdf2, 'r', 'LineWidth', 2);
    title(['Theoretical SNR PDF @ ', num2str(speed*3.6), ' km/h']);
    legend('User 1', 'User 2');
    xlabel('SNR'); ylabel('PDF'); grid on;
end
