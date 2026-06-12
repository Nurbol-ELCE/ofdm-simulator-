% =========================================================================
%  OFDM SIMULATOR — Built from scratch in MATLAB
%  Author : Nurbol Sissenbay
%  Date   : 2026
%
%  Pipeline:
%  Bits → S/P → QPSK Mod → IFFT → Add CP → Channel → Remove CP →
%  FFT → Equalization → QPSK Demod → P/S → BER
%
%  Features:
%   - QPSK modulation / demodulation
%   - OFDM with configurable N_fft and cyclic prefix
%   - AWGN channel
%   - Rayleigh multipath fading channel
%   - Zero-Forcing equalization (perfect CSI + pilot-based)
%   - BER vs Eb/N0 curves with theoretical overlay
%   - Constellation diagram at multiple SNR values
% =========================================================================

clc; clear; close all;

% =========================================================================
%  SECTION 1 — SYSTEM PARAMETERS
% =========================================================================

N_fft       = 64;       % FFT size (number of subcarriers)
N_cp        = 16;       % Cyclic prefix length (must be >= channel length)
N_pilots    = 8;        % Number of pilot subcarriers for channel estimation
N_symbols   = 100;      % Number of OFDM symbols to simulate per SNR point
M           = 4;        % Modulation order (4 = QPSK)
bits_per_sym = log2(M); % Bits per QPSK symbol = 2

% Subcarrier allocation
% We reserve N_pilots subcarriers for pilots, rest carry data
% Pilot indices spaced evenly across the subcarriers
pilot_idx   = round(linspace(1, N_fft, N_pilots + 2));
pilot_idx   = pilot_idx(2:end-1);          % Remove edge subcarriers
data_idx    = setdiff(1:N_fft, pilot_idx); % All remaining subcarriers
N_data      = length(data_idx);            % Number of data subcarriers

bits_per_ofdm_sym = N_data * bits_per_sym; % Bits per OFDM symbol

% SNR range for BER simulation
EbN0_dB     = 0:2:24;  % Eb/N0 in dB
EbN0_linear = 10.^(EbN0_dB / 10);

% Rayleigh channel impulse response (multipath)
% 3-tap channel: tap delays and powers (must be < N_cp)
h_channel   = [0.8, 0.5*exp(1j*pi/4), 0.3*exp(1j*pi/3)]; % complex taps
channel_len = length(h_channel);   % Must be <= N_cp — verified: 3 <= 16 ✓

fprintf('=== OFDM System Parameters ===\n');
fprintf('FFT size         : %d\n', N_fft);
fprintf('Cyclic prefix    : %d samples\n', N_cp);
fprintf('Data subcarriers : %d\n', N_data);
fprintf('Pilot subcarriers: %d\n', N_pilots);
fprintf('Bits per OFDM sym: %d\n', bits_per_ofdm_sym);
fprintf('Channel taps     : %d (CP=%d >= taps=%d ✓)\n\n', channel_len, N_cp, channel_len);

% =========================================================================
%  SECTION 2 — QPSK MODULATOR / DEMODULATOR
% =========================================================================
% QPSK constellation: Gray coded
% 00 → +1+1j,  01 → -1+1j,  11 → -1-1j,  10 → +1-1j  (normalized)

function symbols = qpsk_modulate(bits)
    % Input : bit stream (column vector, even length)
    % Output: complex QPSK symbols
    bits = bits(:);
    assert(mod(length(bits),2)==0, 'Bit stream must have even length');
    bits_paired = reshape(bits, 2, [])';  % Nx2 matrix
    b0 = bits_paired(:,1);
    b1 = bits_paired(:,2);
    % Gray coded mapping
    I = 1 - 2*b0;   % b0=0 → I=+1,  b0=1 → I=-1
    Q = 1 - 2*b1;   % b1=0 → Q=+1,  b1=1 → Q=-1
    symbols = (I + 1j*Q) / sqrt(2);  % Normalize power to 1
end

function bits = qpsk_demodulate(symbols)
    % Input : complex QPSK symbols (possibly noisy)
    % Output: bit stream using nearest-neighbor decision
    symbols = symbols(:);
    I = real(symbols);
    Q = imag(symbols);
    b0 = double(I < 0);  % I<0 → b0=1, I>0 → b0=0
    b1 = double(Q < 0);  % Q<0 → b1=1, Q>0 → b1=0
    bits_matrix = [b0, b1];
    bits = reshape(bits_matrix', [], 1);
end

% =========================================================================
%  SECTION 3 — OFDM TRANSMITTER
% =========================================================================

function tx_signal = ofdm_transmit(data_bits, pilot_idx, data_idx, ...
                                    N_fft, N_cp, N_symbols)
    % Modulate bits to QPSK symbols
    data_syms = qpsk_modulate(data_bits);   % N_data * N_symbols x 1
    data_syms = reshape(data_syms, [], N_symbols); % N_data x N_symbols

    % Known pilot symbol (BPSK pilot: value = 1+0j, known to receiver)
    pilot_val = 1 + 0j;

    tx_signal = [];
    for k = 1:N_symbols
        % Build frequency-domain OFDM symbol
        freq_sym = zeros(N_fft, 1);
        freq_sym(data_idx)  = data_syms(:, k);       % Insert data
        freq_sym(pilot_idx) = pilot_val;              % Insert pilots

        % IFFT: frequency domain → time domain
        time_sym = ifft(freq_sym, N_fft) * sqrt(N_fft); % Unitary IFFT

        % Add cyclic prefix: copy last N_cp samples to the front
        cp        = time_sym(end - N_cp + 1 : end);
        tx_signal = [tx_signal; cp; time_sym];  % Append [CP | OFDM symbol]
    end
end

% =========================================================================
%  SECTION 4 — CHANNEL MODELS
% =========================================================================

function rx_signal = awgn_channel(tx_signal, snr_linear)
    % Add complex AWGN noise
    % snr_linear is per-sample SNR (signal power / noise power)
    signal_power = mean(abs(tx_signal).^2);
    noise_var    = signal_power / (2 * snr_linear); % Factor 2 for complex noise
    noise        = sqrt(noise_var) * (randn(size(tx_signal)) + ...
                    1j * randn(size(tx_signal)));
    rx_signal    = tx_signal + noise;
end

function rx_signal = rayleigh_channel(tx_signal, h, snr_linear)
    % Convolve with multipath channel then add AWGN
    rx_faded     = conv(tx_signal, h(:));    % Linear convolution
    rx_faded     = rx_faded(1:length(tx_signal)); % Trim to original length
    rx_signal    = awgn_channel(rx_faded, snr_linear);
end

% =========================================================================
%  SECTION 5 — OFDM RECEIVER
% =========================================================================

function [rx_bits, H_est] = ofdm_receive(rx_signal, pilot_idx, data_idx, ...
                                          N_fft, N_cp, N_symbols, use_pilot_eq)
    % use_pilot_eq: true = pilot-based LS estimation, false = perfect (passed externally)

    rx_bits = [];
    H_est_all = zeros(N_fft, N_symbols);

    for k = 1:N_symbols
        % Extract kth OFDM symbol (strip CP first)
        start_idx = (k-1) * (N_fft + N_cp) + 1;
        sym_with_cp = rx_signal(start_idx : start_idx + N_fft + N_cp - 1);

        % Remove cyclic prefix
        sym_no_cp = sym_with_cp(N_cp + 1 : end);   % Keep last N_fft samples

        % FFT: time domain → frequency domain
        freq_sym = fft(sym_no_cp, N_fft) / sqrt(N_fft); % Unitary FFT

        if use_pilot_eq
            % ── PILOT-BASED LEAST SQUARES CHANNEL ESTIMATION ──
            % At pilot positions: Y[k] = H[k] * X_pilot → H[k] = Y[k] / X_pilot
            H_at_pilots = freq_sym(pilot_idx) / (1 + 0j); % pilot_val = 1

            % Interpolate H across all subcarriers (linear interpolation)
            H_full = interp1(pilot_idx', H_at_pilots, ...
                             (1:N_fft)', 'linear', 'extrap');
        else
            % Perfect CSI placeholder — filled externally
            H_full = ones(N_fft, 1); % Will be overridden in BER loop
        end

        H_est_all(:, k) = H_full;

        % ── ZERO-FORCING EQUALIZATION ──
        % Undo channel effect: X_hat[k] = Y[k] / H_est[k]
        eq_sym = freq_sym ./ H_full;

        % Extract data subcarriers only
        data_eq = eq_sym(data_idx);

        % Demodulate QPSK
        sym_bits = qpsk_demodulate(data_eq);
        rx_bits  = [rx_bits; sym_bits];
    end
    H_est = H_est_all;
end

% =========================================================================
%  SECTION 6 — BER SIMULATION LOOP
% =========================================================================

fprintf('=== Running BER Simulation ===\n');
fprintf('SNR points: %d | OFDM symbols per point: %d\n\n', ...
        length(EbN0_dB), N_symbols);

BER_awgn      = zeros(1, length(EbN0_dB));
BER_rayleigh  = zeros(1, length(EbN0_dB));
BER_pilot     = zeros(1, length(EbN0_dB));

% Store constellation data for 3 SNR values
plot_snr_idx  = [1, 5, length(EbN0_dB)];  % Low / mid / high SNR
const_symbols = cell(3, 2);  % {snr_point, channel_type}

for snr_idx = 1:length(EbN0_dB)

    % Convert Eb/N0 to per-subcarrier SNR
    % For QPSK: SNR = Eb/N0 * bits_per_sym * N_data/N_fft
    snr_linear = EbN0_linear(snr_idx) * bits_per_sym * (N_data / N_fft);

    % Also need Rayleigh channel in frequency domain for perfect EQ
    % H_freq = fft(h_channel, N_fft) — channel frequency response
    H_freq_true = fft(h_channel, N_fft).';  % N_fft x 1

    total_bits      = bits_per_ofdm_sym * N_symbols;
    errors_awgn     = 0;
    errors_rayleigh = 0;
    errors_pilot    = 0;

    % Monte Carlo: run enough trials for statistical reliability
    % Target: at least 100 errors per SNR point
    n_trials = max(1, ceil(100 / (total_bits * 0.1))); % adaptive
    n_trials = min(n_trials, 5); % cap at 5 trials for speed

    for trial = 1:n_trials

        % Generate random bits
        tx_bits = randi([0 1], bits_per_ofdm_sym * N_symbols, 1);

        % ── TRANSMIT ──
        tx_signal = ofdm_transmit(tx_bits, pilot_idx, data_idx, ...
                                   N_fft, N_cp, N_symbols);

        % ─────────────────────────────────────────────
        %  SCENARIO A: AWGN channel + perfect CSI
        % ─────────────────────────────────────────────
        rx_awgn = awgn_channel(tx_signal, snr_linear);

        rx_bits_awgn = [];
        for k = 1:N_symbols
            start_idx   = (k-1)*(N_fft+N_cp)+1;
            sym_cp      = rx_awgn(start_idx:start_idx+N_fft+N_cp-1);
            sym         = sym_cp(N_cp+1:end);
            freq_sym    = fft(sym, N_fft) / sqrt(N_fft);
            % Perfect equalization: H=1 for AWGN (flat channel)
            data_eq     = freq_sym(data_idx);
            rx_bits_awgn = [rx_bits_awgn; qpsk_demodulate(data_eq)];
        end

        errors_awgn = errors_awgn + sum(tx_bits ~= rx_bits_awgn);

        % Save constellation at selected SNR points
        if any(snr_idx == plot_snr_idx) && trial == 1
            plot_pos = find(plot_snr_idx == snr_idx);
            const_symbols{plot_pos, 1} = reshape(rx_bits_awgn, 2, [])'; % store for plot
            % Actually store frequency-domain received symbols directly:
            k_mid = round(N_symbols/2);
            start_idx = (k_mid-1)*(N_fft+N_cp)+1;
            sym_cp = rx_awgn(start_idx:start_idx+N_fft+N_cp-1);
            sym = sym_cp(N_cp+1:end);
            freq_sym = fft(sym,N_fft)/sqrt(N_fft);
            const_symbols{plot_pos,1} = freq_sym(data_idx);
        end

        % ─────────────────────────────────────────────
        %  SCENARIO B: Rayleigh fading + perfect CSI
        % ─────────────────────────────────────────────
        rx_rayleigh = rayleigh_channel(tx_signal, h_channel, snr_linear);

        rx_bits_ray = [];
        for k = 1:N_symbols
            start_idx   = (k-1)*(N_fft+N_cp)+1;
            sym_cp      = rx_rayleigh(start_idx:start_idx+N_fft+N_cp-1);
            sym         = sym_cp(N_cp+1:end);
            freq_sym    = fft(sym, N_fft) / sqrt(N_fft);
            % Perfect ZF equalization using known channel
            freq_sym_eq = freq_sym ./ H_freq_true;
            data_eq     = freq_sym_eq(data_idx);
            rx_bits_ray = [rx_bits_ray; qpsk_demodulate(data_eq)];
        end

        errors_rayleigh = errors_rayleigh + sum(tx_bits ~= rx_bits_ray);

        % Save Rayleigh constellation
        if any(snr_idx == plot_snr_idx) && trial == 1
            plot_pos = find(plot_snr_idx == snr_idx);
            k_mid = round(N_symbols/2);
            start_idx = (k_mid-1)*(N_fft+N_cp)+1;
            sym_cp = rx_rayleigh(start_idx:start_idx+N_fft+N_cp-1);
            sym = sym_cp(N_cp+1:end);
            freq_sym = fft(sym,N_fft)/sqrt(N_fft);
            freq_sym_eq = freq_sym ./ H_freq_true;
            const_symbols{plot_pos,2} = freq_sym_eq(data_idx);
        end

        % ─────────────────────────────────────────────
        %  SCENARIO C: Rayleigh fading + pilot-based EQ
        % ─────────────────────────────────────────────
        rx_bits_pilot = [];
        for k = 1:N_symbols
            start_idx   = (k-1)*(N_fft+N_cp)+1;
            sym_cp      = rx_rayleigh(start_idx:start_idx+N_fft+N_cp-1);
            sym         = sym_cp(N_cp+1:end);
            freq_sym    = fft(sym, N_fft) / sqrt(N_fft);
            % Pilot-based LS estimation
            H_pilots    = freq_sym(pilot_idx) ./ (1+0j); % known pilot=1
            H_interp    = interp1(pilot_idx(:), H_pilots(:), ...
                                  (1:N_fft)', 'linear', 'extrap');
            freq_sym_eq = freq_sym ./ H_interp;
            data_eq     = freq_sym_eq(data_idx);
            rx_bits_pilot = [rx_bits_pilot; qpsk_demodulate(data_eq)];
        end

        errors_pilot = errors_pilot + sum(tx_bits ~= rx_bits_pilot);

    end % trial loop

    total = total_bits * n_trials;
    BER_awgn(snr_idx)     = errors_awgn     / total;
    BER_rayleigh(snr_idx) = errors_rayleigh  / total;
    BER_pilot(snr_idx)    = errors_pilot     / total;

    fprintf('Eb/N0 = %2d dB | BER_AWGN = %.2e | BER_Rayleigh = %.2e | BER_Pilot = %.2e\n', ...
            EbN0_dB(snr_idx), BER_awgn(snr_idx), BER_rayleigh(snr_idx), BER_pilot(snr_idx));
end

% =========================================================================
%  SECTION 7 — THEORETICAL BER
% =========================================================================
% QPSK in AWGN: BER = Q(sqrt(2*Eb/N0)) = 0.5*erfc(sqrt(Eb/N0))
BER_theory_awgn = 0.5 * erfc(sqrt(EbN0_linear));

% QPSK over flat Rayleigh fading (theoretical):
% BER = 0.5 * (1 - sqrt(EbN0 / (1 + EbN0)))
BER_theory_rayleigh = 0.5 * (1 - sqrt(EbN0_linear ./ (1 + EbN0_linear)));

% =========================================================================
%  SECTION 8 — PLOTS
% =========================================================================

%% --- PLOT 1: BER vs Eb/N0 ---
figure('Name', 'BER vs Eb/N0', 'NumberTitle', 'off', ...
       'Position', [100, 100, 800, 500]);

semilogy(EbN0_dB, BER_theory_awgn,    'k--',  'LineWidth', 2); hold on;
semilogy(EbN0_dB, BER_theory_rayleigh,'r--',  'LineWidth', 2);
semilogy(EbN0_dB, BER_awgn,           'bo-',  'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB, BER_rayleigh,       'rs-',  'LineWidth', 1.5, 'MarkerSize', 7);
semilogy(EbN0_dB, BER_pilot,          'g^-',  'LineWidth', 1.5, 'MarkerSize', 7);

% Replace zero BERs (can't plot on log scale) with NaN
BER_awgn(BER_awgn == 0)         = NaN;
BER_rayleigh(BER_rayleigh == 0) = NaN;
BER_pilot(BER_pilot == 0)       = NaN;

grid on;
xlabel('E_b/N_0 (dB)', 'FontSize', 13);
ylabel('Bit Error Rate (BER)', 'FontSize', 13);
title('OFDM BER Performance — QPSK', 'FontSize', 14);
legend({'Theory: AWGN', ...
        'Theory: Rayleigh', ...
        'Sim: AWGN + perfect EQ', ...
        'Sim: Rayleigh + perfect EQ', ...
        'Sim: Rayleigh + pilot EQ'}, ...
       'Location', 'southwest', 'FontSize', 11);
xlim([EbN0_dB(1), EbN0_dB(end)]);
ylim([1e-5, 1]);
set(gca, 'FontSize', 11);

% Annotation box
annotation('textbox', [0.55 0.7 0.35 0.18], ...
    'String', {sprintf('N_{fft} = %d', N_fft), ...
               sprintf('N_{cp}  = %d', N_cp), ...
               sprintf('N_{pilots} = %d', N_pilots), ...
               sprintf('Modulation: QPSK')}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', ...
    'EdgeColor', [0.7 0.7 0.7], 'FontSize', 10);

%% --- PLOT 2: Constellation diagrams (AWGN at 3 SNR levels) ---
figure('Name', 'Constellation Diagrams', 'NumberTitle', 'off', ...
       'Position', [100, 650, 900, 320]);

snr_labels = {sprintf('Eb/N0 = %d dB (Low)', EbN0_dB(plot_snr_idx(1))), ...
              sprintf('Eb/N0 = %d dB (Mid)', EbN0_dB(plot_snr_idx(2))), ...
              sprintf('Eb/N0 = %d dB (High)', EbN0_dB(plot_snr_idx(end)))};

for p = 1:3
    subplot(1, 3, p);
    syms = const_symbols{p, 1};
    if ~isempty(syms)
        scatter(real(syms), imag(syms), 15, 'b.', 'MarkerEdgeAlpha', 0.5);
    end
    hold on;
    % Ideal QPSK points
    ideal = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);
    scatter(real(ideal), imag(ideal), 80, 'r', 'filled', 'Marker', 'x', ...
            'LineWidth', 2);
    grid on;
    axis equal;
    xlim([-2 2]); ylim([-2 2]);
    xlabel('In-phase (I)', 'FontSize', 10);
    ylabel('Quadrature (Q)', 'FontSize', 10);
    title(snr_labels{p}, 'FontSize', 10);
    xline(0, 'k--', 'LineWidth', 0.5, 'Alpha', 0.4);
    yline(0, 'k--', 'LineWidth', 0.5, 'Alpha', 0.4);
end
sgtitle('QPSK Constellation — AWGN Channel (red × = ideal points)', 'FontSize', 13);

%% --- PLOT 3: Channel frequency response ---
figure('Name', 'Channel Response', 'NumberTitle', 'off', ...
       'Position', [950, 100, 700, 400]);

subplot(2,1,1);
H_full = fft(h_channel, N_fft);
plot(0:N_fft-1, 20*log10(abs(H_full)), 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Subcarrier index', 'FontSize', 11);
ylabel('|H(k)| (dB)', 'FontSize', 11);
title('Rayleigh Channel — Frequency Response (Magnitude)', 'FontSize', 12);
xlim([0 N_fft-1]);

subplot(2,1,2);
stem(0:channel_len-1, abs(h_channel), 'b', 'filled', 'LineWidth', 1.5);
grid on;
xlabel('Tap index', 'FontSize', 11);
ylabel('|h(n)|', 'FontSize', 11);
title('Channel Impulse Response (3 taps)', 'FontSize', 12);
xlim([-0.5 channel_len-0.5]);

%% --- PLOT 4: OFDM Time-domain signal (one symbol) ---
figure('Name', 'Time Domain Signal', 'NumberTitle', 'off', ...
       'Position', [950, 550, 700, 350]);

% Generate one OFDM symbol for visualization
bits_demo = randi([0 1], bits_per_ofdm_sym, 1);
tx_demo   = ofdm_transmit(bits_demo, pilot_idx, data_idx, N_fft, N_cp, 1);
sym_len   = N_fft + N_cp;

t = 0:sym_len-1;
subplot(2,1,1);
plot(t, real(tx_demo(1:sym_len)), 'b-', 'LineWidth', 1);
hold on;
xline(N_cp, 'r--', 'LineWidth', 1.5);
text(N_cp/2, max(real(tx_demo))*0.8, 'CP', 'Color', 'r', ...
     'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
grid on;
xlabel('Sample index', 'FontSize', 11);
ylabel('Amplitude', 'FontSize', 11);
title('OFDM Symbol Time Domain — Real Part (with Cyclic Prefix)', 'FontSize', 12);

subplot(2,1,2);
[pxx, f] = pwelch(tx_demo, [], [], N_fft*4, 1);
plot(f, 10*log10(pxx), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Normalized frequency', 'FontSize', 11);
ylabel('PSD (dB/Hz)', 'FontSize', 11);
title('Power Spectral Density of OFDM Signal', 'FontSize', 12);

% =========================================================================
%  SECTION 9 — SUMMARY STATISTICS
% =========================================================================

fprintf('\n=== Simulation Summary ===\n');
fprintf('Theoretical QPSK AWGN BER at 10 dB Eb/N0 : %.4e\n', ...
        0.5*erfc(sqrt(10^(10/10))));
[~, idx10] = min(abs(EbN0_dB - 10));
fprintf('Simulated   QPSK AWGN BER at 10 dB Eb/N0 : %.4e\n', BER_awgn(idx10));
fprintf('\nGap between pilot EQ and perfect EQ at 10 dB: %.2f dB\n', ...
        10*log10(BER_pilot(idx10) / max(BER_rayleigh(idx10), 1e-10)));
fprintf('\nAll plots generated. Simulation complete.\n');
