# OFDM Simulator from Scratch (MATLAB)

A complete OFDM transceiver chain implemented from first principles in MATLAB — no Communications Toolbox shortcuts. Built to understand and validate the physical-layer building blocks behind 4G/5G wireless systems.

## Pipeline

```
Bits → QPSK Modulator → OFDM Symbol Mapping → IFFT → Add Cyclic Prefix
     → Channel (AWGN / Rayleigh Multipath) → Remove CP → FFT
     → Equalization (Perfect CSI / Pilot-based LS) → QPSK Demodulator
     → Recovered Bits → BER Calculation
```

## Features

- QPSK modulation and demodulation with Gray coding, implemented manually
- Configurable OFDM parameters: FFT size (64), cyclic prefix length (16), pilot subcarriers (8)
- AWGN channel model with exact SNR-to-noise-variance conversion
- 3-tap Rayleigh multipath fading channel
- Two equalization strategies: zero-forcing with perfect channel knowledge, and pilot-based least-squares channel estimation with interpolation
- BER vs Eb/N0 curves benchmarked against theoretical QPSK formulas (AWGN and Rayleigh)
- Constellation diagrams at low / mid / high SNR
- Channel frequency response and time-domain OFDM waveform visualization

## Results

**BER vs Eb/N0** — simulated AWGN curve closely tracks the theoretical Q-function curve, validating the modulation/IFFT/FFT chain. Rayleigh fading shows the expected slower BER decay, and pilot-based equalization shows a realistic gap versus perfect channel knowledge.

![BER plot](results/ber_vs_ebn0.png)

**Constellation diagrams** — QPSK symbols visibly tighten around ideal points as SNR increases from 0 dB to 24 dB.

![Constellation](results/constellation.png)

**Channel frequency response** — the 3-tap Rayleigh channel produces frequency-selective fading across subcarriers, illustrating why per-subcarrier equalization is necessary.

![Channel response](results/channel_response.png)

## How to Run

Requires MATLAB (tested on R2021a+) with no additional toolboxes.

```matlab
ofdm_simulator
```

Running the script generates all four figures and prints BER values per SNR point to the console.

## System Parameters

| Parameter | Value |
|---|---|
| FFT size (N_fft) | 64 |
| Cyclic prefix (N_cp) | 16 |
| Pilot subcarriers | 8 |
| Modulation | QPSK |
| Channel | AWGN, 3-tap Rayleigh |
| Eb/N0 range | 0–24 dB |

## Background

This project was built to deepen understanding of OFDM — the foundation of WiFi, LTE, and 5G NR — as part of ongoing research on Integrated Sensing and Communication (ISAC) systems.

## Possible Extensions

- 16-QAM / 64-QAM modulation
- PAPR analysis (CCDF)
- Carrier frequency offset and synchronization
- Convolutional coding with Viterbi decoding
- 2x2 MIMO with Alamouti space-time coding
