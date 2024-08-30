# Digital-Comm---Matched-Filters
A digital communication project that aims to explore matched filters.

# Matched Filters, Correlators, ISI, and Raised Cosine Filters in Digital Communications
## Overview
This project focuses on the simulation and analysis of key concepts in digital communications, including matched filters, correlators, intersymbol interference (ISI), and raised cosine filters. The project involves designing and analyzing systems that utilize these components to optimize signal transmission and reception in a noise environment.

## Project Objectives
### 1. Matched Filters and Correlators:

- Implement and compare the performance of matched filters and correlators in a noise-free environment.
- Analyze the outputs of these filters at the sampling instants to evaluate their effectiveness.
### 2. Noise Analysis:

- Study the impact of Additive White Gaussian Noise (AWGN) on the system's performance.
- Calculate the Bit Error Rate (BER) as a function of the energy-per-bit to noise power spectral density ratio (Eb/N0) for different filters.
- Compare the simulation results with theoretical BER values.
### 3. Intersymbol Interference (ISI) and Raised Cosine Filters:

- Investigate the effect of ISI by generating and analyzing eye patterns.
- Utilize square root raised cosine filters at both the transmitter and receiver to observe the impact on the overall system performance.
## Implementation Details
### 1. Matched Filters and Correlators in a Noise-Free Environment
- Generate a random binary signal using Pulse Amplitude Modulation (PAM) with binary polar signaling (+1, -1).
- Simulate the transmission system using a given pulse shaping function, sampled at discrete intervals.
- Compare the outputs of a matched filter and a correlator using MATLAB.
### 2. Noise Analysis
- Simulate the system with a binary polar signaling PAM signal in the presence of AWGN.
- Calculate and plot the BER against different values of Eb/N0 using MATLAB.
- Compare the BER results with theoretical values.
### 3. ISI and Raised Cosine Filters
- Simulate the transmission and reception of signals using square root raised cosine filters with different roll-off factors and delays.
- Generate and analyze eye patterns for various scenarios to study the effect of ISI.
