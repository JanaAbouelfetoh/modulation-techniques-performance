# Digital Communications Final Project
**Performance of Different Modulation Techniques**

> Begad Mohamed (8584) · Hadeer Hany (8657) · Jana Waleed (8881) — May 2026

---

## Project Structure

```
├── part2_MASK.m                    # MATLAB simulation — Part II (M-ASK)
├── Part2_MASK_BER.png              # BER figure for Part II
├── Digital_Comms_Report_Complete.tex  # Full LaTeX report (all 3 parts)
└── README.md
```

> **Note:** Part I and Part III figures (`Fig1_Manual_Performance.png`, `Fig2_Built-in.png`, `Screenshot 2026-05-15 023329.png`) are generated separately by teammates and must be placed in the same folder as the `.tex` file before compiling.

---

## Parts Overview

### Part I — OOK, PRK, and BFSK
Compares BER performance of On-Off Keying, Phase-Reversal Keying, and Binary FSK using both a manual signal-space implementation and MATLAB built-in verification functions.

### Part II — M-ASK (M = 2, 4, 8)
Simulates and compares BER for 2-ASK, 4-ASK, and 8-ASK with a fixed minimum distance constraint ($d_{\min} = 2$). Includes both simulated and theoretical curves.

### Part III — BPSK, QPSK, 4QAM, and 4ASK
Compares rectangular BPSK, QPSK, and 4QAM under the same $d_{\min} = 2$ constraint, with the 4-ASK curve from Part II included for reference.

---

## Running Part II

### Requirements
- MATLAB (any recent version)
- **No toolboxes required** — base MATLAB only

### How to Run
1. Open `part2_MASK.m` in MATLAB.
2. Run the script (`F5` or click **Run**).
3. The simulation will print progress to the console and display the BER figure when complete.

### Simulation Parameters
| Parameter | Value |
|-----------|-------|
| Bits per SNR point | 10⁶ |
| SNR range | 0 – 60 dB (step 3 dB) |
| M values | 2, 4, 8 |
| Minimum distance | $d_{\min} = 2$ (all M) |
| Channel | AWGN |
| Detection | Nearest-neighbour (min. Euclidean distance) |

### Constellation Levels & Symbol Energies
| Modulation | Levels | Avg. Symbol Energy $E_s$ |
|------------|--------|--------------------------|
| 2-ASK | {−1, +1} | 1 |
| 4-ASK | {−3, −1, +1, +3} | 5 |
| 8-ASK | {−7, −5, −3, −1, +1, +3, +5, +7} | 21 |

---

## Compiling the Report

The report is written in LaTeX and intended for **Overleaf**.

1. Upload all files to an Overleaf project.
2. Make sure all figure images are in the **same folder** as the `.tex` file.
3. Compile with **pdfLaTeX**.

---

## Key Results Summary

| Modulation | SNR for BER < 10⁻⁵ |
|------------|----------------------|
| 2-ASK | ~12 dB |
| 4-ASK | ~18 dB |
| 8-ASK | ~24 dB |

Each step up in M costs ~6 dB more SNR — a direct consequence of the fixed $d_{\min}$ constraint increasing average symbol energy.