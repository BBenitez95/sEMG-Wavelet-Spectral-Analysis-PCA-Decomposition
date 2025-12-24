# sEMG Wavelet Spectral Analysis & PCA Decomposition

**Version 1.0**  
**Author:** Brian Benitez  
**GitHub:** https://github.com/BBenitez95

---

## Overview

This R Shiny application performs time-frequency decomposition of surface electromyography (sEMG) signals using the von Tscharner wavelet methodology (von Tscharner, 2000). The tool provides:

1. **Wavelet Transform**: Decomposes sEMG signals into spectral intensity distributions across 11 frequency bands (~7-430 Hz)
2. **PCA Decomposition**: Identifies population-level patterns of spectral variation and quantifies where individual recordings fall along the primary axis of spectral variation

---

## Interpretive Framework

### What This Tool Does (Reliably)

This application performs **signal characterization**; it describes the spectral content of sEMG recordings and identifies consistent patterns of spectral variation across a dataset. The mathematical operations are valid and reproducible:

- Time-frequency decomposition using validated wavelet methodology
- Principal component analysis to identify dominant modes of spectral variation  
- Quantification of spectral balance (high vs. low frequency content)
- Tracking of spectral shifts within subjects across conditions or time

### What This Tool Does NOT Do

This tool **cannot** reliably infer:

- Motor unit recruitment strategies
- Muscle fiber type composition
- The relative activation of "fast-twitch" vs. "slow-twitch" motor units

### The Scientific Context

The relationship between sEMG spectral properties and underlying motor unit characteristics is **not established** and remains scientifically contested. Farina, Merletti, and Enoka (2014) provide a comprehensive analysis of why spectral analysis of interference EMG signals cannot reliably extract neural information.

The fundamental problem is that the surface EMG signal contains **inseparably mixed information** about:

1. **Neural drive** (motor neuron discharge times)
2. **Action potential shapes** (determined by fiber membrane properties, conduction velocity)
3. **Volume conductor effects** (tissue filtering, electrode position, subcutaneous fat thickness)

A spectral difference between two recordings could reflect *any* combination of these factors. For example, a motor unit with higher conduction velocity may produce a surface EMG signal spectrally indistinguishable from a lower-velocity motor unit that happens to be positioned closer to the recording electrodes.

---

## Appropriate Terminology

Throughout this application and in your own reporting, use **descriptive signal terminology** rather than physiological inference:

| Instead of... | Use... |
|---------------|--------|
| "Fast-twitch fiber contribution" | "High-frequency spectral content" |
| "Slow-twitch fiber activation" | "Low-frequency spectral content" |
| "Motor unit recruitment shift" | "Spectral centroid shift" |
| "Fiber type composition" | "Spectral distribution" |
| "The muscle recruited more fast-twitch fibers" | "The recording showed increased high-frequency content" |

The NNLS coefficients labeled "High" and "Low" in this application refer to **spectral characteristics of the signal**, not properties of the underlying motor units or muscle fibers.

---

## Installation

### Requirements
- R ≥ 4.0
- RStudio (recommended)

### Dependencies
The application uses `pacman` for automatic package management. On first run, required packages will be installed automatically:

```r
# Core packages installed automatically:
shiny, shinyjs, plotly, readxl, openxlsx, data.table, signal, dplyr, nnls
```

### Running the Application

```r
# From RStudio: Open sEMG_Spectral_Analysis_App.R and click "Run App"

# Or from R console:
shiny::runApp("path/to/sEMG_Spectral_Analysis_App.R")
```

---

## Methodology

### The von Tscharner Wavelet

The wavelet is defined in the frequency domain as:

```
ψ(f) = (f/fc)^η × exp[(1 - f/fc) × η]
```

where `fc` is the center frequency and `η = s × fc` is the shape parameter (default s = 0.3).

**Key properties:**
- **Asymmetric frequency response**: Sharper low-frequency roll-off reduces motion artifact contamination
- **Constant relative bandwidth**: Uniform resolution on logarithmic frequency scale
- **Near-optimal uncertainty**: BW × Δt products of 0.75-0.92

### Filter Bank Configuration

Center frequencies follow a power-law distribution:

```
fc(j) = (q + j)^r / scale,  j = 0, 1, ..., J-1
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| J | 11 | Number of wavelets |
| q | 1.45 | Frequency offset |
| r | 1.959 | Frequency exponent |
| scale | 0.3 | Global scale factor |

With defaults, the filter bank spans approximately 7-400 Hz.

### Processing Pipeline

1. **Bandpass Filtering** (optional): 4th-order Butterworth, default 10-500 Hz
2. **Wavelet Transform**: Convolution with 11-band filter bank in frequency domain
3. **Center-Band Renormalization**: Prevents double-counting in overlapping frequency regions
4. **Gaussian Temporal Smoothing**: Band-specific smoothing (σ = 3/8 × Δt)
5. **Intensity Extraction**: Mean power within analysis window

### PCA Decomposition

- **Non-centered PCA**: Preserves non-negative spectral structure
- **Eigenvector orientation**: Bro method for consistent sign assignment
- **Boundary spectra**: Extreme spectral shapes constrained to remain non-negative
- **NNLS decomposition**: Non-negative coefficients summing to 1

---

## Output Metrics

### Wavelet Transform Tab

| Metric | Description | Notes |
|--------|-------------|-------|
| Total Intensity | Sum of spectral intensities | Proportional to signal power; affected by electrode placement, gain |
| Mean Frequency | Intensity-weighted average frequency | Σ(fc × I) / Σ(I) |

### PCA & Decomposition Tab

| Metric | Description | Notes |
|--------|-------------|-------|
| Theta | Angle in PC1-PC2 plane | Describes position along primary axis of spectral variation |
| High Coefficient | NNLS weight for high-frequency component | Describes spectral balance, NOT fiber type |
| Low Coefficient | NNLS weight for low-frequency component | Describes spectral balance, NOT fiber type |
| S_rec variance | Variance explained by 2 PCs | Should be >0.95 |
| S_gen variance | Variance explained by wavelet fit | Should be close to S_rec |

---

## Quality Control

### Indicators of Valid Analysis

- S_rec variance > 0.95 (two-component model is appropriate)
- |S_rec - S_gen| < 0.05 (wavelet fits capture boundary spectra)
- fc_low < fc_high (proper component orientation)
- Spectra appear smooth without obvious artifacts

---

## Recommended Reporting Practices

When publishing results from this tool:

1. **Describe spectral features, not physiological mechanisms**
   - ✓ "High-frequency spectral content increased with contraction intensity"
   - ✗ "Fast-twitch motor units were increasingly recruited"

2. **Acknowledge limitations explicitly**
   - Note that spectral characteristics are influenced by multiple factors beyond motor unit properties
   - Cite the Farina et al. (2014) review regarding interpretive constraints

3. **Report as one component of a multi-method approach**
   - Combine with other measures (force, kinematics, perceived effort) for context
   - Consider within-subject designs to control anatomical variability

4. **Provide methodological details**
   - Sample rate, filter settings, analysis window parameters
   - Wavelet configuration (or note if defaults used)

---

## Limitations

### Methodological Constraints

- **Cross-subject variability**: Electrode placement, subcutaneous tissue thickness, and muscle geometry differ between individuals, confounding between-subject comparisons
- **Volume conductor effects**: The same motor unit activity can produce different spectral signatures depending on electrode position relative to active fibers
- **Two-component assumption**: Spectral variation may be more complex than a single high-low axis

### Interpretive Constraints

- **Fiber type inference is not established**: The relationship between sEMG spectral content and muscle fiber composition lacks convincing empirical support (Farina et al., 2014)
- **Causal attribution requires caution**: Spectral shifts could reflect changes in motor unit recruitment, discharge rate modulation, fatigue-related conduction velocity changes, or technical factors
- **Non-unique inverse problem**: Many different combinations of motor unit activity can produce similar sEMG signals

---

## Alternative Approaches

For researchers seeking to extract neural information from surface EMG, Farina et al. (2014) recommend **high-density surface EMG decomposition**:

- Uses arrays of tens to hundreds of electrodes
- Decomposes interference signal into individual motor unit discharge trains
- Separates neural code (discharge times) from action potential shapes
- Requires validated decomposition algorithms
- Enables direct measurement of motor unit behavior rather than inference from global signal features

---

## References

### Primary Methodological Source

1. von Tscharner, V. (2000). Intensity analysis in time-frequency space of surface myoelectric signals by wavelets of specified resolution. *Journal of Electromyography and Kinesiology*, 10(6), 433-445.

### Critical Perspectives on EMG Spectral Interpretation

2. **Farina, D., Merletti, R., & Enoka, R. M. (2014). The extraction of neural strategies from the surface EMG: an update. *Journal of Applied Physiology*, 117(11), 1215-1230.**
   - Comprehensive review of limitations in extracting neural information from surface EMG
   - Concludes that spectral analysis for motor unit characterization is "meaningless for the vast majority of situations"

3. Farina, D., Merletti, R., & Enoka, R. M. (2004). The extraction of neural strategies from the surface EMG. *Journal of Applied Physiology*, 96(4), 1486-1495.

4. Enoka, R. M., & Duchateau, J. (2015). Inappropriate interpretation of surface EMG signals and muscle fiber characteristics impedes understanding of the control of neuromuscular function. *Journal of Applied Physiology*, 119(12), 1516-1518.

### The Point-Counterpoint Debate

5. von Tscharner, V., & Nigg, B. M. (2008). Point: Spectral properties of the surface EMG can characterize motor unit recruitment strategies and muscle fiber type. *Journal of Applied Physiology*, 105(5), 1671-1673.

6. Farina, D. (2008). Counterpoint: Spectral properties of the surface EMG do not provide information about motor unit recruitment and muscle fiber type. *Journal of Applied Physiology*, 105(5), 1673-1674.

### Additional Methodological References

7. von Tscharner, V., & Goepfert, B. (2006). Estimation of the interplay between groups of fast and slow muscle fibers of the tibialis anterior and gastrocnemius muscle while running. *Journal of Electromyography and Kinesiology*, 16(2), 188-197.

8. Bro, R., Acar, E., & Kolda, T. G. (2008). Resolving the sign ambiguity in the singular value decomposition. *Journal of Chemometrics*, 22(2), 135-140.

9. Wakeling, J. M., Pascual, S. A., Nigg, B. M., & von Tscharner, V. (2001). Surface EMG shows distinct populations of muscle activity when measured during sustained sub-maximal exercise. *European Journal of Applied Physiology*, 86(1), 40-47.

---

## License

MIT License

---

## Citation

If using this tool in published research, please cite:

1. The original methodology: von Tscharner (2000)
2. This implementation: Benitez, B. (2024). sEMG Wavelet Spectral Analysis & PCA Decomposition. GitHub repository.
3. Consider citing Farina et al. (2014) when discussing interpretive limitations

---

*This tool is provided for research purposes. Users are responsible for appropriate interpretation of results within the constraints described above.*
