#!/usr/bin/env Rscript
################################################################################
# sEMG Wavelet Spectral Analysis & PCA Decomposition
# Version 1.0
################################################################################
# Author: Brian Benitez
# GitHub: https://github.com/Brianben95
################################################################################

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  shiny, shinyjs,
  plotly, readxl, openxlsx, 
  data.table, signal, dplyr, nnls
)

# =============================================================================
# CONSTANTS & CONFIGURATION
# =============================================================================

# Default wavelet parameters (von Tscharner 2000)
DEFAULTS <- list(
  J = 11L,        # Number of wavelets
  q = 1.45,       # Frequency offset parameter
  r = 1.959,      # Frequency scaling exponent  
  scale = 0.3     # Global scale factor
)

# Gaussian equalizer uncertainty constant from von Tscharner (2000)
# Section 3.5.4: "a Gauss filter with a width of 3/8 the time-resolution"
UNCERTAINTY_CONST <- 0.375

# Professional monochromatic color palette
COLORS <- list(
  # Primary grays
  text_primary   = "#212529",
  text_secondary = "#495057",
  text_muted     = "#6c757d",
  
  # Backgrounds
  bg_page       = "#f8f9fa",
  bg_card       = "#ffffff",
  bg_input      = "#ffffff",
  
  # Borders
  border_light  = "#dee2e6",
  border_medium = "#adb5bd",
  
  # Accents (monochromatic)
  accent_dark   = "#343a40",
  accent_medium = "#495057",
  accent_light  = "#6c757d",
  
  # Plot colors (black/gray scale)
  plot_primary   = "#212529",
  plot_secondary = "#6c757d",
  plot_tertiary  = "#adb5bd",
  
  # High/Low frequency indicators
  high_freq = "#212529",
  low_freq  = "#868e96",
  
  # PC colors
  pc1 = "#212529",
  pc2 = "#adb5bd"
)

# =============================================================================
# MATHEMATICAL FUNCTIONS - VON TSCHARNER WAVELET TRANSFORM
# =============================================================================

#' Next power of 2
next_pow2 <- function(n) 2^ceiling(log2(n))

#' Time resolution (delta-t) for a Cauchy wavelet
#' @param fc Center frequency (Hz)
#' @param s Scale parameter
#' @return Time resolution in seconds
compute_dt <- function(fc, s) {
  eta <- s * fc
  if (eta <= 0) return(Inf)
  sqrt(eta * (exp(1/eta) - 1)) / (2 * pi * fc)
}

#' Cauchy wavelet kernel (frequency domain)
#' The von Tscharner wavelet is defined as:
#'   ψ(f) = (f/fc)^η × exp[(1 - f/fc) × η]
#' where η = s × fc
#' @param f Frequency vector (Hz)
#' @param fc Center frequency (Hz)
#' @param s Scale parameter (typically 0.3 for sEMG)
#' @return Wavelet amplitude at each frequency
cauchy_kernel <- function(f, fc, s) {
  if (fc <= 0) return(rep(0, length(f)))
  eta <- s * fc
  x <- f / fc
  result <- rep(0, length(f))
  valid <- (f >= 0) & (x > 0)
  result[valid] <- (x[valid]^eta) * exp((1 - x[valid]) * eta)
  result
}

#' Cauchy wavelet amplitude at single frequency
#' @param f Single frequency value (Hz)
#' @param fc Center frequency (Hz)
#' @param s Scale parameter
#' @return Wavelet amplitude
cauchy_amplitude <- function(f, fc, s) {
  if (f < 0 || fc <= 0) return(0)
  x <- f / fc
  eta <- s * fc
  (x^eta) * exp((1 - x) * eta)
}

#' Unit-normalized Cauchy wavelet
#' @param freq Frequency vector (Hz)
#' @param fc Center frequency (Hz)
#' @param s Scale parameter
#' @return Normalized wavelet (sums to 1)
cauchy_normalized <- function(freq, fc, s) {
  w <- sapply(freq, function(f) cauchy_amplitude(f, fc, s))
  w <- pmax(w, 0)
  total <- sum(w)
  if (total > 1e-14) w / total else rep(0, length(freq))
}

#' Gaussian equalizer kernel (frequency domain)
#' Used to smooth the intensity signal in time
#' @param f Frequency vector (Hz)
#' @param sigma_ms Sigma value in milliseconds (from empirical table)
#' @return Gaussian filter response
gaussian_kernel <- function(f, sigma_ms) {
  if (sigma_ms < 1e-3) return(rep(1, length(f)))
  sigma_s <- sigma_ms / 1000  # Convert to seconds
  exp(-(2 * pi * f * sigma_s)^2)
}

#' Band edges for a Cauchy wavelet
#' 
#' From von Tscharner (2000) Section 3.5.3:
#' "The bandwidth of the wavelet was defined as the width at 1/e of the 
#' maximum amplitude of the power spectrum."
#' 
#' Since power = amplitude², bandwidth at 1/e of power corresponds to
#' amplitude = sqrt(1/e) ≈ 0.6065
#'
#' @param fc Center frequency (Hz)
#' @param s Scale parameter
#' @return List with f_lo, f_hi, and bandwidth
band_edges <- function(fc, s) {
  eta <- fc * s
  if (eta <= 0) return(list(f_lo = 0, f_hi = 0, bw = 0))
  
  # Find where amplitude = sqrt(1/e) (power = 1/e)
  # Wavelet amplitude: W(x) = x^eta * exp((1-x)*eta), x = f/fc
  # Setting log(W) = log(sqrt(1/e)) = -0.5:
  # eta*log(x) + (1-x)*eta = -0.5
  g <- function(x) eta * log(x) + (1 - x) * eta + 0.5
  
  x1 <- tryCatch(
    uniroot(g, c(1e-6, 1 - 1e-9))$root,
    error = function(e) 0.5
  )
  x2 <- tryCatch(
    uniroot(g, c(1 + 1e-9, 10))$root,
    error = function(e) 1.5
  )
  
  list(f_lo = x1 * fc, f_hi = x2 * fc, bw = (x2 - x1) * fc)
}

#' Compute Gaussian equalizer sigma value
#' 
#' From von Tscharner (2000):
#' - Section 3.5.3: "The bandwidth of the wavelet was defined as the width at 
#'   1/e of the maximum amplitude of the power spectrum."
#' - The time-resolution dt is related to bandwidth through the uncertainty 
#'   principle: BW × dt ≈ constant (empirically ~0.75-0.92)
#' - Section 3.5.4: "a Gauss filter with a width of 3/8 the time-resolution"
#'
#' For default parameters (J=11, q=1.45, r=1.959, scale=0.3), we use the 
#' empirical uncertainty products from von Tscharner's Table 1.
#'
#' @param j Wavelet index (1-based, so j=1 corresponds to wavelet 0 in paper)
#' @param J Total number of wavelets
#' @param q Frequency offset parameter
#' @param r Frequency exponent
#' @param scale Scale parameter
#' @return List with sigma (ms), dt (ms), and bw (Hz)
compute_sigma_empirical <- function(j, J = 11, q = 1.45, r = 1.959, scale = 0.3) {
  # Empirical BW×dt products from von Tscharner (2000) Table 1
  # These are the uncertainty products for the default parameters
  UNCERTAINTY_PRODUCTS <- c(0.747, 0.922, 0.870, 0.861, 0.914, 
                            0.882, 0.914, 0.870, 0.879, 0.897, 0.867)
  
  # Handle NA/NULL inputs
  if (is.null(j) || is.na(j) || is.null(J) || is.na(J) ||
      is.null(q) || is.na(q) || is.null(r) || is.na(r) ||
      is.null(scale) || is.na(scale)) {
    return(list(sigma = NA, dt = NA, bw = NA))
  }
  
  # Center frequency for this wavelet
  fc <- (q + (j - 1))^r / scale
  
  # Compute bandwidth at 1/e of POWER (sqrt(1/e) of amplitude)
  # This matches the paper's definition in Section 3.5.3
  eta <- scale * fc
  
  # Find band edges where amplitude = sqrt(1/e)
  # Solving: x^eta * exp((1-x)*eta) = sqrt(1/e)
  # Taking log: eta*log(x) + (1-x)*eta = -0.5
  find_edge <- function(x) {
    if (x <= 0) return(-Inf)
    eta * log(x) + (1 - x) * eta + 0.5
  }
  
  x_lo <- tryCatch(uniroot(find_edge, c(1e-6, 1 - 1e-9))$root, error = function(e) 0.1)
  x_hi <- tryCatch(uniroot(find_edge, c(1 + 1e-9, 10))$root, error = function(e) 2)
  
  bw <- (x_hi - x_lo) * fc
  
  # Determine uncertainty product
  is_default <- isTRUE(abs(q - 1.45) < 0.001) && 
    isTRUE(abs(r - 1.959) < 0.001) && 
    isTRUE(abs(scale - 0.3) < 0.001) &&
    isTRUE(J == 11)
  
  if (is_default && j >= 1 && j <= 11) {
    # Use empirical uncertainty product from Table 1
    product <- UNCERTAINTY_PRODUCTS[j]
  } else {
    # For non-default parameters, use average uncertainty product
    product <- 0.875
  }
  
  # dt = product / BW (convert to ms)
  dt <- product / bw * 1000
  
  # sigma = 3/8 × dt (Section 3.5.4)
  sigma <- 0.375 * dt
  
  list(
    sigma = sigma,
    dt = dt,
    bw = bw
  )
}

#' Compute all sigma values for a wavelet filter bank
#' @param fs Sampling frequency (Hz) - not used but kept for API compatibility
#' @param J Number of wavelets
#' @param q Frequency offset parameter
#' @param r Frequency exponent
#' @param scale Scale parameter
#' @return Vector of sigma values in milliseconds
compute_all_sigmas <- function(fs = 2000, J = 11, q = 1.45, r = 1.959, scale = 0.3) {
  sapply(1:J, function(j) compute_sigma_empirical(j, J, q, r, scale)$sigma)
}

#' Legacy compute_sigma for backwards compatibility
#' Uses analytical approximation
compute_sigma <- function(fc, s) {
  UNCERTAINTY_CONST * compute_dt(fc, s) * 1000
}

#' von Tscharner Wavelet Transform
#' 
#' Implements the complete wavelet filter bank as described in:
#' von Tscharner V. (2000) J Electromyogr Kinesiol, 10(6):433-445
#'
#' @param sig Input signal vector
#' @param fs Sampling frequency (Hz)
#' @param J Number of wavelets (default 11)
#' @param q Frequency offset parameter (default 1.45)
#' @param r Frequency exponent (default 1.959)
#' @param scale Global scale factor (default 0.3)
#' @param renorm Apply center-band renormalization (default TRUE)
#' @param apply_gauss Apply Gaussian equalizer (default TRUE)
#' @return List containing:
#'   - powc: Power matrix (time x wavelets)
#'   - powspec: Total power spectrum (sum over time for each wavelet)
#'   - freqc: Center frequencies (Hz)
#'   - sigma_used: Sigma values used (ms)
#'   - dt_values: Time resolution for each wavelet (seconds)
von_tscharner_transform <- function(sig, fs, 
                                    J = DEFAULTS$J, 
                                    q = DEFAULTS$q, 
                                    r = DEFAULTS$r, 
                                    scale = DEFAULTS$scale,
                                    renorm = TRUE, 
                                    apply_gauss = TRUE) {
  n <- length(sig)
  n2 <- next_pow2(n)
  
  # Center frequencies: fc(j) = (q + j)^r / scale, j = 0, 1, ..., J-1
  fc <- (q + 0:(J-1))^r / scale
  
  # Zero-pad and remove DC
  sig_padded <- c(sig - mean(sig), rep(0, n2 - n))
  sig_fft <- fft(sig_padded)
  
  # Frequency axis
  freqs <- seq(0, fs - fs/n2, length.out = n2)
  freqs[freqs > fs/2] <- freqs[freqs > fs/2] - fs
  
  # Pre-compute wavelet kernels
  Q_mat <- matrix(0, n2, J)
  for (j in 1:J) {
    Q_mat[, j] <- cauchy_kernel(abs(freqs), fc[j], scale)
  }
  
  # Compute sigma values using iterative method (von Tscharner Section 3.5.3)
  sigma_vec <- compute_all_sigmas(fs, J, q, r, scale)
  
  # Pre-compute Gaussian equalizers
  G_mat <- matrix(1, n2, J)
  if (apply_gauss) {
    for (j in 1:J) {
      G_mat[, j] <- gaussian_kernel(abs(freqs), sigma_vec[j])
    }
  }
  
  # Time resolution derived from empirical sigma: dt = sigma / 0.375
  # (sigma is in ms, dt_values should be in seconds)
  dt_values <- sigma_vec / 0.375 / 1000
  
  # Initialize output
  powc <- matrix(0, n, J)
  
  # Process each wavelet
  for (j in 1:J) {
    Q <- Q_mat[, j]
    
    # Center-band renormalization
    # Divides each wavelet by sqrt of sum of squared amplitudes
    # from adjacent wavelets in the overlapping frequency region
    if (renorm && j > 1 && j < J) {
      mid_lo <- (fc[j-1] + fc[j]) / 2
      mid_hi <- (fc[j] + fc[j+1]) / 2
      mask <- (abs(freqs) > mid_lo) & (abs(freqs) <= mid_hi)
      
      if (any(mask)) {
        denom <- sqrt(
          Q_mat[mask, j-1]^2 + 
            Q_mat[mask, j]^2 + 
            Q_mat[mask, j+1]^2
        )
        denom[denom == 0] <- 1
        Q[mask] <- Q[mask] / denom
      }
    }
    
    # Band-limited signal
    band_fft <- sig_fft * Q
    band_time <- abs(fft(band_fft, inverse = TRUE)) / n2 * n2
    
    # Power (intensity) computation
    # The factor of 4 accounts for the one-sided spectrum
    band_power <- 4 * band_time^2 / n2^2
    
    # Apply Gaussian equalizer (temporal smoothing)
    if (apply_gauss) {
      power_fft <- fft(band_power)
      power_filtered <- Re(fft(power_fft * G_mat[, j], inverse = TRUE)) / n2
      powc[, j] <- power_filtered[1:n]
    } else {
      powc[, j] <- band_power[1:n]
    }
  }
  
  # Ensure non-negative power values
  powc <- pmax(powc, 0)
  
  list(
    powc = powc,
    powspec = colSums(powc),
    freqc = fc,
    sigma_used = sigma_vec,
    dt_values = dt_values
  )
}

# =============================================================================
# SIGNAL PROCESSING FUNCTIONS
# =============================================================================

#' Butterworth bandpass filter
#' @param x Input signal
#' @param fs Sampling frequency (Hz)
#' @param lo Low cutoff frequency (Hz)
#' @param hi High cutoff frequency (Hz)
#' @param order Filter order
#' @return Filtered signal
bandpass_filter <- function(x, fs, lo, hi, order = 4) {
  if (length(x) < 3 * (order + 1)) return(x)
  if (lo >= hi) return(x)
  if (hi >= fs/2) hi <- fs/2 * 0.99
  if (lo <= 0) lo <- 0.5
  
  tryCatch({
    bf <- signal::butter(order, c(lo/(fs/2), hi/(fs/2)), type = "pass")
    signal::filtfilt(bf, x)
  }, error = function(e) {
    warning("Filter failed: ", e$message)
    x
  })
}

#' Remove DC offset from signal
#' @param x Input signal
#' @param detrend_linear If TRUE, remove linear trend; otherwise just mean
#' @return Centered signal
remove_dc <- function(x, detrend_linear = FALSE) {
  if (length(x) < 2) return(x)
  if (detrend_linear) {
    # Linear detrend
    t <- seq_along(x)
    fit <- lm(x ~ t)
    x - fitted(fit)
  } else {
    x - mean(x, na.rm = TRUE)
  }
}

# =============================================================================
# PCA AND DECOMPOSITION FUNCTIONS
# =============================================================================

#' PCA without mean centering
#' Used for spectral data where the baseline (DC) is meaningful
#' @param A Data matrix (features x observations)
#' @return List with eigvec, eigval, and scores
pca_no_centering <- function(A) {
  N <- ncol(A)
  B <- (A %*% t(A)) / (N - 1)
  eig <- eigen(B, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  list(
    eigvec = eig$vectors[, ord],
    eigval = eig$values[ord],
    scores = t(eig$vectors[, ord]) %*% A
  )
}

#' Orient PCs using Bro's signed squared projections method
#' Applies Bro, Acar & Kolda (2008) to resolve eigenvector sign ambiguity
#' @param eigvec Eigenvector matrix
#' @param eigval Eigenvalue vector
#' @param scores Score matrix
#' @return Oriented PCA results
orient_pcs <- function(eigvec, eigval, scores) {
  # Bro et al. (2008): signed sum of squared projections
  # S = sum( sign(score_j) * score_j^2 )
  # Flip if S < 0
  
  # Orient PC1
  pc1_scores <- scores[1, ]
  S1 <- sum(sign(pc1_scores) * pc1_scores^2)
  if (S1 < 0) {
    eigvec[, 1] <- -eigvec[, 1]
    scores[1, ] <- -scores[1, ]
  }
  
  # Orient PC2
  pc2_scores <- scores[2, ]
  S2 <- sum(sign(pc2_scores) * pc2_scores^2)
  if (S2 < 0) {
    eigvec[, 2] <- -eigvec[, 2]
    scores[2, ] <- -scores[2, ]
  }
  
  list(eigvec = eigvec, eigval = eigval, scores = scores)
}

#' Fit Cauchy wavelets to PC-derived spectra
#' @param eigvec Eigenvector matrix (frequency x PCs)
#' @param freq Frequency vector
#' @return List of fitted wavelet parameters
fit_wavelets_to_pcs <- function(eigvec, freq) {
  PC1 <- eigvec[, 1]
  PC2 <- eigvec[, 2]
  
  pos_mask <- PC2 > 0
  neg_mask <- PC2 < 0
  
  a_lo <- if (any(pos_mask)) max(-PC1[pos_mask] / PC2[pos_mask]) else 0
  a_hi <- if (any(neg_mask)) min(-PC1[neg_mask] / PC2[neg_mask]) else 0
  
  # Create target spectra
  lo_raw <- PC1 + a_lo * PC2
  hi_raw <- PC1 + a_hi * PC2
  i_low <- lo_raw / sum(lo_raw)
  i_high <- hi_raw / sum(hi_raw)
  
  # Fit Cauchy wavelets using NLS
  fit_cauchy <- function(target) {
    tryCatch({
      fit <- nls(y ~ cauchy_normalized(freq, fc, s), 
                 data = data.frame(y = target, freq = freq),
                 start = list(fc = 100, s = 0.5), 
                 lower = c(1e-6, 0), 
                 upper = c(1e4, 10), 
                 algorithm = "port")
      coef(fit)
    }, error = function(e) {
      c(fc = sum(freq * abs(target)) / sum(abs(target)), s = 0.5)
    })
  }
  
  p_lo <- fit_cauchy(i_low)
  p_hi <- fit_cauchy(i_high)
  
  # Ensure low < high
  if (p_lo["fc"] > p_hi["fc"]) {
    tmp <- p_lo
    p_lo <- p_hi
    p_hi <- tmp
    tmp <- i_low
    i_low <- i_high
    i_high <- tmp
  }
  
  # Band edges
  e_lo <- band_edges(p_lo["fc"], p_lo["s"])
  e_hi <- band_edges(p_hi["fc"], p_hi["s"])
  
  list(
    fc_lo = unname(p_lo["fc"]),
    fc_hi = unname(p_hi["fc"]),
    s_lo = unname(p_lo["s"]),
    s_hi = unname(p_hi["s"]),
    wave_lo = cauchy_normalized(freq, p_lo["fc"], p_lo["s"]),
    wave_hi = cauchy_normalized(freq, p_hi["fc"], p_hi["s"]),
    i_low = i_low,
    i_high = i_high,
    a_lo = a_lo,
    a_hi = a_hi,
    flo_lo = e_lo$f_lo,
    fhi_lo = e_lo$f_hi,
    bw_lo = e_lo$bw,
    flo_hi = e_hi$f_lo,
    fhi_hi = e_hi$f_hi,
    bw_hi = e_hi$bw,
    dt_lo = compute_dt(p_lo["fc"], p_lo["s"]),
    dt_hi = compute_dt(p_hi["fc"], p_hi["s"])
  )
}

#' NNLS decomposition of spectra into wavelet basis
#' @param A Data matrix (frequencies x observations)
#' @param wave_hi High-frequency wavelet
#' @param wave_lo Low-frequency wavelet
#' @return List with coefficients, reconstructed spectra, and variance explained
nnls_decomposition <- function(A, wave_hi, wave_lo) {
  basis <- cbind(wave_hi, wave_lo)
  N <- ncol(A)
  coeffs <- matrix(0, 2, N)
  S_gen <- matrix(0, nrow(A), N)
  totals <- colSums(A)
  
  for (j in 1:N) {
    spectrum <- A[, j]
    if (totals[j] > 1e-12) {
      fit <- nnls::nnls(basis, spectrum / totals[j])
      S_gen[, j] <- totals[j] * (basis %*% fit$x)
    } else {
      fit <- nnls::nnls(basis, spectrum)
      S_gen[, j] <- basis %*% fit$x
    }
    c_sum <- sum(fit$x)
    coeffs[, j] <- if (c_sum > 1e-12) fit$x / c_sum else c(0.5, 0.5)
  }
  
  rownames(coeffs) <- c("High", "Low")
  var_gen <- 1 - sum((A - S_gen)^2) / sum(A^2)
  
  list(coeffs = coeffs, S_gen = S_gen, var_gen = var_gen)
}

# =============================================================================
# LABEL PARSING FUNCTIONS
# =============================================================================

#' Extract condition from observation name
#' @param name Observation name string
#' @return Condition string
extract_condition <- function(name) {
  parts <- strsplit(trimws(name), "_")[[1]]
  cond_parts <- parts[!grepl("\\d", parts)]
  if (length(cond_parts) == 0) return("Default")
  paste(cond_parts, collapse = "_")
}

#' Extract timepoint from observation name
#' @param name Observation name string
#' @return Numeric timepoint or NA
extract_timepoint <- function(name) {
  parts <- strsplit(trimws(name), "_")[[1]]
  for (p in parts) {
    if (grepl("^[Rr]ep\\d+$|^[Tt]ime\\d+$", p)) {
      return(as.numeric(gsub("\\D", "", p)))
    }
    if (grepl("^\\d+$", p) && as.numeric(p) < 100) {
      return(as.numeric(p))
    }
  }
  NA
}

#' Extract numeric value from string
#' @param val Input value (string or numeric)
#' @return Numeric value or NA
extract_numeric <- function(val) {
  if (is.numeric(val)) return(val)
  m <- regmatches(as.character(val), regexpr("[-+]?[0-9]*\\.?[0-9]+", as.character(val)))
  if (length(m) > 0) as.numeric(m) else NA
}

# =============================================================================
# PLOTTING HELPERS
# =============================================================================

#' Create empty placeholder plot
#' @param msg Message to display
#' @return Plotly object
empty_plot <- function(msg = "No data") {
  plot_ly(
    type = "scatter", 
    mode = "markers", 
    x = 0, y = 0, 
    marker = list(opacity = 0),
    hoverinfo = "none", 
    showlegend = FALSE
  ) %>%
    layout(
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE),
      annotations = list(
        list(
          text = msg, 
          x = 0.5, y = 0.5, 
          xref = "paper", yref = "paper",
          showarrow = FALSE, 
          font = list(size = 14, color = COLORS$text_muted)
        )
      ),
      plot_bgcolor = COLORS$bg_card,
      paper_bgcolor = COLORS$bg_card
    ) %>%
    config(displaylogo = FALSE)
}

#' Standard plot layout for consistency
plot_layout <- function(p, xlab = "", ylab = "", title = "") {
  p %>% layout(
    title = list(
      text = title,
      font = list(size = 12, color = COLORS$text_primary, family = "Inter, sans-serif"),
      x = 0.02,
      xanchor = "left"
    ),
    xaxis = list(
      title = list(text = xlab, font = list(size = 11, color = COLORS$text_secondary)),
      tickfont = list(size = 10, color = COLORS$text_secondary),
      gridcolor = COLORS$border_light,
      linecolor = COLORS$border_medium,
      zerolinecolor = COLORS$border_medium
    ),
    yaxis = list(
      title = list(text = ylab, font = list(size = 11, color = COLORS$text_secondary)),
      tickfont = list(size = 10, color = COLORS$text_secondary),
      gridcolor = COLORS$border_light,
      linecolor = COLORS$border_medium,
      zerolinecolor = COLORS$border_medium
    ),
    legend = list(
      font = list(size = 10, color = COLORS$text_secondary),
      bgcolor = "rgba(255,255,255,0.9)"
    ),
    plot_bgcolor = COLORS$bg_card,
    paper_bgcolor = COLORS$bg_card,
    margin = list(t = 40, b = 40, l = 60, r = 20)
  )
}

# =============================================================================
# APPLICATION STYLES
# =============================================================================

app_css <- '
@import url("https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap");

* { box-sizing: border-box; }

body {
  font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  background: #f8f9fa;
  color: #212529;
  font-size: 13px;
  line-height: 1.5;
}

.container-fluid { padding: 16px 24px; }

/* Title Panel */
.title-panel {
  background: #ffffff;
  padding: 16px 24px;
  margin: -16px -24px 20px -24px;
  border-bottom: 1px solid #dee2e6;
}

.title-panel h2 {
  font-size: 18px;
  font-weight: 600;
  color: #212529;
  margin: 0 0 4px 0;
}

.title-panel .subtitle {
  font-size: 12px;
  color: #6c757d;
}

/* Section Headers */
h4 {
  font-size: 11px;
  font-weight: 600;
  color: #212529;
  text-transform: uppercase;
  letter-spacing: 0.75px;
  margin: 0 0 12px 0;
  padding-bottom: 8px;
  border-bottom: 1px solid #dee2e6;
}

h5 {
  font-size: 14px;
  font-weight: 600;
  color: #212529;
  margin: 0 0 10px 0;
}

/* Cards/Wells */
.well {
  background: #ffffff;
  border: 1px solid #dee2e6;
  border-radius: 6px;
  box-shadow: 0 1px 3px rgba(0,0,0,0.04);
  padding: 16px;
  margin-bottom: 12px;
}

/* Form Controls */
.control-label {
  font-size: 11px;
  font-weight: 500;
  color: #495057;
  margin-bottom: 4px;
  display: block;
}

.form-control {
  font-size: 12px;
  padding: 8px 10px;
  border: 1px solid #ced4da;
  border-radius: 4px;
  background: #ffffff;
  transition: border-color 0.15s, box-shadow 0.15s;
}

.form-control:focus {
  border-color: #495057;
  box-shadow: 0 0 0 2px rgba(73,80,87,0.1);
  outline: none;
}

/* Buttons */
.btn {
  font-size: 11px;
  font-weight: 500;
  padding: 8px 14px;
  border-radius: 4px;
  border: none;
  transition: all 0.15s;
  text-transform: uppercase;
  letter-spacing: 0.3px;
}

.btn-primary {
  background: #212529;
  color: #ffffff;
}

.btn-primary:hover {
  background: #343a40;
}

.btn-success {
  background: #495057;
  color: #ffffff;
}

.btn-success:hover {
  background: #343a40;
}

.btn-info {
  background: #6c757d;
  color: #ffffff;
}

.btn-info:hover {
  background: #5a6268;
}

.btn-block {
  display: block;
  width: 100%;
  margin-bottom: 8px;
}

/* Plot Containers */
.plot-container {
  background: #ffffff;
  border: 1px solid #dee2e6;
  border-radius: 6px;
  padding: 16px;
  margin-bottom: 12px;
}

.plot-header {
  font-size: 12px;
  font-weight: 600;
  color: #212529;
  margin-bottom: 12px;
  padding-bottom: 8px;
  border-bottom: 1px solid #e9ecef;
}

/* Status Panel */
.status-panel {
  background: #f8f9fa;
  border: 1px solid #dee2e6;
  border-radius: 4px;
  padding: 12px;
  font-size: 11px;
  font-family: "SF Mono", "Consolas", monospace;
  color: #495057;
}

/* Info Panel */
.info-panel {
  background: #f8f9fa;
  border-left: 3px solid #6c757d;
  padding: 10px 12px;
  font-size: 11px;
  color: #495057;
  margin-bottom: 12px;
  border-radius: 0 4px 4px 0;
}

/* Observation Table */
.obs-table {
  max-height: 200px;
  overflow-y: auto;
  border: 1px solid #dee2e6;
  border-radius: 4px;
}

.obs-table table {
  width: 100%;
  font-size: 11px;
  border-collapse: collapse;
}

.obs-table th {
  background: #f8f9fa;
  padding: 8px 10px;
  text-align: left;
  font-weight: 600;
  border-bottom: 1px solid #dee2e6;
  position: sticky;
  top: 0;
}

.obs-table td {
  padding: 6px 10px;
  border-bottom: 1px solid #f1f3f5;
}

.obs-table tr:hover td {
  background: #f8f9fa;
}

/* Clickable observation list */
.obs-row {
  display: flex;
  padding: 8px 12px;
  border-bottom: 1px solid #f1f3f5;
  cursor: pointer;
  font-size: 11px;
  transition: background 0.1s;
}

.obs-row:hover {
  background: #f8f9fa;
}

.obs-row.selected {
  background: #495057;
  color: #ffffff;
}

.obs-row .obs-name {
  flex: 2;
  font-weight: 500;
}

.obs-row .obs-source {
  flex: 2;
  color: #6c757d;
}

.obs-row.selected .obs-source {
  color: #adb5bd;
}

.obs-row .obs-points {
  flex: 1;
  text-align: right;
  color: #6c757d;
}

.obs-row.selected .obs-points {
  color: #adb5bd;
}

.obs-header {
  display: flex;
  padding: 8px 12px;
  background: #f8f9fa;
  border-bottom: 1px solid #dee2e6;
  font-size: 10px;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.5px;
  color: #495057;
  position: sticky;
  top: 0;
}

.obs-header .obs-name { flex: 2; }
.obs-header .obs-source { flex: 2; }
.obs-header .obs-points { flex: 1; text-align: right; }

/* Tab Styling */
.nav-tabs {
  border-bottom: 1px solid #dee2e6;
}

.nav-tabs > li > a {
  font-size: 12px;
  font-weight: 500;
  color: #6c757d;
  padding: 10px 16px;
  border: none;
  border-bottom: 2px solid transparent;
  transition: color 0.15s, border-color 0.15s;
}

.nav-tabs > li > a:hover {
  color: #212529;
  background: transparent;
  border-color: transparent;
}

.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover,
.nav-tabs > li.active > a:focus {
  color: #212529;
  background: transparent;
  border: none;
  border-bottom: 2px solid #212529;
}

/* Documentation Styles */
.doc-section {
  margin-bottom: 24px;
}

.doc-section h5 {
  color: #212529;
  margin-bottom: 12px;
}

.doc-section p {
  font-size: 13px;
  line-height: 1.6;
  color: #495057;
  margin-bottom: 12px;
}

.doc-section ul {
  font-size: 13px;
  line-height: 1.6;
  color: #495057;
  padding-left: 24px;
}

.formula {
  background: #f8f9fa;
  padding: 12px 16px;
  border-left: 3px solid #495057;
  margin: 16px 0;
  font-family: "SF Mono", "Consolas", monospace;
  font-size: 13px;
}

.doc-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 12px;
  border: 1px solid #dee2e6;
  border-radius: 4px;
  overflow: hidden;
  margin: 12px 0;
}

.doc-table th,
.doc-table td {
  border: 1px solid #dee2e6;
  padding: 10px 12px;
  text-align: left;
}

.doc-table th {
  background: #495057;
  color: #ffffff;
  font-weight: 600;
}

.doc-table tr:nth-child(even) {
  background: #f8f9fa;
}

.doc-note {
  background: #f8f9fa;
  border-left: 3px solid #6c757d;
  padding: 12px 16px;
  margin: 16px 0;
  font-size: 12px;
  border-radius: 0 4px 4px 0;
}

.reference {
  font-size: 11px;
  color: #6c757d;
  padding-left: 24px;
  text-indent: -24px;
  margin-bottom: 8px;
  line-height: 1.5;
}

/* Checkbox Styling */
.checkbox label,
.shiny-input-checkboxgroup label {
  font-size: 12px;
  font-weight: 400;
}
'

# =============================================================================
# USER INTERFACE
# =============================================================================

ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(HTML(app_css))),
  
  div(class = "title-panel",
      h2("sEMG Wavelet Spectral Analysis"),
      div(class = "subtitle", "Signal Characterization Tool — von Tscharner Wavelet Transform & PCA Decomposition"),
      div(style = "font-size: 10px; color: #6c757d; margin-top: 4px;", 
          "Note: Spectral metrics describe signal characteristics, not motor unit or fiber type properties. See Documentation tab.")
  ),
  
  tabsetPanel(
    id = "tabs",
    
    # =========================================================================
    # TAB 1: WAVELET TRANSFORM
    # =========================================================================
    tabPanel("Wavelet Transform",
             sidebarLayout(
               sidebarPanel(width = 3,
                            
                            # Data Input
                            wellPanel(
                              h4("Data Input"),
                              fileInput("wt_file", "Excel File", accept = c(".xls", ".xlsx")),
                              div(class = "info-panel", 
                                  "Each sheet = one observation. Column 1 = time (s), Column 2 = signal.")
                            ),
                            
                            # Signal Processing
                            wellPanel(
                              h4("Signal Processing"),
                              numericInput("wt_fs", "Sample Rate (Hz)", 2000, 100, 50000),
                              numericInput("wt_offset", "Offset (samples)", 0, 0),
                              numericInput("wt_npts", "Window Length (samples)", 2048, 256, 32768),
                              checkboxInput("wt_filter", "Bandpass Filter", TRUE),
                              conditionalPanel(
                                "input.wt_filter",
                                fluidRow(
                                  column(6, numericInput("wt_flo", "Low (Hz)", 10, 1)),
                                  column(6, numericInput("wt_fhi", "High (Hz)", 500, 50))
                                ),
                                numericInput("wt_order", "Filter Order", 4, 2, 8)
                              ),
                              checkboxInput("wt_dc", "Remove DC Offset", TRUE)
                            ),
                            
                            # Wavelet Parameters
                            wellPanel(
                              h4("Wavelet Parameters"),
                              fluidRow(
                                column(6, numericInput("wt_J", "J (wavelets)", DEFAULTS$J, 2, 30)),
                                column(6, numericInput("wt_q", "q", DEFAULTS$q, 0.1, 10, 0.01))
                              ),
                              fluidRow(
                                column(6, numericInput("wt_r", "r", DEFAULTS$r, 0.1, 5, 0.001)),
                                column(6, numericInput("wt_scale", "scale", DEFAULTS$scale, 0.01, 2, 0.01))
                              ),
                              actionButton("wt_reset", "Reset to Defaults", class = "btn-info btn-block")
                            ),
                            
                            # Execute
                            wellPanel(
                              h4("Execute"),
                              actionButton("wt_run", "Run Transform", class = "btn-success btn-block")
                            ),
                            
                            # Status
                            wellPanel(
                              h4("Status"),
                              div(class = "status-panel", textOutput("wt_status"))
                            )
               ),
               
               mainPanel(width = 9,
                         # Wavelet Parameters Table
                         div(class = "plot-container",
                             div(class = "plot-header", "Wavelet Filter Bank Parameters"),
                             div(class = "obs-table", tableOutput("wt_params"))),
                         
                         # Observations Table (clickable)
                         div(class = "plot-container",
                             div(class = "plot-header", "Loaded Observations"),
                             div(class = "obs-table", uiOutput("wt_obs_list"))
                         ),
                         
                         # Signal and Spectrum
                         fluidRow(
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", textOutput("wt_sig_header", inline = TRUE)),
                                      plotlyOutput("wt_sig", height = "220px"))),
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", "Spectral Intensity"),
                                      plotlyOutput("wt_spec", height = "220px")))
                         ),
                         
                         # All Spectra
                         div(class = "plot-container",
                             div(class = "plot-header", "All Spectra Overlay"),
                             fluidRow(
                               column(12, 
                                      selectInput("wt_display", NULL, width = "200px",
                                                  c("Raw Intensities" = "raw", 
                                                    "Unit Area Normalized" = "unit", 
                                                    "Power Normalized" = "power"),
                                                  selected = "raw")
                               )
                             ),
                             plotlyOutput("wt_all", height = "260px")),
                         
                         # Export
                         wellPanel(
                           h4("Export Results"),
                           fluidRow(
                             column(6, textInput("wt_fname", "Output Filename", placeholder = "spectral_intensities")),
                             column(6, selectInput("wt_norm", "Data Format",
                                                   c("Raw Intensities" = "raw", 
                                                     "Unit Area Normalized" = "unit", 
                                                     "Total Intensity Normalized" = "subj_ti")))
                           ),
                           checkboxInput("wt_drop", "Exclude First Wavelet (low-frequency noise)", TRUE),
                           checkboxInput("wt_summary", "Include Summary Sheet (Total Intensity & Mean Frequency)", TRUE),
                           downloadButton("wt_download", "Download Excel", class = "btn-primary btn-block")
                         )
               )
             )
    ),
    
    # =========================================================================
    # TAB 2: PCA & DECOMPOSITION
    # =========================================================================
    tabPanel("PCA & Decomposition",
             sidebarLayout(
               sidebarPanel(width = 3,
                            
                            wellPanel(
                              h4("Load Spectral Data"),
                              fileInput("pca_file", "Spectral Intensities File", 
                                        accept = c(".csv", ".txt", ".xls", ".xlsx")),
                              uiOutput("pca_sheet_ui"),
                              div(class = "info-panel", 
                                  "Format: First column = frequencies (Hz), remaining columns = observations.")
                            ),
                            
                            wellPanel(
                              h4("Analysis Pipeline"),
                              actionButton("pca_run", "Run PCA", class = "btn-primary btn-block"),
                              actionButton("pca_decomp", "Fit Wavelets & Decompose", class = "btn-success btn-block")
                            ),
                            
                            wellPanel(
                              h4("Export Results"),
                              textInput("pca_fname", "Output Filename", placeholder = "pca_decomposition"),
                              downloadButton("pca_download", "Download Excel", class = "btn-primary btn-block")
                            ),
                            
                            wellPanel(
                              h4("Analysis Status"),
                              div(class = "status-panel", uiOutput("pca_status"))
                            )
               ),
               
               mainPanel(width = 9,
                         fluidRow(
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", "Principal Component Loadings"),
                                      plotlyOutput("pca_pc", height = "280px"))),
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", "Population Mean Spectrum"),
                                      plotlyOutput("pca_mean", height = "280px")))
                         ),
                         
                         fluidRow(
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", "Optimized Wavelet Fits"),
                                      plotlyOutput("pca_wave", height = "280px"))),
                           column(6, 
                                  div(class = "plot-container",
                                      div(class = "plot-header", "PCA Biplot by Condition"),
                                      plotlyOutput("pca_biplot", height = "280px")))
                         )
               )
             )
    ),
    
    # =========================================================================
    # TAB 3: DOCUMENTATION
    # =========================================================================
    tabPanel("Documentation",
             fluidRow(
               column(10, offset = 1,
                      style = "padding: 24px; background: #ffffff; border-radius: 6px; margin-top: 20px;
                   border: 1px solid #dee2e6;",
                      
                      # --- Section 1: Introduction ---
                      div(class = "doc-section",
                          h5("1. Introduction"),
                          p("This application implements the von Tscharner wavelet transform for surface 
               electromyography (sEMG) spectral analysis, following the methodology described 
               in von Tscharner (2000). The wavelet-based intensity analysis provides 
               time-frequency decomposition optimized for the non-stationary characteristics 
               of myoelectric signals during dynamic muscle contraction."),
                          p("The analysis pipeline consists of two stages: (1) wavelet transformation of 
               individual sEMG recordings to obtain spectral intensity distributions, and 
               (2) population-level Principal Component Analysis (PCA) to identify patterns of 
               spectral variation and quantify where each observation falls along the primary 
               axis of high-frequency versus low-frequency spectral content."),
                          div(class = "doc-note", style = "background: #fff3cd; border-left: 4px solid #856404; padding: 12px; margin: 12px 0;",
                              tags$strong("Important:"), " This tool performs ", tags$em("signal characterization:"), 
                              "It describes the spectral content of sEMG recordings. The 'High' and 'Low' 
                              frequency components refer to ", tags$strong("spectral characteristics of the signal"), 
                              ", not to muscle fiber types or motor unit properties. See Section 2 for 
                              critical interpretive guidance."
                          )
                      ),
                      
                      # --- Section 2: Critical Interpretive Framework ---
                      div(class = "doc-section",
                          h5("2. Interpretive Framework"),
                          div(class = "doc-note", style = "background: #fff3cd; border-left: 4px solid #856404; padding: 12px; margin: 12px 0;",
                              tags$strong("Important:"), 
                              " The relationship between sEMG spectral properties and underlying motor unit 
                              characteristics remains an area of active scientific debate. Spectral metrics 
                              should be interpreted as ", tags$em("signal characteristics"), " rather than 
                              direct measurements of motor unit recruitment or fiber type composition."
                          ),
                          
                          p(tags$strong("What this tool does (reliably):")),
                          tags$ul(
                            tags$li("Performs valid time-frequency decomposition of sEMG signals"),
                            tags$li("Identifies patterns of spectral variation across a dataset"),
                            tags$li("Quantifies spectral balance (high vs. low frequency content)"),
                            tags$li("Tracks spectral shifts within subjects across conditions")
                          ),
                          
                          p(tags$strong("Appropriate applications:")),
                          tags$ul(
                            tags$li(tags$strong("Fatigue monitoring:"), " Spectral compression during sustained 
                      contractions is an empirically robust phenomenon related to conduction velocity changes"),
                            tags$li(tags$strong("Within-subject comparisons:"), " When electrode placement, 
                      tissue properties, and geometry are controlled"),
                            tags$li(tags$strong("Exploratory analysis:"), " Identifying spectral patterns 
                      for further investigation"),
                            tags$li(tags$strong("Signal quality assessment:"), " Detecting artifacts or noise")
                          ),
                          
                          p(tags$strong("Recommended Terminology:")),
                          tags$table(class = "doc-table",
                                     tags$tr(tags$th("Instead of..."), tags$th("Use...")),
                                     tags$tr(tags$td("Fast-twitch fiber contribution"), 
                                             tags$td("High-frequency spectral content")),
                                     tags$tr(tags$td("Slow-twitch fiber activation"), 
                                             tags$td("Low-frequency spectral content")),
                                     tags$tr(tags$td("Motor unit recruitment shift"), 
                                             tags$td("Spectral centroid shift")),
                                     tags$tr(tags$td("Fiber type composition"), 
                                             tags$td("Spectral distribution"))
                          )
                      ),
                      
                      # --- Section 3: The von Tscharner Wavelet ---
                      div(class = "doc-section",
                          h5("3. The von Tscharner Wavelet"),
                          p("The wavelet is defined in the frequency domain as a modified Cauchy distribution:"),
                          div(class = "formula", HTML("&psi;(f) = (f / f<sub>c</sub>)<sup>&eta;</sup> &times; 
                                         exp[(1 - f / f<sub>c</sub>) &times; &eta;]")),
                          p(HTML("where f<sub>c</sub> is the center frequency and &eta; = s &times; f<sub>c</sub> 
                   is the shape parameter. The scale factor s (default 0.3) controls the 
                   time-frequency trade-off.")),
                          p("Key properties of this wavelet:"),
                          tags$ul(
                            tags$li(tags$strong("Asymmetric frequency response:"), " The sharper roll-off on 
                      the low-frequency side reduces contamination from motion artifacts and 
                      baseline drift, which predominantly occupy frequencies below 20 Hz."),
                            tags$li(tags$strong("Constant relative bandwidth:"), " Each wavelet spans 
                      approximately the same proportion of its center frequency, providing 
                      uniform resolution on a logarithmic frequency scale."),
                            tags$li(tags$strong("Near-optimal uncertainty:"), " The wavelet approaches the 
                      theoretical minimum for joint time-frequency uncertainty, with empirical 
                      BW x dt products ranging from 0.75 to 0.92.")
                          )
                      ),
                      
                      # --- Section 4: Filter Bank Configuration ---
                      div(class = "doc-section",
                          h5("4. Filter Bank Configuration"),
                          p("The J wavelets are distributed according to a power law that provides denser 
               sampling at lower frequencies:"),
                          div(class = "formula", HTML("f<sub>c</sub>(j) = (q + j)<sup>r</sup> / scale, 
                                         &nbsp;&nbsp;&nbsp; j = 0, 1, ..., J-1")),
                          tags$table(class = "doc-table",
                                     tags$tr(tags$th("Parameter"), tags$th("Default"), tags$th("Description")),
                                     tags$tr(tags$td("J"), tags$td("11"), 
                                             tags$td("Number of wavelets in the filter bank")),
                                     tags$tr(tags$td("q"), tags$td("1.45"), 
                                             tags$td("Offset parameter; controls the lowest center frequency")),
                                     tags$tr(tags$td("r"), tags$td("1.959"), 
                                             tags$td("Exponent; controls the spacing between center frequencies")),
                                     tags$tr(tags$td("scale"), tags$td("0.3"), 
                                             tags$td("Global scale factor; affects both frequency placement and bandwidth"))
                          ),
                          p("With default parameters, the filter bank spans approximately 7-400 Hz, covering 
               the primary frequency content of sEMG signals.")
                      ),
                      
                      # --- Section 5: Gaussian Equalizer ---
                      div(class = "doc-section",
                          h5("5. Gaussian Temporal Smoothing"),
                          p("Following wavelet convolution, a Gaussian equalizer is applied in the frequency 
               domain to smooth the instantaneous intensity signal in time. This reduces 
               high-frequency fluctuations in the power estimate while preserving the underlying 
               temporal structure of muscle excitation."),
                          p("The Gaussian width (sigma) is proportional to the time-resolution of each wavelet:"),
                          div(class = "formula", "sigma = (3/8) x time_resolution"),
                          p("For the default parameters (scale = 0.3), the empirical sigma values from 
               von Tscharner (2000) Table 1 are used directly. When parameters are modified, 
               sigma values are computed from the bandwidth-time resolution relationship.")
                      ),
                      
                      # --- Section 6: Center-Band Renormalization ---
                      div(class = "doc-section",
                          h5("6. Center-Band Renormalization"),
                          p("Adjacent wavelets in the filter bank overlap in frequency. To prevent 
               double-counting of spectral energy in overlapping regions, center-band 
               renormalization divides each wavelet amplitude by the root-sum-square of 
               itself and its two neighbors:"),
                          div(class = "formula", 
                              "Q'(f) = Q(f) / sqrt[Q(j-1)^2 + Q(j)^2 + Q(j+1)^2]"),
                          p("This procedure is applied only within the frequency region bounded by the 
               midpoints between adjacent center frequencies. Edge wavelets (j = 0 and j = J-1) 
               are not renormalized.")
                      ),
                      
                      # --- Section 7: PCA Decomposition ---
                      div(class = "doc-section",
                          h5("7. Principal Component Analysis"),
                          p("The PCA stage operates on the population of spectral intensity vectors. Unlike 
               standard PCA, the data are ", tags$em("not"), " mean-centered prior to 
               eigendecomposition. This preserves the non-negative nature of spectral intensities 
               and ensures that the principal components represent actual spectral shapes rather 
               than deviations from a mean spectrum."),
                          p("The covariance matrix is computed as:"),
                          div(class = "formula", "B = A * A' / (N - 1)"),
                          p("where A is the (frequencies x observations) data matrix. Eigenvectors of B 
               provide the principal component loadings, and scores are computed as the 
               projection of each observation onto these components."),
                          p("Typically, the first two principal components explain >95% of the total 
               variance in sEMG spectral data, reflecting the dominant modes of spectral 
               variation within the population.")
                      ),
                      
                      # --- Section 8: Eigenvector Orientation ---
                      div(class = "doc-section",
                          h5("8. Resolving Eigenvector Sign Ambiguity"),
                          p("Eigenvectors computed from eigendecomposition have an inherent ", 
                            tags$em("sign ambiguity."), "Mathematically, if v is an eigenvector, 
               then -v is equally valid."),
                          p("This application resolves the ambiguity using the method of Bro, Acar & Kolda (2008). 
               For each principal component, the sign is chosen to maximize:"),
                          div(class = "formula", "S = Σ sign(score_j) × score_j²"),
                          p("This weighs each observation's 'vote' by the squared magnitude of its projection, 
               so observations with larger projections have more influence on the orientation.")
                      ),
                      
                      # --- Section 9: Boundary Spectra ---
                      div(class = "doc-section",
                          h5("9. Boundary Spectra and Wavelet Fitting"),
                          p("The PC loadings define a two-dimensional subspace within which all observed 
               spectra approximately lie. The ", tags$em("boundary spectra"), " are the 
               extreme points of this subspace that remain non-negative (a physical requirement 
               for spectral intensities)."),
                          p("These boundary spectra are found by computing linear combinations of PC1 and PC2:"),
                          div(class = "formula", "S = PC1 + a * PC2"),
                          p("where the boundary coefficients a_lo and a_hi are determined by the constraint 
               that S >= 0 for all frequencies. The resulting boundary spectra represent the 
               extreme spectral shapes observed in the dataset—operationally labeled 
               'low-frequency' and 'high-frequency' based on their spectral centroids."),
                          p("Cauchy wavelets are then fitted to these boundary spectra using nonlinear least 
               squares, yielding characteristic center frequencies and bandwidths for each 
               spectral component."),
                          div(class = "doc-note", style = "background: #fff3cd; border-left: 4px solid #856404; padding: 10px; margin: 10px 0;",
                              tags$strong("Note:"), " These labels describe ", tags$em("signal characteristics"), 
                              ", not muscle fiber types or motor unit properties."
                          )
                      ),
                      
                      # --- Section 10: NNLS Decomposition ---
                      div(class = "doc-section",
                          h5("10. Non-Negative Least Squares Decomposition"),
                          p("Each observation's spectrum is decomposed into contributions from the two fitted 
               wavelets using non-negative least squares (NNLS):"),
                          div(class = "formula", "S(f) = c_high * W_high(f) + c_low * W_low(f)"),
                          p("The non-negativity constraint ensures that coefficients remain physically 
               meaningful (spectral intensity cannot be negative). The coefficients are 
               normalized to sum to 1, yielding the relative proportion of spectral content 
               in the high-frequency versus low-frequency ranges."),
                          div(class = "doc-note", style = "background: #e2e3e5; border-left: 4px solid #6c757d; padding: 10px; margin: 10px 0;",
                              tags$strong("Interpretation:"), " The 'High' and 'Low' coefficients describe 
                              where each spectrum falls along the axis of spectral variation identified 
                              by PCA. They quantify ", tags$em("spectral balance"), ", not fiber type 
                              composition or motor unit recruitment."
                          )
                      ),
                      
                      # --- Section 11: Theta Angle ---
                      div(class = "doc-section",
                          h5("11. The Theta Angle"),
                          p("The theta angle provides a single-value summary of each observation's position 
               in the PC1-PC2 plane:"),
                          div(class = "formula", "theta = atan2(PC1_score, PC2_score)"),
                          p("Interpretation of theta (as a descriptive signal metric):"),
                          tags$ul(
                            tags$li(tags$strong("Higher theta values"), " correspond to spectra with 
                      relatively greater low-frequency content."),
                            tags$li(tags$strong("Lower theta values"), " correspond to spectra with 
                      relatively greater high-frequency content."),
                            tags$li(tags$strong("Changes in theta"), " across conditions or over time 
                      indicate shifts in the spectral balance of the recorded signal.")
                          ),
                          div(class = "doc-note", style = "background: #e2e3e5; border-left: 4px solid #6c757d; padding: 10px; margin: 10px 0;",
                              "Theta describes ", tags$em("where a spectrum falls"), " along the primary 
                              axis of variation. It is a descriptive metric of signal characteristics, 
                              not a measure of underlying physiology."
                          )
                      ),
                      
                      # --- Section 12: Quality Metrics ---
                      div(class = "doc-section",
                          h5("12. Quality Metrics"),
                          p("Two variance-explained metrics assess the quality of the decomposition:"),
                          tags$ul(
                            tags$li(tags$strong("S_rec variance:"), " Proportion of total variance explained 
                      by the first two principal components. Values >0.95 indicate that 
                      two-component representation is appropriate."),
                            tags$li(tags$strong("S_gen variance:"), " Proportion of variance explained by 
                      the NNLS reconstruction using fitted wavelets. Should be close to S_rec; 
                      large discrepancies suggest the fitted wavelets poorly approximate the 
                      boundary spectra.")
                          ),
                          p("A successful decomposition typically shows S_rec > 0.95 and |S_rec - S_gen| < 0.05.")
                      ),
                      
                      # --- Section 13: Workflow ---
                      div(class = "doc-section",
                          h5("13. Recommended Workflow"),
                          tags$ol(
                            tags$li(tags$strong("Wavelet Transform Tab:"), 
                                    tags$ul(
                                      tags$li("Load Excel file with one sheet per observation (time in column 1, 
                          signal in column 2)"),
                                      tags$li("Verify sample rate matches your acquisition settings"),
                                      tags$li("OPTIONAL: Apply bandpass filtering (10-500 Hz typical for sEMG)"),
                                      tags$li("Set analysis window offset and length as needed"),
                                      tags$li("Run transform and verify spectra appear reasonable"),
                                      tags$li("Export spectral intensities (Raw recommended for PCA)")
                                    )
                            ),
                            tags$li(tags$strong("PCA & Decomposition Tab:"),
                                    tags$ul(
                                      tags$li("Load exported spectral intensities file"),
                                      tags$li("Run PCA - verify variance explained is adequate (>95%)"),
                                      tags$li("Run wavelet fitting and NNLS decomposition"),
                                      tags$li("Check quality metrics in status panel"),
                                      tags$li("Export results for further statistical analysis")
                                    )
                            )
                          )
                      ),
                      
                      # --- Section 14: Output Files ---
                      div(class = "doc-section",
                          h5("14. Output File Descriptions"),
                          p(tags$strong("Wavelet Transform Tab Export"), " generates an Excel file with:"),
                          tags$table(class = "doc-table",
                                     tags$tr(tags$th("Sheet"), tags$th("Contents")),
                                     tags$tr(
                                       tags$td("Spectral_Intensities"),
                                       tags$td("Matrix of spectral intensities. First column contains center frequencies 
                        (Hz). Each subsequent column represents one observation, with values 
                        corresponding to wavelet intensity at each frequency band. Format 
                        (raw/normalized) depends on export settings.")
                                     ),
                                     tags$tr(
                                       tags$td("Summary"),
                                       tags$td(HTML("Per-observation summary metrics:<br>
                        &bull; <b>Observation</b>: Name identifier<br>
                        &bull; <b>Total_Intensity</b>: Sum of spectral intensities across all 
                        frequency bands (proportional to signal power)<br>
                        &bull; <b>Mean_Frequency_Hz</b>: Intensity-weighted mean frequency, 
                        calculated as &Sigma;(f<sub>c</sub> &times; I) / &Sigma;(I)"))
                                     )
                          ),
                          tags$br(),
                          p(tags$strong("PCA & Decomposition Tab Export"), " generates an Excel file with:"),
                          tags$table(class = "doc-table",
                                     tags$tr(tags$th("Sheet"), tags$th("Contents")),
                                     tags$tr(
                                       tags$td("PC_Loadings"),
                                       tags$td("Principal component loading vectors. Columns: Frequency_Hz, PC_I, PC_II. 
                        These define the spectral shapes that capture the primary modes of 
                        variation in the dataset.")
                                     ),
                                     tags$tr(
                                       tags$td("Theta"),
                                       tags$td(HTML("Per-observation theta angles:<br>
                        &bull; <b>Observation</b>: Name identifier<br>
                        &bull; <b>Condition</b>: Extracted condition label<br>
                        &bull; <b>Theta_rad</b>: Angle in PC space (radians). Higher values 
                        indicate relatively greater low-frequency spectral content. This is a 
                        descriptive signal metric, not a physiological measurement."))
                                     ),
                                     tags$tr(
                                       tags$td("Wavelet_Parameters"),
                                       tags$td(HTML("Fitted Cauchy wavelet parameters for Low and High frequency 
                        components:<br>
                        &bull; <b>fc_Hz</b>: Center frequency<br>
                        &bull; <b>scale_s</b>: Scale parameter<br>
                        &bull; <b>f_low_Hz, f_high_Hz</b>: Band edges at 1/e power<br>
                        &bull; <b>bandwidth_Hz</b>: Spectral bandwidth<br>
                        &bull; <b>dt_ms</b>: Time resolution"))
                                     ),
                                     tags$tr(
                                       tags$td("Coefficients"),
                                       tags$td("NNLS decomposition coefficients for each observation. 'High' and 'Low' 
                        columns represent the relative contribution (0-1) of high-frequency 
                        and low-frequency spectral components. These describe SPECTRAL BALANCE 
                        of the signal, not fiber type composition. Values sum to 1 for each 
                        observation.")
                                     ),
                                     tags$tr(
                                       tags$td("S_reconstructed"),
                                       tags$td("Spectra reconstructed from the first two principal components. 
                        Useful for assessing how well the 2-PC model captures each 
                        observation's spectral shape.")
                                     ),
                                     tags$tr(
                                       tags$td("S_generated"),
                                       tags$td("Spectra generated from NNLS coefficients and fitted wavelets. 
                        Comparison with original spectra indicates decomposition quality.")
                                     ),
                                     tags$tr(
                                       tags$td("Quality_Metrics"),
                                       tags$td("Summary quality indicators: S_rec_variance (variance explained by 
                        2 PCs), S_gen_variance (variance explained by wavelet reconstruction), 
                        and fitted center frequencies.")
                                     )
                          )
                      ),
                      
                      # --- Section 15: Limitations ---
                      div(class = "doc-section",
                          h5("15. Limitations and Considerations"),
                          
                          p(tags$strong("Methodological Considerations:")),
                          tags$ul(
                            tags$li(tags$strong("Cross-subject variability:"), " Electrode placement, subcutaneous 
                      tissue thickness, and muscle geometry differ between individuals. Cross-subject 
                      comparisons should use appropriate statistical modeling and be interpreted 
                      with caution."),
                            tags$li(tags$strong("Volume conductor effects:"), " The same motor unit activity can 
                      produce different spectral signatures depending on electrode position 
                      relative to active fibers."),
                            tags$li(tags$strong("Two-component assumption:"), " Spectral variation may be more 
                      complex than a single high-low axis. Check S_rec variance to verify."),
                            tags$li(tags$strong("Stationarity:"), " Analysis assumes relatively stable muscle 
                      activation within the analysis window."),
                            tags$li(tags$strong("Parameter sensitivity:"), " Default wavelet parameters are 
                      optimized for human locomotion studies. Other applications may benefit 
                      from adjustment.")
                          ),
                          
                          p(tags$strong("Interpretive Considerations:")),
                          tags$ul(
                            tags$li(tags$strong("Signal vs. source:"), " This tool characterizes the ", 
                                    tags$em("signal"), "; spectral features are influenced by multiple 
                      factors beyond motor unit properties alone."),
                            tags$li(tags$strong("Causal attribution:"), " Spectral shifts could reflect 
                      changes in conduction velocity, discharge rate modulation, 
                      technical factors, or combinations thereof."),
                            tags$li(tags$strong("Appropriate framing:"), " Report results as spectral 
                      characteristics rather than direct physiological measurements.")
                          ),
                          
                          p(tags$strong("Recommended Practices:")),
                          tags$ul(
                            tags$li("Use within-subject designs when possible"),
                            tags$li("Report spectral features descriptively"),
                            tags$li("Consider as one component of a multi-method approach"),
                            tags$li("Acknowledge the ongoing scientific debate in publications")
                          )
                      ),
                      
                      # --- Section 16: References ---
                      div(class = "doc-section",
                          h5("16. References"),
                          
                          p(tags$strong("Primary Methodological Source:")),
                          div(class = "reference", 
                              "[1] von Tscharner, V. (2000). Intensity analysis in time-frequency space of surface 
                              myoelectric signals by wavelets of specified resolution. Journal of 
                              Electromyography and Kinesiology, 10(6), 433-445."),
                          
                          p(tags$strong("Reviews on EMG Signal Interpretation:"), style = "margin-top: 15px;"),
                          div(class = "reference",
                              "[2] Farina, D., Merletti, R., & Enoka, R. M. (2014). The extraction of neural 
                              strategies from the surface EMG: an update. Journal of Applied Physiology, 
                              117(11), 1215-1230."),
                          div(class = "reference", 
                              "[3] Farina, D., Merletti, R., & Enoka, R. M. (2004). The extraction of neural 
                              strategies from the surface EMG. Journal of Applied Physiology, 96(4), 1486-1495."),
                          div(class = "reference", 
                              "[4] Enoka, R. M., & Duchateau, J. (2015). Inappropriate interpretation of surface 
                              EMG signals and muscle fiber characteristics impedes understanding of the control 
                              of neuromuscular function. Journal of Applied Physiology, 119(12), 1516-1518."),
                          
                          p(tags$strong("The Point-Counterpoint Debate:"), style = "margin-top: 15px;"),
                          div(class = "reference", 
                              "[5] von Tscharner, V., & Nigg, B. M. (2008). Point: Spectral properties of the surface 
                              EMG can characterize motor unit recruitment strategies and muscle fiber type. 
                              Journal of Applied Physiology, 105(5), 1671-1673."),
                          div(class = "reference", 
                              "[6] Farina, D. (2008). Counterpoint: Spectral properties of the surface EMG do not 
                              provide information about motor unit recruitment and muscle fiber type. 
                              Journal of Applied Physiology, 105(5), 1673-1674."),
                          div(class = "reference", 
                              "[7] von Tscharner, V., & Nigg, B. M. (2008). Last word on point:counterpoint: spectral 
                              properties of the surface EMG can characterize/do not provide information about 
                              motor unit recruitment strategies and muscle fiber type. Journal of Applied Physiology."),
                          
                          p(tags$strong("Additional Methodological References:"), style = "margin-top: 15px;"),
                          div(class = "reference", 
                              "[8] von Tscharner, V., & Goepfert, B. (2006). Estimation of the interplay between 
                              groups of fast and slow muscle fibers of the tibialis anterior and 
                              gastrocnemius muscle while running. J Electromyogr Kinesiol, 16(2), 188-197."),
                          div(class = "reference", 
                              "[9] Bro, R., Acar, E., & Kolda, T. G. (2008). Resolving the sign ambiguity in the 
                              singular value decomposition. Journal of Chemometrics, 22(2), 135-140."),
                          div(class = "reference", 
                              "[10] Wakeling, J. M., Pascual, S. A., Nigg, B. M., & von Tscharner, V. (2001). Surface 
                              EMG shows distinct populations of muscle activity when measured during sustained 
                              sub-maximal exercise. European Journal of Applied Physiology, 86(1), 40-47."),
                          
                          tags$hr(),
                          p(style = "font-size: 11px; color: #6c757d;", 
                            "Author: Brian Benitez | ",
                            tags$a(href = "https://github.com/Brianben95", target = "_blank", 
                                   "github.com/Brianben95"))
                      )
               )
             )
    )
  )
)

# =============================================================================
# SERVER LOGIC
# =============================================================================

server <- function(input, output, session) {
  
  # ---------------------------------------------------------------------------
  # Reactive Values
  # ---------------------------------------------------------------------------
  
  wt <- reactiveValues(
    obs = data.frame(),
    data = list(),
    results = list(),
    current_idx = 1L,
    status = "Ready"
  )
  
  pca <- reactiveValues(
    A = NULL,
    freq = NULL,
    names = NULL,
    eigvec = NULL,
    eigval = NULL,
    scores = NULL,
    theta = NULL,
    fit = NULL,
    S_rec = NULL,
    S_gen = NULL,
    var_rec = NULL,
    var_gen = NULL,
    coeffs = NULL,
    status = "Load spectral data to begin"
  )
  
  # ---------------------------------------------------------------------------
  # Wavelet Transform Tab - Event Handlers
  # ---------------------------------------------------------------------------
  
  # Reset defaults button
  observeEvent(input$wt_reset, {
    updateNumericInput(session, "wt_J", value = DEFAULTS$J)
    updateNumericInput(session, "wt_q", value = DEFAULTS$q)
    updateNumericInput(session, "wt_r", value = DEFAULTS$r)
    updateNumericInput(session, "wt_scale", value = DEFAULTS$scale)
  })
  
  # Load and scan data when file is uploaded
  observeEvent(input$wt_file, {
    req(input$wt_file)
    wt$status <- "Loading..."
    
    tryCatch({
      sheets <- readxl::excel_sheets(input$wt_file$datapath)
      
      obs_list <- list()
      data_list <- list()
      
      for (s in sheets) {
        df <- suppressMessages(readxl::read_excel(
          input$wt_file$datapath, sheet = s, 
          col_names = FALSE, .name_repair = "minimal"
        ))
        names(df) <- c("Time", "Signal")[1:min(2, ncol(df))]
        df <- suppressWarnings(as.data.frame(lapply(df, as.numeric)))
        df <- df[complete.cases(df), ]
        if (nrow(df) > 1) df$Time <- df$Time - df$Time[1]
        
        obs_list[[length(obs_list) + 1]] <- data.frame(
          Name = s,
          Source = basename(input$wt_file$name),
          Points = nrow(df),
          stringsAsFactors = FALSE
        )
        data_list[[s]] <- df
      }
      
      wt$obs <- do.call(rbind, obs_list)
      wt$data <- data_list
      wt$current_idx <- 1L
      wt$results <- list()
      wt$status <- sprintf("Loaded: %d observations", length(sheets))
      
    }, error = function(e) {
      wt$status <- paste("Load error:", e$message)
    })
  })
  
  # Click handler for observation selection
  observeEvent(input$wt_obs_click, {
    req(input$wt_obs_click)
    idx <- input$wt_obs_click$idx
    if (!is.null(idx) && idx >= 1L && idx <= length(wt$data)) {
      wt$current_idx <- idx
    }
  })
  
  # Run transform
  observeEvent(input$wt_run, {
    req(length(wt$data) > 0)
    wt$status <- "Processing..."
    
    tryCatch({
      withProgress(message = "Running wavelet transform", value = 0, {
        n_obs <- length(wt$data)
        wt$results <- list()
        
        for (i in seq_along(wt$data)) {
          incProgress(1/n_obs, detail = names(wt$data)[i])
          
          df <- wt$data[[i]]
          if (nrow(df) < 10) next
          
          sig <- df$Signal
          
          # Remove DC
          if (input$wt_dc) {
            sig <- remove_dc(sig)
          }
          
          # Bandpass filter
          if (input$wt_filter) {
            sig <- bandpass_filter(sig, input$wt_fs, input$wt_flo, input$wt_fhi, input$wt_order)
          }
          
          # Extract window
          start_idx <- min(input$wt_offset + 1, length(sig))
          if (input$wt_npts > 0) {
            end_idx <- min(start_idx + input$wt_npts - 1, length(sig))
          } else {
            end_idx <- length(sig)
          }
          sig_win <- sig[start_idx:end_idx]
          
          if (length(sig_win) < 64) next
          
          # Compute time range for shading
          t_start <- (start_idx - 1) / input$wt_fs
          t_end <- (end_idx - 1) / input$wt_fs
          
          # Run transform
          result <- von_tscharner_transform(
            sig_win, 
            input$wt_fs,
            J = input$wt_J,
            q = input$wt_q,
            r = input$wt_r,
            scale = input$wt_scale,
            renorm = TRUE,  # Always use center-band renormalization
            apply_gauss = TRUE  # Always apply Gaussian equalizer per von Tscharner
          )
          
          # Store segment info for visualization
          result$segment_start <- t_start
          result$segment_end <- t_end
          
          wt$results[[names(wt$data)[i]]] <- result
        }
      })
      
      wt$current_idx <- 1
      wt$status <- sprintf("Complete: %d observations processed", length(wt$results))
      
    }, error = function(e) {
      wt$status <- paste("Error:", e$message)
    })
  })
  
  # ---------------------------------------------------------------------------
  # Wavelet Transform Tab - Outputs
  # ---------------------------------------------------------------------------
  
  output$wt_status <- renderText(wt$status)
  
  output$wt_params <- renderTable({
    # Compute wavelet parameters based on current settings
    J <- input$wt_J
    q <- input$wt_q
    r <- input$wt_r
    scale <- input$wt_scale
    
    # Center frequencies
    fc <- (q + 0:(J-1))^r / scale
    
    # Get empirical values (sigma, dt, bw) for each wavelet
    params <- lapply(1:J, function(j) compute_sigma_empirical(j, J, q, r, scale))
    
    sigma_ms <- sapply(params, function(p) p$sigma)
    dt_ms <- sapply(params, function(p) p$dt)
    bw_hz <- sapply(params, function(p) p$bw)
    
    # Bandwidth × time resolution product
    bw_dt <- bw_hz * dt_ms / 1000
    
    # Create table (0-indexed j to match von Tscharner paper)
    data.frame(
      j = 0:(J-1),
      `fc (Hz)` = round(fc, 2),
      `Δt (ms)` = round(dt_ms, 2),
      `BW (Hz)` = round(bw_hz, 2),
      `BW×Δt` = round(bw_dt, 2),
      `σ (ms)` = round(sigma_ms, 2),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE, width = "100%", align = "c")
  
  # Clickable observation list
  output$wt_obs_list <- renderUI({
    if (nrow(wt$obs) == 0) {
      return(div(style = "padding: 12px; color: #6c757d; text-align: center;",
                 "No observations loaded"))
    }
    
    current <- wt$current_idx
    
    rows <- lapply(seq_len(nrow(wt$obs)), function(i) {
      is_selected <- (i == current)
      div(
        class = paste("obs-row", if (is_selected) "selected" else ""),
        onclick = sprintf("Shiny.setInputValue('wt_obs_click', {idx: %d, ts: Date.now()});", i),
        div(class = "obs-name", wt$obs$Name[i]),
        div(class = "obs-source", wt$obs$Source[i]),
        div(class = "obs-points", paste0(wt$obs$Points[i], " pts"))
      )
    })
    
    tagList(
      div(class = "obs-header",
          div(class = "obs-name", "Name"),
          div(class = "obs-source", "Source"),
          div(class = "obs-points", "Points")
      ),
      rows
    )
  })
  
  # Dynamic header for signal plot
  output$wt_sig_header <- renderText({
    if (length(wt$data) == 0) {
      "Representative Signal"
    } else {
      idx <- as.integer(wt$current_idx)
      if (idx >= 1 && idx <= length(wt$data)) {
        paste0("Signal: ", names(wt$data)[idx])
      } else {
        "Representative Signal"
      }
    }
  })
  
  output$wt_sig <- renderPlotly({
    if (length(wt$data) == 0) return(empty_plot("Load data first"))
    
    idx <- as.integer(wt$current_idx)
    if (idx < 1 || idx > length(wt$data)) return(empty_plot("Invalid observation index"))
    
    df <- wt$data[[idx]]
    if (is.null(df) || nrow(df) < 10) return(empty_plot("No valid data"))
    
    # Get segment info if available (check bounds first)
    res <- NULL
    if (length(wt$results) >= idx) {
      res <- wt$results[[idx]]
    }
    has_segment <- !is.null(res) && !is.null(res$segment_start) && !is.null(res$segment_end)
    
    # Pre-compute hover text
    hover_text <- sprintf("Time: %.3f s<br>Amplitude: %.4f", df$Time, df$Signal)
    
    # Create base plot
    p <- plot_ly(x = df$Time, y = df$Signal, type = "scatter", mode = "lines",
                 line = list(color = COLORS$plot_primary, width = 0.8),
                 hoverinfo = "text",
                 text = hover_text)
    
    # Add shaded region as a shape (doesn't block hover)
    shapes_list <- list()
    if (has_segment) {
      shapes_list[[1]] <- list(
        type = "rect",
        x0 = res$segment_start,
        x1 = res$segment_end,
        y0 = 0,
        y1 = 1,
        yref = "paper",
        fillcolor = "#6c757d",
        opacity = 0.15,
        line = list(color = "#495057", width = 1, dash = "dot")
      )
    }
    
    p %>%
      layout(shapes = shapes_list) %>%
      plot_layout(xlab = "Time (s)", ylab = "Amplitude") %>%
      config(displaylogo = FALSE)
  })
  
  output$wt_spec <- renderPlotly({
    if (length(wt$results) == 0) return(empty_plot("Run transform first"))
    
    idx <- as.integer(wt$current_idx)
    if (idx < 1 || idx > length(wt$results)) return(empty_plot("No results for this observation"))
    
    res <- wt$results[[idx]]
    if (is.null(res)) return(empty_plot("No results for this observation"))
    
    hover_text <- sprintf("fc: %.1f Hz<br>Intensity: %.4f", res$freqc, res$powspec)
    
    plot_ly(x = res$freqc, y = res$powspec, type = "scatter", mode = "lines+markers",
            line = list(color = COLORS$plot_primary, width = 2),
            marker = list(color = COLORS$plot_primary, size = 6),
            hoverinfo = "text",
            text = hover_text) %>%
      plot_layout(xlab = "Center Frequency (Hz)", ylab = "Spectral Intensity") %>%
      config(displaylogo = FALSE)
  })
  
  output$wt_all <- renderPlotly({
    if (length(wt$results) == 0) return(empty_plot("Run transform first"))
    
    p <- plot_ly()
    
    gray_scale <- colorRampPalette(c("#adb5bd", "#212529"))(length(wt$results))
    
    display_mode <- input$wt_display
    
    for (i in seq_along(wt$results)) {
      res <- wt$results[[i]]
      if (is.null(res)) next
      
      # Apply display normalization
      spec <- res$powspec
      if (display_mode == "unit") {
        # Unit area normalization (sum = 1)
        spec <- spec / sum(spec)
      } else if (display_mode == "power") {
        # Power normalization (sum of squares = 1)
        spec <- spec / sqrt(sum(spec^2))
      }
      # "raw" = no normalization
      
      p <- add_trace(p, x = res$freqc, y = spec, 
                     type = "scatter", mode = "lines",
                     name = names(wt$results)[i],
                     line = list(color = gray_scale[i], width = 1),
                     showlegend = FALSE)
    }
    
    ylab <- switch(display_mode,
                   "raw" = "Spectral Intensity",
                   "unit" = "Unit Normalized Intensity",
                   "power" = "Power Normalized Intensity"
    )
    
    p %>%
      plot_layout(xlab = "Center Frequency (Hz)", ylab = ylab) %>%
      config(displaylogo = FALSE)
  })
  
  # Download handler for wavelet results
  output$wt_download <- downloadHandler(
    filename = function() {
      base <- if (nchar(input$wt_fname) > 0) input$wt_fname else "spectral_intensities"
      paste0(base, ".xlsx")
    },
    content = function(file) {
      req(length(wt$results) > 0)
      
      # Build spectral data frame
      fc <- wt$results[[1]]$freqc
      fc_full <- fc  # Keep full frequencies for mean frequency calculation
      if (input$wt_drop && length(fc) > 1) fc <- fc[-1]
      
      df_out <- data.frame(Frequency_Hz = fc)
      
      # Summary metrics storage
      obs_names <- c()
      total_intensities <- c()
      mean_frequencies <- c()
      
      for (name in names(wt$results)) {
        spec_full <- wt$results[[name]]$powspec
        spec <- spec_full
        if (input$wt_drop && length(spec) > 1) spec <- spec[-1]
        
        # Compute summary metrics (using full spectrum or trimmed based on setting)
        fc_calc <- if (input$wt_drop && length(fc_full) > 1) fc_full[-1] else fc_full
        spec_calc <- if (input$wt_drop && length(spec_full) > 1) spec_full[-1] else spec_full
        
        ti <- sum(spec_calc)
        mnf <- if (ti > 0) sum(fc_calc * spec_calc) / ti else NA
        
        obs_names <- c(obs_names, name)
        total_intensities <- c(total_intensities, ti)
        mean_frequencies <- c(mean_frequencies, mnf)
        
        # Apply normalization for spectral output
        if (input$wt_norm == "unit") {
          spec <- spec / sum(spec)
        }
        
        df_out[[name]] <- spec
      }
      
      wb <- openxlsx::createWorkbook()
      
      # Sheet 1: Spectral Intensities
      openxlsx::addWorksheet(wb, "Spectral_Intensities")
      openxlsx::writeData(wb, "Spectral_Intensities", df_out)
      
      # Sheet 2: Summary metrics (if requested)
      if (input$wt_summary) {
        summary_df <- data.frame(
          Observation = obs_names,
          Total_Intensity = total_intensities,
          Mean_Frequency_Hz = round(mean_frequencies, 2)
        )
        openxlsx::addWorksheet(wb, "Summary")
        openxlsx::writeData(wb, "Summary", summary_df)
      }
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
  
  # ---------------------------------------------------------------------------
  # PCA Tab - Event Handlers
  # ---------------------------------------------------------------------------
  
  # Sheet selector for Excel files
  output$pca_sheet_ui <- renderUI({
    req(input$pca_file)
    ext <- tools::file_ext(input$pca_file$name)
    if (ext %in% c("xls", "xlsx")) {
      sheets <- readxl::excel_sheets(input$pca_file$datapath)
      selectInput("pca_sheet", "Sheet", choices = sheets)
    }
  })
  
  # Load spectral data
  observeEvent(input$pca_file, {
    req(input$pca_file)
    
    # Reset all PCA results
    pca$eigvec <- NULL
    pca$eigval <- NULL
    pca$scores <- NULL
    pca$theta <- NULL
    pca$fit <- NULL
    pca$S_rec <- NULL
    pca$S_gen <- NULL
    pca$var_rec <- NULL
    pca$var_gen <- NULL
    pca$coeffs <- NULL
    
    tryCatch({
      ext <- tools::file_ext(input$pca_file$name)
      
      df <- if (ext %in% c("xls", "xlsx")) {
        sheet <- if (!is.null(input$pca_sheet)) input$pca_sheet else 1
        readxl::read_excel(input$pca_file$datapath, sheet = sheet)
      } else if (ext == "csv") {
        read.csv(input$pca_file$datapath)
      } else {
        read.table(input$pca_file$datapath, header = TRUE)
      }
      
      # First column should be frequency
      pca$freq <- as.numeric(df[[1]])
      pca$freq <- pca$freq[!is.na(pca$freq)]
      
      # Extract numeric values from frequency column if needed
      if (any(is.na(pca$freq))) {
        pca$freq <- sapply(df[[1]], extract_numeric)
      }
      
      # Remaining columns are observations
      A_df <- df[, -1, drop = FALSE]
      pca$names <- names(A_df)
      pca$A <- as.matrix(A_df)
      rownames(pca$A) <- NULL
      
      # Ensure matching dimensions
      if (length(pca$freq) != nrow(pca$A)) {
        pca$freq <- pca$freq[1:nrow(pca$A)]
      }
      
      pca$status <- sprintf("Loaded: %d frequencies, %d observations", 
                            length(pca$freq), ncol(pca$A))
      
    }, error = function(e) {
      pca$status <- paste("Load error:", e$message)
    })
  })
  
  # Handle sheet changes for Excel files
  observeEvent(input$pca_sheet, {
    req(input$pca_file, input$pca_sheet)
    
    ext <- tools::file_ext(input$pca_file$name)
    if (!(ext %in% c("xls", "xlsx"))) return()
    
    # Reset all PCA results
    pca$eigvec <- NULL
    pca$eigval <- NULL
    pca$scores <- NULL
    pca$theta <- NULL
    pca$fit <- NULL
    pca$S_rec <- NULL
    pca$S_gen <- NULL
    pca$var_rec <- NULL
    pca$var_gen <- NULL
    pca$coeffs <- NULL
    
    tryCatch({
      df <- readxl::read_excel(input$pca_file$datapath, sheet = input$pca_sheet)
      
      # First column should be frequency
      pca$freq <- as.numeric(df[[1]])
      pca$freq <- pca$freq[!is.na(pca$freq)]
      
      # Extract numeric values from frequency column if needed
      if (any(is.na(pca$freq))) {
        pca$freq <- sapply(df[[1]], extract_numeric)
      }
      
      # Remaining columns are observations
      A_df <- df[, -1, drop = FALSE]
      pca$names <- names(A_df)
      pca$A <- as.matrix(A_df)
      rownames(pca$A) <- NULL
      
      # Ensure matching dimensions
      if (length(pca$freq) != nrow(pca$A)) {
        pca$freq <- pca$freq[1:nrow(pca$A)]
      }
      
      pca$status <- sprintf("Loaded: %d frequencies, %d observations", 
                            length(pca$freq), ncol(pca$A))
      
    }, error = function(e) {
      pca$status <- paste("Load error:", e$message)
    })
  }, ignoreInit = TRUE)
  
  # Run PCA
  observeEvent(input$pca_run, {
    req(pca$A)
    
    tryCatch({
      # Run PCA without centering
      pca_result <- pca_no_centering(pca$A)
      
      # Orient using Bro's method
      oriented <- orient_pcs(
        pca_result$eigvec, 
        pca_result$eigval, 
        pca_result$scores
      )
      
      pca$eigvec <- oriented$eigvec
      pca$eigval <- oriented$eigval
      pca$scores <- oriented$scores
      
      # Reconstruct using first 2 PCs
      pca$S_rec <- pca$eigvec[, 1:2] %*% pca$scores[1:2, ]
      pca$var_rec <- sum(pca$eigval[1:2]) / sum(pca$eigval)
      
      # Compute theta (frequency shift angle)
      pca$theta <- atan2(pca$scores[1, ], pca$scores[2, ])
      
      pca$status <- sprintf("PCA complete: PC1+PC2 explain %.1f%% variance", 
                            pca$var_rec * 100)
      
    }, error = function(e) {
      pca$status <- paste("PCA error:", e$message)
    })
  })
  
  # Run decomposition
  observeEvent(input$pca_decomp, {
    req(pca$eigvec, pca$freq)
    
    tryCatch({
      # Fit wavelets
      pca$fit <- fit_wavelets_to_pcs(pca$eigvec, pca$freq)
      
      # NNLS decomposition
      decomp <- nnls_decomposition(pca$A, pca$fit$wave_hi, pca$fit$wave_lo)
      pca$coeffs <- decomp$coeffs
      colnames(pca$coeffs) <- pca$names
      pca$S_gen <- decomp$S_gen
      pca$var_gen <- decomp$var_gen
      
      # Quality check
      quality_ok <- (pca$var_gen >= 0.80) && 
        (abs(pca$var_rec - pca$var_gen) <= 0.05) &&
        (pca$fit$fc_lo < pca$fit$fc_hi)
      
      pca$status <- sprintf(
        "Decomposition complete\nS_rec: %.1f%%  |  S_gen: %.1f%%\nfc_lo: %.1f Hz  |  fc_hi: %.1f Hz\nStatus: %s",
        pca$var_rec * 100, pca$var_gen * 100,
        pca$fit$fc_lo, pca$fit$fc_hi,
        if (quality_ok) "Pass" else "Check results"
      )
      
    }, error = function(e) {
      pca$status <- paste("Decomposition error:", e$message)
    })
  })
  
  # ---------------------------------------------------------------------------
  # PCA Tab - Outputs
  # ---------------------------------------------------------------------------
  
  output$pca_status <- renderUI({
    tags$pre(pca$status, style = "margin: 0; white-space: pre-wrap;")
  })
  
  output$pca_pc <- renderPlotly({
    if (is.null(pca$eigvec)) return(empty_plot("Run PCA first"))
    
    var1 <- pca$eigval[1] / sum(pca$eigval) * 100
    var2 <- pca$eigval[2] / sum(pca$eigval) * 100
    
    plot_ly() %>%
      add_trace(x = pca$freq, y = pca$eigvec[, 1], type = "scatter", mode = "lines",
                name = sprintf("PC I (%.1f%%)", var1),
                line = list(color = COLORS$pc1, width = 2)) %>%
      add_trace(x = pca$freq, y = pca$eigvec[, 2], type = "scatter", mode = "lines",
                name = sprintf("PC II (%.1f%%)", var2),
                line = list(color = COLORS$pc2, width = 2)) %>%
      plot_layout(xlab = "Frequency (Hz)", ylab = "Loading") %>%
      config(displaylogo = FALSE)
  })
  
  output$pca_mean <- renderPlotly({
    if (is.null(pca$A)) return(empty_plot("Load data first"))
    
    mean_spec <- rowMeans(pca$A)
    
    plot_ly(x = pca$freq, y = mean_spec, type = "scatter", mode = "lines",
            line = list(color = COLORS$plot_primary, width = 2)) %>%
      plot_layout(xlab = "Frequency (Hz)", ylab = "Mean Intensity") %>%
      config(displaylogo = FALSE)
  })
  
  output$pca_wave <- renderPlotly({
    if (is.null(pca$fit)) return(empty_plot("Run decomposition first"))
    
    # Smooth frequency axis for fitted wavelets
    f <- seq(min(pca$freq), max(pca$freq), length.out = 500)
    
    # Compute UNNORMALIZED wavelet amplitudes for display
    y_lo <- sapply(f, function(x) cauchy_amplitude(x, pca$fit$fc_lo, pca$fit$s_lo))
    y_hi <- sapply(f, function(x) cauchy_amplitude(x, pca$fit$fc_hi, pca$fit$s_hi))
    
    plot_ly() %>%
      add_trace(x = f, y = y_lo, type = "scatter", mode = "lines",
                name = sprintf("Low (%.0f Hz)", pca$fit$fc_lo),
                line = list(color = COLORS$low_freq, width = 3)) %>%
      add_trace(x = f, y = y_hi, type = "scatter", mode = "lines",
                name = sprintf("High (%.0f Hz)", pca$fit$fc_hi),
                line = list(color = COLORS$high_freq, width = 3)) %>%
      add_trace(x = pca$freq, y = pca$fit$i_low, type = "scatter", mode = "lines",
                name = "Template (low)",
                line = list(color = COLORS$low_freq, width = 1.5, dash = "dash"),
                showlegend = FALSE) %>%
      add_trace(x = pca$freq, y = pca$fit$i_high, type = "scatter", mode = "lines",
                name = "Template (high)",
                line = list(color = COLORS$high_freq, width = 1.5, dash = "dash"),
                showlegend = FALSE) %>%
      plot_layout(xlab = "Frequency (Hz)", ylab = "Intensity (a.u.)") %>%
      config(displaylogo = FALSE)
  })
  
  output$pca_biplot <- renderPlotly({
    if (is.null(pca$scores) || is.null(pca$names)) return(empty_plot("Run PCA first"))
    
    conditions <- sapply(pca$names, extract_condition)
    cond_data <- list()
    
    for (i in seq_along(pca$names)) {
      cond <- conditions[i]
      if (is.null(cond_data[[cond]])) {
        cond_data[[cond]] <- list(pc1 = c(), pc2 = c())
      }
      cond_data[[cond]]$pc1 <- c(cond_data[[cond]]$pc1, pca$scores[1, i])
      cond_data[[cond]]$pc2 <- c(cond_data[[cond]]$pc2, pca$scores[2, i])
    }
    
    p <- plot_ly()
    gray_scale <- colorRampPalette(c("#6c757d", "#212529"))(length(cond_data))
    cidx <- 1
    
    for (cond in names(cond_data)) {
      mean_pc1 <- mean(cond_data[[cond]]$pc1)
      mean_pc2 <- mean(cond_data[[cond]]$pc2)
      col <- gray_scale[cidx]
      
      p <- add_trace(p, x = c(0, mean_pc2), y = c(0, mean_pc1),
                     type = "scatter", mode = "lines+markers",
                     name = cond,
                     line = list(color = col, width = 2),
                     marker = list(size = 8, color = col))
      cidx <- cidx + 1
    }
    
    var1 <- pca$eigval[1] / sum(pca$eigval) * 100
    var2 <- pca$eigval[2] / sum(pca$eigval) * 100
    
    p %>%
      plot_layout(
        xlab = sprintf("PC II (%.1f%%)", var2), 
        ylab = sprintf("PC I (%.1f%%)", var1)
      ) %>%
      layout(
        xaxis = list(zeroline = TRUE, zerolinecolor = COLORS$border_medium),
        yaxis = list(zeroline = TRUE, zerolinecolor = COLORS$border_medium)
      ) %>%
      config(displaylogo = FALSE)
  })
  
  # Download handler for PCA results
  output$pca_download <- downloadHandler(
    filename = function() {
      base <- if (nchar(input$pca_fname) > 0) input$pca_fname else "pca_decomposition"
      paste0(base, ".xlsx")
    },
    content = function(file) {
      req(pca$freq)
      
      wb <- openxlsx::createWorkbook()
      
      # PC Loadings
      if (!is.null(pca$eigvec)) {
        openxlsx::addWorksheet(wb, "PC_Loadings")
        openxlsx::writeData(wb, "PC_Loadings", data.frame(
          Frequency_Hz = pca$freq,
          PC_I = pca$eigvec[, 1],
          PC_II = pca$eigvec[, 2]
        ))
      }
      
      # Theta values
      if (!is.null(pca$theta)) {
        openxlsx::addWorksheet(wb, "Theta")
        openxlsx::writeData(wb, "Theta", data.frame(
          Observation = pca$names,
          Condition = sapply(pca$names, extract_condition),
          Theta_rad = pca$theta
        ))
      }
      
      # Wavelet parameters
      if (!is.null(pca$fit)) {
        openxlsx::addWorksheet(wb, "Wavelet_Parameters")
        openxlsx::writeData(wb, "Wavelet_Parameters", data.frame(
          Band = c("Low", "High"),
          fc_Hz = c(pca$fit$fc_lo, pca$fit$fc_hi),
          scale_s = c(pca$fit$s_lo, pca$fit$s_hi),
          f_low_Hz = c(pca$fit$flo_lo, pca$fit$flo_hi),
          f_high_Hz = c(pca$fit$fhi_lo, pca$fit$fhi_hi),
          bandwidth_Hz = c(pca$fit$bw_lo, pca$fit$bw_hi),
          dt_ms = c(pca$fit$dt_lo * 1000, pca$fit$dt_hi * 1000)
        ))
      }
      
      # Coefficients
      if (!is.null(pca$coeffs)) {
        openxlsx::addWorksheet(wb, "Coefficients")
        df <- as.data.frame(t(pca$coeffs))
        df$Observation <- colnames(pca$coeffs)
        df$Condition <- sapply(colnames(pca$coeffs), extract_condition)
        df <- df[, c("Observation", "Condition", "High", "Low")]
        openxlsx::writeData(wb, "Coefficients", df)
      }
      
      # Reconstructed spectra
      if (!is.null(pca$S_rec)) {
        openxlsx::addWorksheet(wb, "S_reconstructed")
        df <- data.frame(Frequency_Hz = pca$freq)
        for (i in seq_along(pca$names)) {
          df[[pca$names[i]]] <- pca$S_rec[, i]
        }
        openxlsx::writeData(wb, "S_reconstructed", df)
      }
      
      # Generated spectra
      if (!is.null(pca$S_gen)) {
        openxlsx::addWorksheet(wb, "S_generated")
        df <- data.frame(Frequency_Hz = pca$freq)
        for (i in seq_along(pca$names)) {
          df[[pca$names[i]]] <- pca$S_gen[, i]
        }
        openxlsx::writeData(wb, "S_generated", df)
      }
      
      # Quality metrics
      openxlsx::addWorksheet(wb, "Quality_Metrics")
      openxlsx::writeData(wb, "Quality_Metrics", data.frame(
        Metric = c("S_rec_variance", "S_gen_variance", "fc_low_Hz", "fc_high_Hz"),
        Value = c(
          if (!is.null(pca$var_rec)) pca$var_rec else NA,
          if (!is.null(pca$var_gen)) pca$var_gen else NA,
          if (!is.null(pca$fit)) pca$fit$fc_lo else NA,
          if (!is.null(pca$fit)) pca$fit$fc_hi else NA
        )
      ))
      
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

# =============================================================================
# RUN APPLICATION
# =============================================================================

shinyApp(ui = ui, server = server)