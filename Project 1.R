# PURPOSE:
#   Before modelling we need to verify:
#     (1) The raw CPI level is non-stationary (trending upward)
#     (2) The YoY % transformation removes the trend
#     (3) There is seasonal structure worth acknowledging
#     (4) STL is more robust than Classical decomposition for CPI
#         because of the COVID outlier spike (2021–2022)
#
# STRUCTURE:
#   Section 1 — Raw CPI level: plot + ADF + KPSS
#   Section 2 — ACF / PACF: raw level vs YoY %  (before/after detrending)
#   Section 3 — Classical additive decomposition
#   Section 4 — STL decomposition (robust to outliers)
#   Section 5 — Side-by-side comparison: Classical vs STL remainder
#   Section 6 — ACF / PACF of STL-adjusted series (model-ready)
#   Section 7 — Stationarity confirmation on final model series
# ══════════════════════════════════════════════════════════════════════════════


# ─────────────────────────────────────────────────────────────────────────────
# 0.  PACKAGES
# ─────────────────────────────────────────────────────────────────────────────
pkgs <- c("fredr", "dplyr", "ggplot2", "lubridate",
          "forecast", "tseries", "seasonal", "patchwork")
new  <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

fredr_set_key("54fcfe2acde60a42116b51b9dca79680")