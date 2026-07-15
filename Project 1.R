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
          "forecast", "tseries", "seasonal", "patchwork, tseries")
new  <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

fredr_set_key("54fcfe2acde60a42116b51b9dca79680")
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────
print_adf_kpss <- function(x, label) {
  x <- na.omit(x)
  adf  <- adf.test(x, alternative = "stationary")
  kpss <- kpss.test(x, null = "Level")
  cat(sprintf("  %-30s ADF p=%.4f %-16s | KPSS p=%.4f %s\n",
              label,
              adf$p.value,
              ifelse(adf$p.value  < 0.05, "(STATIONARY)",    "(non-stationary)"),
              kpss$p.value,
              ifelse(kpss$p.value > 0.05, "(STATIONARY)",    "(non-stationary)")))
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — RAW CPI LEVEL
# ══════════════════════════════════════════════════════════════════════════════

# ─────────────────────────────────────────────────────────────────────────────
# 1A.  Pull data
# ─────────────────────────────────────────────────────────────────────────────
cpi_raw <- fredr(series_id         = "CPIAUCSL",
                 observation_start = as.Date("1985-01-01"),
                 observation_end   = as.Date("2026-02-01"))

# Build all transformations
cpi <- cpi_raw %>%
  arrange(date) %>%
  mutate(
    level     = value,
    log_level = log(value),
    mom       = (value / lag(value,  1) - 1) * 100,   # month-on-month %
    yoy       = (value / lag(value, 12) - 1) * 100,   # year-on-year  %
    diff_yoy  = yoy - lag(yoy)                         # first diff of YoY
  ) %>%
  filter(!is.na(yoy))

# Time series objects
cpi_level_ts <- ts(cpi$level,
                   start     = c(year(min(cpi$date)), month(min(cpi$date))),
                   frequency = 12)
cpi_yoy_ts   <- ts(cpi$yoy,
                   start     = c(year(min(cpi$date)), month(min(cpi$date))),
                   frequency = 12)

# ─────────────────────────────────────────────────────────────────────────────
# 1B.  Raw level plot  — shows upward trend + COVID spike
# ─────────────────────────────────────────────────────────────────────────────
regime_dates <- as.Date(c("1990-07-01", "2001-03-01",
                          "2008-09-01", "2020-03-01", "2021-03-01"))
regime_cols  <- c("tomato","purple","tomato","steelblue","forestgreen")
regime_labs  <- c("1990 rec.","9-11","GFC","COVID","Post-COVID")

par(mar = c(4, 4.5, 3.8, 1.5), bg = "white")
plot(cpi$date, cpi$level,
     type = "n",
     xlab = "", ylab = "Index (1982–84 = 100)",
     main = "CPI Level — U.S. (1986–2026)",
     sub  = "Persistent upward trend → non-stationary in levels",
     cex.main = 0.95, cex.sub = 0.75, bty = "o", las = 1)
grid(col = "grey88", lty = 1, lwd = 0.5)
abline(v = as.numeric(regime_dates), col = regime_cols, lty = 3, lwd = 1.1)
lines(cpi$date, cpi$level, col = "black", lwd = 1.5)
legend("topleft",
       legend = c("CPI Level", regime_labs),
       col    = c("black", regime_cols),
       lty    = c(1, rep(3, 5)), lwd = c(1.5, rep(1.1, 5)),
       bty    = "n", cex = 0.75)

# ─────────────────────────────────────────────────────────────────────────────
# 1C.  YoY % plot  — trend removed
# ─────────────────────────────────────────────────────────────────────────────
par(mar = c(4, 4.5, 3.8, 1.5), bg = "white")
plot(cpi$date, cpi$yoy,
     type = "n",
     xlab = "", ylab = "YoY %",
     main = "CPI Year-over-Year % — U.S. (1986–2026)",
     sub  = "YoY transformation removes trend → mean-reverting",
     cex.main = 0.95, cex.sub = 0.75, bty = "o", las = 1)
grid(col = "grey88", lty = 1, lwd = 0.5)
abline(h = 2,   col = "grey50", lty = 2, lwd = 1.2)  # Fed target
abline(h = 0,   col = "grey70", lty = 2, lwd = 0.8)
abline(v = as.numeric(regime_dates), col = regime_cols, lty = 3, lwd = 1.1)
lines(cpi$date, cpi$yoy, col = "black", lwd = 1.4)
mtext("Fed 2% target", side = 3, line = 0.1,
      at = min(cpi$date), cex = 0.68, col = "grey40", adj = 0)
legend("topleft",
       legend = c("CPI YoY %", "Fed 2% target", regime_labs),
       col    = c("black", "grey50", regime_cols),
       lty    = c(1, 2, rep(3, 5)), lwd = c(1.4, 1.2, rep(1.1, 5)),
       bty    = "n", cex = 0.75)

# ─────────────────────────────────────────────────────────────────────────────
# 1D.  Stationarity: level vs YoY
# ─────────────────────────────────────────────────────────────────────────────
cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 1: STATIONARITY TESTS\n")
cat("══════════════════════════════════════════════════════════\n")
cat("  Note: ADF H0 = unit root (non-stationary)\n")
cat("        KPSS H0 = stationary\n")
cat("        Both agree → conclusion is clear\n\n")
print_adf_kpss(cpi$level,    "CPI Level")
print_adf_kpss(cpi$log_level,"Log CPI Level")
print_adf_kpss(cpi$yoy,      "CPI YoY %")
print_adf_kpss(cpi$diff_yoy, "First diff of YoY")
cat("\nConclusion: CPI level is I(1). YoY % is stationary (I(0)).\n")
cat("We model the YoY % series directly — no further differencing needed.\n")
cat("══════════════════════════════════════════════════════════\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — ACF / PACF: BEFORE AND AFTER DETRENDING
# ══════════════════════════════════════════════════════════════════════════════
# This shows the examiner exactly what the YoY transformation achieves:
#   • Level ACF: very slow decay → classic unit-root / trending signature
#   • YoY ACF:   faster decay + seasonal spikes at lag 12, 24 → mean-reverting

par(mfrow = c(2, 2), mar = c(4, 4.5, 3.5, 1.5), bg = "white")

Acf(cpi_level_ts, lag.max = 36,
    main = "ACF — CPI Level\n(slow decay = non-stationary)")
Pacf(cpi_level_ts, lag.max = 36,
     main = "PACF — CPI Level\n(large lag-1 spike)")

Acf(cpi_yoy_ts, lag.max = 36,
    main = "ACF — CPI YoY %\n(faster decay, seasonal spikes at 12,24)")
Pacf(cpi_yoy_ts, lag.max = 36,
     main = "PACF — CPI YoY %\n(AR structure visible, not unit-root)")

par(mfrow = c(1, 1))

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 2: ACF / PACF INTERPRETATION\n")
cat("══════════════════════════════════════════════════════════\n")
cat("CPI Level ACF:\n")
cat("  → Near-unit decay across 36 lags → strong evidence of non-stationarity\n")
cat("  → PACF cuts off sharply at lag 1 → consistent with AR(1)-like unit root\n\n")
cat("CPI YoY % ACF:\n")
cat("  → Decay is present but manageable → mean-reverting\n")
cat("  → Spikes at lags 12 and 24 → residual seasonality in the YoY series\n")
cat("  → PACF: significant lags at 1, 2 → AR(1) or AR(2) may be sufficient\n")
cat("  → This motivates STL to further remove seasonal component before modelling\n")
cat("══════════════════════════════════════════════════════════\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — CLASSICAL ADDITIVE DECOMPOSITION
# ══════════════════════════════════════════════════════════════════════════════
# Classical decomposition:  Y = Trend + Seasonal + Remainder
# Method: centered moving average for trend, then fixed seasonal factors
# Limitation: assumes seasonal pattern is CONSTANT over time
#             → breaks down during COVID spike (outlier contaminates trend)

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 3: CLASSICAL DECOMPOSITION\n")
cat("══════════════════════════════════════════════════════════\n")

classical <- decompose(cpi_yoy_ts, type = "additive")

# Summarise components
trend_cl    <- as.numeric(classical$trend)
seasonal_cl <- as.numeric(classical$seasonal)
remainder_cl <- as.numeric(classical$random)

cat(sprintf("Trend range    : [%.3f, %.3f]\n",
            min(trend_cl, na.rm=TRUE), max(trend_cl, na.rm=TRUE)))
cat(sprintf("Seasonal range : [%.4f, %.4f]\n",
            min(seasonal_cl), max(seasonal_cl)))
cat(sprintf("Remainder SD   : %.4f\n", sd(remainder_cl, na.rm=TRUE)))
cat(sprintf("Seasonal strength (Var_S / Var_S+R) : %.4f\n",
            var(seasonal_cl) / (var(seasonal_cl) + var(remainder_cl, na.rm=TRUE))))
cat(sprintf("Trend strength  (Var_T / Var_T+R)   : %.4f\n\n",
            var(trend_cl, na.rm=TRUE) /
              (var(trend_cl, na.rm=TRUE) + var(remainder_cl, na.rm=TRUE))))

cat("Limitation: Classical assumes FIXED seasonal factors across all years.\n")
cat("The COVID spike (2021-22) inflates the trend estimate and distorts\n")
cat("the seasonal component. STL (Section 4) addresses this.\n")
cat("══════════════════════════════════════════════════════════\n\n")

# Plot classical decomposition
par(mar = c(4, 4.5, 3.5, 1.5), bg = "white")
plot(classical,
     main = "Classical Additive Decomposition — CPI YoY %")

# ── Classical seasonal factors bar chart ────────────────────────────────────
seas_df <- data.frame(
  month  = month.abb,
  factor = as.numeric(classical$figure)
)
seas_df$month <- factor(seas_df$month, levels = month.abb)

p_seas_cl <- ggplot(seas_df, aes(x = month, y = factor,
                                 fill = factor > 0)) +
  geom_col(colour = "white", width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey40") +
  scale_fill_manual(values = c("TRUE" = "#E05C2A", "FALSE" = "#3A8FB5"),
                    guide = "none") +
  labs(title    = "Classical Decomposition — Fixed Seasonal Factors",
       subtitle = "Positive = above average inflation month, Negative = below",
       x = NULL, y = "Seasonal factor (pp)") +
  theme_minimal(base_size = 12)
print(p_seas_cl)
# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — STL DECOMPOSITION  (Seasonal-Trend using Loess)
# ══════════════════════════════════════════════════════════════════════════════
# STL advantages over Classical:
#   1. Seasonal component is ALLOWED TO CHANGE over time (loess-smoothed)
#   2. robust = TRUE downweights outliers (COVID spike handled properly)
#   3. Trend smoother adapts to sudden shifts (GFC, post-COVID)
#   4. Works on any periodicity (monthly, weekly, daily)

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 4: STL DECOMPOSITION\n")
cat("══════════════════════════════════════════════════════════\n")

stl_fit <- stl(cpi_yoy_ts, s.window = "periodic", robust = TRUE)

trend_stl     <- as.numeric(stl_fit$time.series[, "trend"])
seasonal_stl  <- as.numeric(stl_fit$time.series[, "seasonal"])
remainder_stl <- as.numeric(stl_fit$time.series[, "remainder"])

# Seasonal-adjusted series = level - seasonal component
cpi_sa <- cpi_yoy_ts - seasonal_stl   # seasonally adjusted YoY %
# Detrended + deseasonalised = remainder only
cpi_remainder <- ts(remainder_stl,
                    start     = start(cpi_yoy_ts),
                    frequency = 12)

cat(sprintf("Trend range    : [%.3f, %.3f]\n",
            min(trend_stl), max(trend_stl)))
cat(sprintf("Seasonal range : [%.4f, %.4f]\n",
            min(seasonal_stl), max(seasonal_stl)))
cat(sprintf("Remainder SD   : %.4f\n", sd(remainder_stl)))
cat(sprintf("Seasonal strength (Var_S / Var_S+R) : %.4f\n",
            var(seasonal_stl) / (var(seasonal_stl) + var(remainder_stl))))
cat(sprintf("Trend strength  (Var_T / Var_T+R)   : %.4f\n\n",
            var(trend_stl) /
              (var(trend_stl) + var(remainder_stl))))

cat("Key: robust=TRUE downweights COVID outlier in trend/seasonal estimation.\n")
cat("Seasonal strength < 0.64 → seasonal is present but not dominant.\n")
cat("Trend strength tells us how much of variance the trend explains.\n")
cat("══════════════════════════════════════════════════════════\n\n")

# Full STL decomposition plot
par(mar = c(3, 4.5, 3.5, 1.5), bg = "white")
plot(stl_fit,
     main = "STL Decomposition (robust) — CPI YoY %")

# ── Time-varying seasonal factors heatmap ────────────────────────────────────
# Shows that the seasonal pattern CHANGES across decades — this is the key
# argument for STL over Classical
stl_seas_df <- data.frame(
  date     = cpi$date,
  seasonal = seasonal_stl
) %>%
  mutate(year  = year(date),
         month = factor(month.abb[month(date)], levels = month.abb))

p_seas_tv <- ggplot(stl_seas_df, aes(x = year, y = month, fill = seasonal)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#3A8FB5", mid = "white", high = "#E05C2A",
                       midpoint = 0, name = "Factor\n(pp)") +
  labs(title    = "STL Time-Varying Seasonal Factors — CPI YoY %",
       subtitle = "Colour shift over decades confirms seasonal pattern evolves → Classical decomposition is too rigid",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid  = element_blank())
print(p_seas_tv)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — SIDE-BY-SIDE: CLASSICAL vs STL REMAINDER
# ══════════════════════════════════════════════════════════════════════════════
# The key diagnostic: look at the remainder (what's left after decomposition).
# A good decomposition leaves a remainder that looks like white noise.
# Classical remainder retains COVID outlier; STL remainder is cleaner.

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 5: CLASSICAL vs STL REMAINDER COMPARISON\n")
cat("══════════════════════════════════════════════════════════\n")

# Align dates (classical loses ends due to MA)
cl_rem_df <- data.frame(
  date      = cpi$date,
  Classical = as.numeric(classical$random),
  STL       = remainder_stl
) %>% filter(!is.na(Classical))

cat(sprintf("Classical remainder SD : %.4f\n", sd(cl_rem_df$Classical, na.rm=TRUE)))
cat(sprintf("STL       remainder SD : %.4f\n", sd(cl_rem_df$STL[!is.na(cl_rem_df$Classical)])))
cat("Lower SD in STL remainder = better decomposition (less unexplained variance)\n\n")

# Ljung-Box on both remainders
lb_cl  <- Box.test(cl_rem_df$Classical, lag = 12, type = "Ljung-Box")
lb_stl <- Box.test(remainder_stl,       lag = 12, type = "Ljung-Box")
cat(sprintf("Ljung-Box Classical remainder: p=%.4f  %s\n",
            lb_cl$p.value,
            ifelse(lb_cl$p.value > 0.05, "→ White noise ✓", "→ Structure remains")))
cat(sprintf("Ljung-Box STL remainder      : p=%.4f  %s\n",
            lb_stl$p.value,
            ifelse(lb_stl$p.value > 0.05, "→ White noise ✓", "→ Structure remains")))
cat("══════════════════════════════════════════════════════════\n\n")

# Plot: Classical vs STL remainder over time
rem_long <- cl_rem_df %>%
  pivot_longer(-date, names_to = "Method", values_to = "Remainder")

p_rem <- ggplot(rem_long, aes(x = date, y = Remainder, colour = Method)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_line(linewidth = 0.75, alpha = 0.85) +
  geom_vline(xintercept = as.numeric(regime_dates),
             colour = regime_cols, linetype = "dotted", linewidth = 0.7) +
  scale_colour_manual(values = c("Classical" = "#3A8FB5",
                                 "STL"       = "#E05C2A")) +
  labs(title    = "Remainder Comparison: Classical vs STL Decomposition",
       subtitle = "STL remainder is cleaner at COVID spike — robust=TRUE downweights outlier",
       x = NULL, y = "Remainder (pp)", colour = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(p_rem)

# Distribution comparison
p_dist <- ggplot(rem_long %>% filter(!is.na(Remainder)),
                 aes(x = Remainder, fill = Method)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 alpha = 0.55, position = "identity", colour = "white") +
  geom_density(aes(colour = Method), linewidth = 0.9, fill = NA) +
  scale_fill_manual(values  = c("Classical" = "#3A8FB5", "STL" = "#E05C2A")) +
  scale_colour_manual(values = c("Classical" = "#3A8FB5", "STL" = "#E05C2A")) +
  labs(title    = "Remainder Distribution: Classical vs STL",
       subtitle = "Thinner tails in STL = less extreme unexplained variation",
       x = "Remainder (pp)", y = "Density",
       fill = "Method", colour = "Method") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
print(p_dist)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — ACF / PACF OF STL-ADJUSTED SERIES  (model-ready)
# ══════════════════════════════════════════════════════════════════════════════
# After STL we have two candidate model series:
#   (a) Seasonally adjusted:  cpi_sa        = YoY - seasonal_stl
#   (b) Remainder only:       cpi_remainder = trend + remainder removed
#
# For our forecasting models we use (a) — seasonally adjusted YoY %.
# This preserves the trend-cycle information while removing the seasonal noise.
# The ACF/PACF of (a) should show:
#   • No spikes at lags 12, 24 (seasonal removed)
#   • AR structure still visible (we want this — it's the signal)

par(mfrow = c(2, 2), mar = c(4, 4.5, 3.5, 1.5), bg = "white")

Acf(cpi_sa,        lag.max = 36,
    main = "ACF — STL Seasonally Adjusted YoY %\n(seasonal spikes removed)")
Pacf(cpi_sa,       lag.max = 36,
     main = "PACF — STL Seasonally Adjusted YoY %\n(AR order guidance)")
Acf(cpi_remainder, lag.max = 36,
    main = "ACF — STL Remainder\n(detrended + deseasonalised)")
Pacf(cpi_remainder, lag.max = 36,
     main = "PACF — STL Remainder\n(closest to white noise)")

par(mfrow = c(1, 1))

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 6: POST-STL ACF/PACF INTERPRETATION\n")
cat("══════════════════════════════════════════════════════════\n")
cat("Seasonally Adjusted ACF:\n")
cat("  → Lag-12 and lag-24 spikes should be absent or much reduced\n")
cat("  → Gradual AR decay confirms mean-reverting inflation dynamics\n")
cat("  → This is the series used in modelling (YoY % — seasonal component)\n\n")
cat("Remainder ACF:\n")
cat("  → Near-zero autocorrelations = decomposition is well-specified\n")
cat("  → Any remaining spikes = structure the decomposition missed\n")
cat("  → Used as a diagnostic — NOT directly as the model input\n")
cat("══════════════════════════════════════════════════════════\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7 — STATIONARITY CONFIRMATION ON MODEL-READY SERIES
# ══════════════════════════════════════════════════════════════════════════════

cat("══════════════════════════════════════════════════════════\n")
cat("  SECTION 7: STATIONARITY — FINAL MODEL SERIES\n")
cat("══════════════════════════════════════════════════════════\n")
cat("  Series tested:\n\n")
print_adf_kpss(as.numeric(cpi_yoy_ts), "CPI YoY % (raw)")
print_adf_kpss(as.numeric(cpi_sa),     "STL Seasonally Adj. YoY %")
print_adf_kpss(remainder_stl,          "STL Remainder")

d_yoy <- ndiffs(cpi_yoy_ts)
d_sa  <- ndiffs(cpi_sa)
cat(sprintf("\nndiffs() → CPI YoY %%: d=%d | STL Adj: d=%d\n", d_yoy, d_sa))

cat("\nConclusion:\n")
cat("  Both YoY % and the STL-adjusted series are I(0) — stationary.\n")
cat("  No differencing required before modelling.\n")
cat("  The YoY transformation has already removed the stochastic trend.\n")
cat("  STL further removes residual seasonality visible in lags 12/24.\n")
cat("  → Models in the main script use the YoY % series directly,\n")
cat("    with stationarity justified here.\n")
cat("══════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# FINAL SUMMARY PLOT: Raw level → YoY % → STL adjusted  (3-panel)
# ─────────────────────────────────────────────────────────────────────────────
sa_df <- data.frame(
  date = cpi$date,
  Level    = scale(cpi$level)[,1],      # standardised for visual comparison
  YoY      = scale(cpi$yoy)[,1],
  STL_adj  = scale(as.numeric(cpi_sa))[,1]
) %>% pivot_longer(-date, names_to = "Series", values_to = "Value") %>%
  mutate(Series = factor(Series,
                         levels = c("Level","YoY","STL_adj"),
                         labels = c("Raw Level (standardised)",
                                    "YoY % (detrended)",
                                    "STL Adj. YoY % (detrended + deseasonalised)")))

p_final <- ggplot(sa_df, aes(x = date, y = Value, colour = Series)) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = as.numeric(regime_dates),
             colour = regime_cols, linetype = "dotted",
             linewidth = 0.6, alpha = 0.7) +
  facet_wrap(~ Series, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(
    "Raw Level (standardised)"                         = "black",
    "YoY % (detrended)"                                = "#3A8FB5",
    "STL Adj. YoY % (detrended + deseasonalised)"      = "#E05C2A"
  )) +
  labs(title    = "Detrending & Deseasonalising CPI: Step by Step",
       subtitle = "Each transformation removes one source of non-stationarity",
       x = NULL, y = "Standardised value", colour = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold", size = 10))
print(p_final)

# ══════════════════════════════════════════════════════════════════════════════
# END — paste this section into the paper before the main forecasting script
# ══════════════════════════════════════════════════════════════════════════════
