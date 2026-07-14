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


