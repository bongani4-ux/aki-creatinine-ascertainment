# =============================================================================
# pipeline.R
# Core analysis for: What can be measured determines what can be found.
#
# One pass over each source's creatinine series computes every per-source
# quantity reported in the manuscript:
#   - measurement intensity (results per patient, never measured)
#   - ascertainment cascade under three baseline rules
#   - baseline provenance by consensus tier
#   - undefined vs unmeasured decomposition (where death timing exists)
#   - AKI classification under hierarchy vs proximal baseline (same patients)
#   - hour-of-day and hours-from-t0 distributions (where clock time exists)
#
# Inputs follow the common export specification (docs/data_specification.md):
#   <SOURCE>_creatinine_series.csv   one row per creatinine result
#   <SOURCE>_episode_covariates.csv  one row per episode
# =============================================================================

suppressPackageStartupMessages(library(data.table))

TIERS <- data.table(lo = c(7, 30, 90, 180), hi = c(30, 90, 180, 365),
                    tier = c("T1", "T2", "T3", "T4"), rank = 1:4)

# ---------------------------------------------------------------- helpers
num <- function(x) suppressWarnings(as.numeric(x))

parse_time <- function(x) {
  # Element-wise parser. as.POSIXct(tryFormats=) is NOT used: it returns NA for
  # the whole vector if any single value (e.g. an empty string) fails a format.
  # Accepts "2150-03-02T15:19:31Z", "2023-06-10T07:12:00", "2019-09-16 08:45:00".
  x <- trimws(as.character(x))
  x[x %in% c("", "nan", "NaN", "NA", "NaT", "None", "null")] <- NA
  x <- sub("\\+00:?00$", "", sub("Z$", "", sub("T", " ", x)))
  out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%d %H:%M:%OS")
  for (f in c("%Y-%m-%d %H:%M", "%Y-%m-%d")) {
    miss <- is.na(out) & !is.na(x)
    if (any(miss)) out[miss] <- as.POSIXct(x[miss], tz = "UTC", format = f)
  }
  out
}

days_between <- function(a, b) as.numeric(difftime(a, b, units = "days"))

first_col <- function(dt, names) {
  hit <- intersect(names, names(dt))
  if (length(hit)) num(dt[[hit[1]]]) else rep(NA_real_, nrow(dt))
}

# ---------------------------------------------------------------- episodes
#' Keep one episode per patient.
#'
#' The earliest episode in time is chosen where t0 carries real time (a
#' datetime, or a per-patient offset that varies between episodes). Where every
#' episode has the same t0 (e.g. all 0), time order cannot be recovered from the
#' export and episode_id order is used instead.
first_episode_per_patient <- function(cov) {
  cov <- unique(cov, by = "episode_id")
  tcol <- intersect(c("t0", "t0_min"), names(cov))
  rule <- "episode_id"
  cov[, .ord := NA_real_]
  if (length(tcol)) {
    raw <- cov[[tcol[1]]]
    n <- num(raw)
    if (mean(!is.na(n)) > 0.9) {
      if (uniqueN(n[!is.na(n)]) > 1) { cov[, .ord := n]; rule <- "t0 (numeric offset)" }
    } else {
      dt <- parse_time(raw)
      if (mean(!is.na(dt)) > 0.9) { cov[, .ord := as.numeric(dt)]; rule <- "t0 (datetime)" }
    }
  }
  # C-locale radix ordering of IDs; missing times sort last (as in the paper)
  setorderv(cov, c("patient_id", ".ord", "episode_id"), na.last = TRUE)
  cov <- cov[!duplicated(patient_id)]
  cov[, .ord := NULL]
  list(cov = cov, rule = rule)
}

PER_DAY <- c(seconds = 86400, minutes = 1440, hours = 24, milliseconds = 86400000)

#' Discharge and death in days from t0.
#' Uses *_t_days where supplied. Otherwise derives from discharge_time /
#' death_time: as numeric offsets in `offset_unit` when the source records
#' times that way (eICU-CRD, SICdb, VitalDB), or as datetimes relative to t0.
offsets <- function(cov, offset_unit = NA) {
  disch <- first_col(cov, "discharge_t_days")
  death <- first_col(cov, "death_t_days")
  derive <- function(col) {
    if (!col %in% names(cov)) return(rep(NA_real_, nrow(cov)))
    if (!is.na(offset_unit)) {
      t0 <- if ("t0" %in% names(cov)) num(cov$t0) else 0
      t0[is.na(t0)] <- 0
      return((num(cov[[col]]) - t0) / PER_DAY[[offset_unit]])
    }
    t0 <- if ("t0" %in% names(cov)) parse_time(cov$t0) else rep(as.POSIXct(NA), nrow(cov))
    days_between(parse_time(cov[[col]]), t0)
  }
  if (all(is.na(disch))) disch <- derive("discharge_time")
  if (all(is.na(death))) death <- derive("death_time")
  if (all(is.na(disch))) disch <- first_col(cov, "hosp_los_days")
  list(disch = disch, death = death)
}

# ---------------------------------------------------------------- staging
kdigo_stage <- function(base, peak7, peak48) {
  ratio <- peak7 / base
  rise  <- peak48 - base
  s <- integer(length(base))
  s[(ratio >= 1.5) | (rise >= 0.3)] <- 1L
  s[ratio >= 2.0] <- 2L
  s[(ratio >= 3.0) | (peak7 >= 4.0)] <- 3L
  s
}

nz <- function(x) { x[is.na(x)] <- FALSE; x }   # NA -> FALSE

# ---------------------------------------------------------------- main
analyse_source <- function(name, series_path, covar_path, setting, clock_col = NA,
                           offset_unit = NA) {
  cov <- fread(covar_path, colClasses = "character", showProgress = FALSE)
  setnames(cov, trimws(names(cov)))
  if (!"patient_id" %in% names(cov)) cov[, patient_id := episode_id]
  n_episodes <- uniqueN(cov$episode_id)
  fe <- first_episode_per_patient(cov); cov <- fe$cov
  ids <- cov$episode_id; n <- length(ids)

  sel <- c("episode_id", "t_days", "value", if (!is.na(clock_col)) clock_col)
  s <- fread(series_path, select = sel, colClasses = list(character = c("episode_id",
             if (!is.na(clock_col)) clock_col)), showProgress = FALSE)
  s[, `:=`(t_days = num(t_days), value = num(value))]
  s <- s[episode_id %chin% ids & !is.na(t_days) & !is.na(value)]
  total <- nrow(s)

  has <- function(cond) unique(s$episode_id[cond])
  t <- s$t_days
  pre <- has(t < 0); perm <- has(t <= 1); hier_any <- has(t <= -7 & t >= -365)
  aki_w <- has(t > 0 & t <= 7); akd_w <- has(t > 7 & t <= 90)

  # baseline provenance: best consensus tier, else tier 0 if any value in (0, 7] d before
  lead <- -t
  s[, rk := NA_integer_]
  for (i in seq_len(nrow(TIERS))) s[lead > TIERS$lo[i] & lead <= TIERS$hi[i], rk := TIERS$rank[i]]
  s[is.na(rk) & lead > 0 & lead <= 7, rk := 5L]
  best <- if (any(!is.na(s$rk))) s[!is.na(rk), .(rk = min(rk)), by = episode_id] else data.table(episode_id = character(), rk = integer())

  # most-recent hierarchy and proximal baseline values; post-index peaks
  hb <- s[t_days <= -7 & t_days >= -365][order(episode_id, t_days)][, .SD[.N], by = episode_id]
  pb <- s[t_days < 0 & t_days > -7][order(episode_id, t_days)][, .SD[.N], by = episode_id]
  p7 <- s[t_days > 0 & t_days <= 7, .(p7 = max(value)), by = episode_id]
  p48 <- s[t_days > 0 & t_days <= 2, .(p48 = max(value)), by = episode_id]

  # hourly timing
  h <- t * 24; h <- h[h >= -24 & h <= 168]
  hft <- tabulate(pmin(floor(h) + 24L, 191L) + 1L, nbins = 192)
  hod <- NULL
  if (!is.na(clock_col) && clock_col %in% names(s)) {
    hr <- as.integer(format(parse_time(s[[clock_col]]), "%H"))
    hod <- tabulate(hr[!is.na(hr)] + 1L, nbins = 24)
    if (sum(hod) == 0) hod <- NULL
  }

  cnt <- s[, .N, by = episode_id][data.table(episode_id = ids), on = "episode_id"]$N
  cnt[is.na(cnt)] <- 0L
  L1 <- ids %chin% pre; L2 <- L1 & ids %chin% aki_w; L3 <- L2 & ids %chin% akd_w
  rk <- best[data.table(episode_id = ids), on = "episode_id"]$rk
  tv <- c("T1", "T2", "T3", "T4", "T0")[rk]; tv[is.na(tv)] <- "None"

  # undefined vs unmeasured (AKD window)
  o <- offsets(cov, offset_unit); disch <- o$disch; death <- o$death
  non3 <- !L3
  died7 <- sum(nz(non3 & death <= 7))
  died790 <- sum(nz(non3 & death > 7 & death <= 90))
  alive7 <- !nz(death <= 7)
  dis90 <- sum(nz(non3 & alive7 & !nz(death > 7 & death <= 90) & nz(disch <= 90)))
  other <- max(sum(non3) - died7 - died790 - dis90, 0)
  elig <- sum(L3) + dis90 + other

  # undefined vs unmeasured (AKI window): death truncates rather than prevents
  non2 <- !L2
  d48 <- sum(nz(non2 & death <= 2)); d2_7 <- sum(nz(non2 & death > 2 & death <= 7))
  dis7 <- sum(nz(non2 & !nz(death <= 7) & nz(disch <= 7)))
  un2 <- max(sum(non2) - d48 - d2_7 - dis7, 0); elig2 <- sum(L2) + dis7 + un2

  q <- quantile(cnt, c(.25, .5, .75, .95), type = 7, names = FALSE)
  pct <- function(x) 100 * mean(x)
  measure <- list(
    db = name, setting = setting, first_episode_rule = fe$rule, episodes = n_episodes,
    n = n, results = total, med = as.integer(q[2]), q1 = as.integer(q[1]),
    q3 = as.integer(q[3]), p95 = as.integer(q[4]), mx = max(cnt),
    zero = pct(cnt == 0), L1 = pct(L1), L2 = pct(L2), L3 = pct(L3),
    hier = pct(tv %in% c("T1", "T2", "T3", "T4")),
    T1 = pct(tv == "T1"), T2 = pct(tv == "T2"), T3 = pct(tv == "T3"), T4 = pct(tv == "T4"),
    t0 = pct(tv == "T0"), nopre = pct(tv == "None"),
    died7 = died7, died790 = died790, dis90 = dis90, inh = other, L3n = sum(L3),
    akd_elig = if (elig) 100 * sum(L3) / elig else NA_real_,
    death_avail = pct(!is.na(death)),
    aki_died48 = d48, aki_died2_7 = d2_7, aki_disch7 = dis7, aki_unmeasured = un2,
    aki_elig = if (elig2) 100 * sum(L2) / elig2 else NA_real_)

  rules <- list(db = name, n = n)
  for (r in list(list("perm", perm), list("strict", pre), list("hier", hier_any))) {
    a <- ids %chin% r[[2]]; b <- a & ids %chin% aki_w; cc <- b & ids %chin% akd_w
    rules[[paste0(r[[1]], "_L1")]] <- pct(a); rules[[paste0(r[[1]], "_L2")]] <- pct(b)
    rules[[paste0(r[[1]], "_L3")]] <- pct(cc)
  }

  misclass <- NULL
  both <- Reduce(merge, list(hb[, .(episode_id, hb = value)], pb[, .(episode_id, pb = value)], p7))
  both <- merge(both, p48, by = "episode_id", all.x = TRUE)
  if (nrow(both) >= 100) {
    both[is.na(p48), p48 := p7]
    sh <- kdigo_stage(both$hb, both$p7, both$p48) > 0
    sp <- kdigo_stage(both$pb, both$p7, both$p48) > 0
    misclass <- list(db = name, n = nrow(both), hier_med = median(both$hb), t0_med = median(both$pb),
                     diff_med = median(both$pb - both$hb), aki_hier = pct(sh), aki_t0 = pct(sp),
                     agree = pct(sh == sp), missed = pct(sh & !sp), false = pct(!sh & sp))
  }

  list(measure = measure, rules = rules, misclass = misclass,
       timing = list(db = name, hours_from_t0 = hft, hour_of_day = hod))
}
