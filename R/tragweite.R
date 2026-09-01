#' Detection probability of a study design (simulation-based)
#'
#' The core question of tragweite: *can this design see the effect it
#' is looking for?* Simulates `reps` worlds with the requested effect
#' built in (via the windkanal generator), analyzes each world with a
#' validated, real-data-computable criterion, and returns the share of
#' worlds in which the effect was detected -- the design's detection
#' probability (power), measured rather than derived from a formula.
#' For the questions tragweite targets, no closed-form power formulas
#' exist.
#'
#' Detection criteria (Satterthwaite mixed-model tests, validated
#' against the windkanal simulation engine):
#' * `"ate"`: test of the treatment coefficient in
#'   `score ~ z + (1|therapist) + (1|patient)`
#' * `"moderation"`: test of the `z:x` interaction
#' * `"matching"`: test of the `z:x:c` triple interaction (patient
#'   attribute x therapist attribute x treatment)
#'
#' **Null calibration of the matching criterion:** with
#' `tau_c = 0` and `tau_xc = 0` the generator draws no therapist
#' attribute (it would not act). In the field, C is *measured*
#' whether or not it acts -- so for `effect = "matching"` tragweite
#' attaches an inert standard-normal therapist attribute to such
#' worlds (own seed stream `seed + i + 1e7`, decoupled from the
#' world's draws). The triple test then runs against a measured,
#' truly effect-free attribute, which is exactly what a null world
#' looks like from the analyst's chair.
#'
#' Failed model fits do not abort the run; they are counted and
#' reported (`n_failed`, with the first error message as a warning
#' when any occur). If *every* world fails, the run stops with that
#' message instead of returning a silent `NaN`.
#'
#' @section The v2S detector is not available yet:
#' `detector = "duet_s"` applies an amplitude rule to the posterior of
#' a dyadic sampler, `windkanal::fit_cate_dyade_v2s()`. That sampler is
#' **not part of the public windkanal 0.3.0**; it is planned for
#' windkanal 0.4. Until then `detector = "duet_s"` (and
#' [ceiling_certificate()], which is built on it) stop with an
#' explanatory error at the point of call. The default
#' `detector = "lmer"` needs only `windkanal::sim_stream()` and runs
#' against the public windkanal 0.3.0 in full.
#'
#' @param effect One of `"ate"`, `"moderation"`, `"matching"`.
#' @param tau True average treatment effect (default 0.5).
#' @param tau_x True patient-moderation amplitude (default 0; used
#'   for `"moderation"` and `"matching"`).
#' @param tau_c True therapist-attribute effect (default 0; only
#'   `"matching"`).
#' @param tau_xc True matching amplitude (default 0; only
#'   `"matching"`).
#' @param n_therapists,caseload,n_sessions,icc,z_level Design under
#'   evaluation; passed to [windkanal::sim_stream()] (`caseload` =
#'   `patients_per_therapist`).
#' @param tau_xc_form,window_delta Functional form of the dyadic term,
#'   passed to the generator (see windkanal; `window_delta` is passed
#'   as the generator's `fenster_delta`). Default is the product form.
#' @param rel_x,rel_c Reliability of the OBSERVED patient / therapist
#'   attribute (1 = perfectly measured). Below 1 the analysis sees an
#'   attenuated version (squared correlation with the truth = rel)
#'   while the truth in the generator is untouched -- measurement error
#'   on both sides, which attenuates a matching effect twice.
#' @param detector `"lmer"` (default; Satterthwaite mixed model) or
#'   `"duet_s"` (dyadic-sampler amplitude rule; matching only, and
#'   currently unavailable -- see the section above). The shipped
#'   amplitude rule is calibrated on its home world stream and is
#'   stream-sensitive; for a planning statement use
#'   [ceiling_certificate()], which always recalibrates the cut on the
#'   stream it is given. The sampler behind `"duet_s"` runs on a
#'   constant seed (`.duet_p_sd(mcmc_seed = 1)`), deliberately: the
#'   sampler stream is held fixed across worlds so that a power
#'   difference between two designs is a design difference and not
#'   sampler noise. The user seed therefore moves the worlds, not the
#'   sampler.
#' @param duet_rule Named vector `c(delta_sd, cut_sd)` of the frozen
#'   v2S amplitude rule (strict cut convention).
#' @param reps Number of simulated worlds (default 200; must be >= 1).
#' @param alpha Test level of the detection criterion (default .05;
#'   must be a single number in (0, 1)). Acts on `detector = "lmer"`
#'   only -- the amplitude rule decides via `duet_rule`, so for
#'   `detector = "duet_s"` the returned `alpha` is `NA`.
#' @param seed Base random seed (mandatory; world `i` uses
#'   `seed + i`). The global RNG state is restored on exit, so a
#'   `tragweite()` call does not move the caller's own stream.
#' @return Object of class `"tragweite"`: a list with `power`,
#'   `mcse`, `n_detected`, `n_valid`, `n_failed`, `design`, `effect`,
#'   `effect_size`, `criterion`, `alpha` (`NA` where it does not act),
#'   `fa_note` (`NULL` unless the detector's false-alarm rate is
#'   off-nominal), `reps`, `seed`.
#' @examples
#' \donttest{
#' tragweite("moderation", tau = 0.5, tau_x = 0.5,
#'           n_therapists = 20, caseload = 10, reps = 50, seed = 1)
#' }
#' @export
tragweite <- function(effect = c("ate", "moderation", "matching"),
                      tau = 0.5, tau_x = 0, tau_c = 0, tau_xc = 0,
                      n_therapists, caseload, n_sessions = 4,
                      icc = 0.10, z_level = c("patient", "therapist"),
                      tau_xc_form = "product", window_delta = 0.5,
                      rel_x = 1, rel_c = 1,
                      detector = c("lmer", "duet_s"),
                      duet_rule = c(delta_sd = 0.018, cut_sd = 0.52),
                      reps = 200, alpha = 0.05, seed) {
  effect <- match.arg(effect)
  z_level <- match.arg(z_level)
  detector <- match.arg(detector)
  stopifnot(rel_x > 0, rel_x <= 1, rel_c > 0, rel_c <= 1)
  .check_alpha(alpha)
  .check_reps(reps)
  if (detector == "duet_s" && effect != "matching") {
    stop("detector = \"duet_s\" is built for effect = \"matching\".",
         call. = FALSE)
  }
  if (detector == "duet_s") {
    # Fail fast and honestly: the sampler this detector needs is not
    # in the public windkanal yet (see .require_duet_s).
    .require_duet_s()
  }
  if (detector == "duet_s" && !isTRUE(all.equal(alpha, 0.05))) {
    warning("alpha does not act on the duet_s detector; `duet_rule` ",
            "makes the decision. For a rule calibrated on your own ",
            "world stream see ceiling_certificate().", call. = FALSE)
  }
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  if (missing(n_therapists) || missing(caseload)) {
    stop("`n_therapists` and `caseload` describe the design under ",
         "evaluation and are mandatory.", call. = FALSE)
  }
  rng_old <- .rng_state()
  on.exit(.rng_restore(rng_old), add = TRUE)
  detected <- rep(NA, reps)
  first_error <- NULL
  for (i in seq_len(reps)) {
    s <- windkanal::sim_stream(
      n_therapists = n_therapists, patients_per_therapist = caseload,
      n_sessions = n_sessions, icc = icc, z_level = z_level,
      tau = tau, tau_x = tau_x, tau_c = tau_c, tau_xc = tau_xc,
      tau_xc_form = tau_xc_form, fenster_delta = window_delta,
      seed = seed + i)
    if (effect == "matching" && is.null(s$therapist_c)) {
      # Null calibration: in the field C is MEASURED whether or not it
      # acts, so attach an inert attribute from a decoupled seed
      # stream. Without it the triple test ran against a constant and
      # the resulting subscript error was swallowed into a NaN power.
      set.seed(seed + i + 1e7)
      c_inert <- stats::rnorm(n_therapists)
      s$therapist_c <- c_inert[s$therapist_id]
    }
    # Measurement error on both sides: the analysis sees observed
    # versions (squared correlation with the truth = rel); the truth
    # in the generator is untouched. Own seed streams, so the world
    # draws stay stable.
    if (rel_x < 1) {
      pid <- unique(s$patient_id)
      set.seed(seed + i + 2e7)
      ex <- stats::rnorm(length(pid))
      s$x <- sqrt(rel_x) * s$x +
        sqrt(1 - rel_x) * ex[match(s$patient_id, pid)]
    }
    if (rel_c < 1 && !is.null(s$therapist_c)) {
      tid <- unique(s$therapist_id)
      set.seed(seed + i + 3e7)
      ec <- stats::rnorm(length(tid))
      s$therapist_c <- sqrt(rel_c) * s$therapist_c +
        sqrt(1 - rel_c) * ec[match(s$therapist_id, tid)]
    }
    detected[i] <- tryCatch(
      if (detector == "duet_s") .detect_duet_s(s, duet_rule)
      else .detect(s, effect, alpha),
                            error = function(e) {
                              if (is.null(first_error)) {
                                first_error <<- conditionMessage(e)
                              }
                              NA
                            })
  }
  n_valid <- sum(!is.na(detected))
  if (n_valid == 0) {
    stop("All ", reps, " worlds failed to fit -- first error: ",
         first_error, call. = FALSE)
  }
  if (!is.null(first_error)) {
    warning(reps - n_valid, " of ", reps, " fits failed -- ",
            "first error: ", first_error, call. = FALSE)
  }
  power <- mean(detected, na.rm = TRUE)
  out <- list(
    power = power,
    mcse = sqrt(power * (1 - power) / max(n_valid, 1)),
    n_detected = sum(detected, na.rm = TRUE),
    n_valid = n_valid,
    n_failed = reps - n_valid,
    design = list(n_therapists = n_therapists, caseload = caseload,
                  n_sessions = n_sessions, icc = icc,
                  z_level = z_level, N = n_therapists * caseload,
                  tau_xc_form = tau_xc_form, rel_x = rel_x,
                  rel_c = rel_c, detector = detector),
    effect = effect,
    effect_size = c(tau = tau, tau_x = tau_x, tau_c = tau_c,
                    tau_xc = tau_xc),
    criterion = .criterion_label(effect, detector, duet_rule),
    alpha = if (detector == "duet_s") NA_real_ else alpha,
    fa_note = if (detector == "duet_s") .duet_fa_note() else NULL,
    reps = reps, seed = seed)
  class(out) <- "tragweite"
  out
}

# Argument guards. alpha steers the answer directly: the typo
# alpha = 5 used to make every world a hit and returned power 1.000
# with a positive verdict.
#' @keywords internal
.check_alpha <- function(alpha) {
  stopifnot(
    "`alpha` must be a single number in (0, 1) -- alpha = 5 instead of .05 would count every world as a hit." =
      is.numeric(alpha) && length(alpha) == 1L && !is.na(alpha) &&
        alpha > 0 && alpha < 1)
  invisible(TRUE)
}

#' @keywords internal
.check_reps <- function(reps) {
  stopifnot(
    "`reps` must be a single number >= 1 -- without worlds there is no power." =
      is.numeric(reps) && length(reps) == 1L && !is.na(reps) && reps >= 1)
  invisible(TRUE)
}

# RNG hygiene: the world loop sets seeds in the global stream. Without
# a restore the caller would draw from a foreign stream after a
# tragweite() call.
#' @keywords internal
.rng_state <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

#' @keywords internal
.rng_restore <- function(state) {
  if (is.null(state)) {
    suppressWarnings(rm(list = ".Random.seed", envir = .GlobalEnv))
  } else {
    assign(".Random.seed", state, envir = .GlobalEnv)
  }
  invisible(NULL)
}

# ---- Availability guard for the v2S sampler ----------------------
# detector = "duet_s" and ceiling_certificate() call
# windkanal::fit_cate_dyade_v2s(). The public windkanal 0.3.0 does not
# ship that function (it exports fit_cate_dyade_v2p only), so the
# honest behavior is to say so at the point of call rather than to
# fail deep inside a world loop with a "could not find function"
# message.
# Export-scoped on purpose: the source implementation reached the
# sampler as windkanal::fit_cate_dyade_v2s, which resolves against
# windkanal's EXPORTS. Matching that here keeps tragweite from
# silently picking up an internal or experimental object of the same
# name in some future windkanal.
#' @keywords internal
.duet_s_available <- function() {
  requireNamespace("windkanal", quietly = TRUE) &&
    "fit_cate_dyade_v2s" %in% getNamespaceExports("windkanal")
}

# `what` names the entry point the user actually called:
# ceiling_certificate() has no `detector` argument, so quoting one at
# that call site would name a parameter the caller cannot pass.
#' @keywords internal
.require_duet_s <- function(what = "detector = 'duet_s'",
                            remedy = paste("Use detector = 'lmer' for the",
                                           "Satterthwaite criterion, which",
                                           "runs against windkanal 0.3.0",
                                           "in full.")) {
  if (!.duet_s_available()) {
    stop(what, " requires windkanal's v2S sampler ",
         "(fit_cate_dyade_v2s), which is planned for windkanal 0.4; ",
         "the public windkanal 0.3.0 does not ship it. ", remedy,
         call. = FALSE)
  }
  invisible(TRUE)
}

# The frozen amplitude rule has NO nominal false-alarm rate -- the
# number below is measured, not asserted.
#' @keywords internal
.duet_fa_note <- function() {
  paste("alpha does not govern this detector -- the frozen amplitude",
        "rule does, and its false-alarm rate is not nominal: measured",
        ".083 (calibration stream) and .067 (fresh stream) on a",
        "standard design (20 x 10, 4 sessions, ICC .10, rel .8/.8, 60",
        "null worlds each, seed 4242), not .05. For a planning",
        "statement recalibrate the cut on your own world stream with",
        "ceiling_certificate().")
}

# Real-data-computable detection tests (Satterthwaite mixed models).
#' @keywords internal
.detect <- function(s, effect, alpha) {
  d <- as.data.frame(s)
  if (effect == "matching" && is.null(d$therapist_c)) {
    stop("matching requires a therapist attribute in the stream; ",
         "a constant c = 0 would make the triple test degenerate.",
         call. = FALSE)
  }
  d$c <- if (!is.null(d$therapist_c)) d$therapist_c else 0
  f <- switch(effect,
    ate        = score ~ z + (1 | therapist_id) + (1 | patient_id),
    moderation = score ~ z * x + (1 | therapist_id) + (1 | patient_id),
    matching   = score ~ z * x * c + (1 | therapist_id) + (1 | patient_id))
  fit <- suppressMessages(suppressWarnings(
    lmerTest::lmer(f, data = d)))
  co <- stats::coef(summary(fit))
  term <- switch(effect, ate = "z", moderation = "z:x",
                 matching = "z:x:c")
  co[term, "Pr(>|t|)"] < alpha
}

# The v2S amplitude detector: frozen rule with a strict cut.
# WARNING: the static rule is STREAM- and STRUCTURE-sensitive (its
# false-alarm rate has been observed to move from .06 to .15 across
# world streams) -- for a defensible planning statement recalibrate
# the rule on the target stream; ceiling_certificate() always does,
# without a switch.
#
# mcmc_seed is deliberately CONSTANT and does not hang off the user
# seed: the sampler stream stays fixed across worlds and designs, so
# that a power difference is a design difference and not sampler
# noise (common random numbers on the estimator side). The price is
# that sampler sensitivity itself cannot be probed through the user
# seed -- set this parameter for that instead.
#' @keywords internal
.duet_p_sd <- function(s, delta_sd, mcmc_seed = 1, nburn = 500,
                       nsim = 500) {
  .require_duet_s()
  # Looked up dynamically rather than as windkanal::fit_cate_dyade_v2s:
  # the function does not exist in the public windkanal yet, and a
  # static reference to it would be an unresolvable dependency. The
  # guard above guarantees it is exported before we get here, and
  # getExportedValue keeps the same export scoping the static form had.
  sampler <- getExportedValue("windkanal", "fit_cate_dyade_v2s")
  fit <- sampler(s, nburn = nburn, nsim = nsim, seed = mcmc_seed)
  # "zerlegung" (decomposition) is the attribute name set by the
  # windkanal sampler; it is that package's contract, not ours.
  z <- attr(fit, "zerlegung")
  mean(z$sd_h > delta_sd)
}

#' @keywords internal
.detect_duet_s <- function(s, duet_rule) {
  .duet_p_sd(s, duet_rule[["delta_sd"]]) > duet_rule[["cut_sd"]]
}

#' @keywords internal
.criterion_label <- function(effect, detector = "lmer",
                             duet_rule = c(delta_sd = 0.018,
                                           cut_sd = 0.52)) {
  if (detector == "duet_s") {
    return(sprintf(paste("v2S amplitude rule P{sd(h) > %.3f} > %.3f",
                         "(frozen, strict cut; recalibrate on target",
                         "stream)"),
                   duet_rule[["delta_sd"]], duet_rule[["cut_sd"]]))
  }
  switch(effect,
    ate = "Satterthwaite test of the treatment coefficient (mixed model, therapist + patient intercepts)",
    moderation = "Satterthwaite test of the z:x interaction (mixed model)",
    matching = "Satterthwaite test of the z:x:c triple interaction (mixed model; parametric baseline criterion)")
}

#' @export
print.tragweite <- function(x, target = 0.80, ...) {
  cat(sprintf("tragweite -- detection probability of the design\n"))
  cat(sprintf("  Effect sought : %s (%s)\n", x$effect,
              paste(names(x$effect_size)[x$effect_size != 0],
                    x$effect_size[x$effect_size != 0],
                    sep = " = ", collapse = ", ")))
  cat(sprintf("  Design        : %d therapists x caseload %d (N = %d), %d sessions, ICC %.2f, %s-level assignment\n",
              x$design$n_therapists, x$design$caseload, x$design$N,
              x$design$n_sessions, x$design$icc, x$design$z_level))
  cat(sprintf("  Criterion     : %s\n", x$criterion))
  cat(sprintf("  Power         : %.3f [MCSE %.3f] (%d/%d worlds; %d failed fits)\n",
              x$power, x$mcse, x$n_detected, x$n_valid, x$n_failed))
  if (!is.null(x$fa_note)) {
    cat(strwrap(x$fa_note, width = 78, prefix = "                  ",
                initial = "  False alarms  : "), sep = "\n")
    cat("\n")
  }
  verdict <- if (x$power >= target) {
    sprintf("at or above the %.0f%% convention -- the design can carry this question.", 100 * target)
  } else {
    sprintf("below the %.0f%% convention -- a real effect of this size would usually be missed.", 100 * target)
  }
  cat(sprintf("  Verdict       : %s\n", verdict))
  invisible(x)
}
