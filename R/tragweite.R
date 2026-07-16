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
#' Detection criteria (v0, all Satterthwaite mixed-model tests --
#' validated in the windkanal hardening study; forest-based criteria
#' for matching will follow the dyadic method's selection rule):
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
#' @param reps Number of simulated worlds (default 200).
#' @param alpha Test level of the detection criterion (default .05).
#' @param seed Base random seed (mandatory; world `i` uses
#'   `seed + i`).
#' @return Object of class `"tragweite"`: a list with `power`,
#'   `mcse`, `n_detected`, `n_valid`, `n_failed`, `design`, `effect`,
#'   `criterion`, `alpha`, `seed`.
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
                      reps = 200, alpha = 0.05, seed) {
  effect <- match.arg(effect)
  z_level <- match.arg(z_level)
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  if (missing(n_therapists) || missing(caseload)) {
    stop("`n_therapists` and `caseload` describe the design under ",
         "evaluation and are mandatory.", call. = FALSE)
  }
  detected <- rep(NA, reps)
  first_error <- NULL
  for (i in seq_len(reps)) {
    s <- windkanal::sim_stream(
      n_therapists = n_therapists, patients_per_therapist = caseload,
      n_sessions = n_sessions, icc = icc, z_level = z_level,
      tau = tau, tau_x = tau_x, tau_c = tau_c, tau_xc = tau_xc,
      seed = seed + i)
    if (effect == "matching" && is.null(s$therapist_c)) {
      # Null calibration: in the field C is MEASURED whether or not
      # it acts, so a null world carries a measured, inert attribute.
      # Drawn from its own seed stream so the world's draws are
      # untouched.
      set.seed(seed + i + 1e7)
      c_inert <- stats::rnorm(n_therapists)
      s$therapist_c <- c_inert[s$therapist_id]
    }
    detected[i] <- tryCatch(.detect(s, effect, alpha),
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
                  z_level = z_level, N = n_therapists * caseload),
    effect = effect,
    effect_size = c(tau = tau, tau_x = tau_x, tau_c = tau_c,
                    tau_xc = tau_xc),
    criterion = .criterion_label(effect),
    alpha = alpha, reps = reps, seed = seed)
  class(out) <- "tragweite"
  out
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
  ziel <- switch(effect, ate = "z", moderation = "z:x",
                 matching = "z:x:c")
  co[ziel, "Pr(>|t|)"] < alpha
}

#' @keywords internal
.criterion_label <- function(effect) {
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
  verdict <- if (x$power >= target) {
    sprintf("at or above the %.0f%% convention -- the design can carry this question.", 100 * target)
  } else {
    sprintf("below the %.0f%% convention -- a real effect of this size would usually be missed.", 100 * target)
  }
  cat(sprintf("  Verdict       : %s\n", verdict))
  invisible(x)
}
