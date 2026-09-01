#' A preregistrable detection certificate for a pairing design
#'
#' The prospective counterpart to [design_ceiling()], for studies that
#' randomize patient-therapist *pairings* and ask whether the pairing
#' effect could be detected at all. For one concrete instantiation
#' (design, reliability of the two baseline batteries, congruence
#' form) `ceiling_certificate()` computes the detection certificate
#' under the dyadic amplitude rule, in two steps:
#'
#' 1. **Recalibration on your own world stream** (mandatory, not a
#'    switch): `reps` null worlds of the same design, strict cut
#'    convention `cut* = min{k : FA(k) <= alpha}`.
#' 2. **Power across the effect grid** with the recalibrated cut.
#'
#' Step 1 is not optional because the frozen amplitude rule is
#' stream-sensitive: its false-alarm rate has been observed to move
#' from .06 to .15 across world streams. A certificate that carried a
#' cut over from somewhere else would be self-deception.
#'
#' @section Requires the v2S sampler (windkanal 0.4):
#' This function is built on `windkanal::fit_cate_dyade_v2s()`, which
#' the public windkanal 0.3.0 does **not** ship; it is planned for
#' windkanal 0.4. Until then `ceiling_certificate()` stops with an
#' explanatory error. Everything else in tragweite -- including
#' [design_ceiling()] with its default `detector = "lmer"` -- runs
#' against the public windkanal 0.3.0 in full.
#'
#' @param n_therapists,caseload,n_sessions,icc Design of the planned
#'   study.
#' @param rel_x,rel_c Reliability of the two pre-baseline batteries
#'   (default .8/.8).
#' @param tau_xc_form,window_delta Congruence form of the target
#'   hypothesis, passed to the generator.
#' @param fixed Non-varied effect parameters, as in
#'   [design_ceiling()]. Unknown names abort; missing ones are filled
#'   up from the default.
#' @param effect_grid Ascending grid of candidate effects.
#' @param target Target power defining the ceiling (default .80).
#' @param delta_sd Amplitude threshold of the rule (default .018); the
#'   CUT is always recalibrated.
#' @param reps Worlds per grid point AND for the recalibration
#'   (>= 1; below `1/alpha` the cut rests on a single world and the
#'   function warns).
#' @param alpha False-alarm target of the recalibration (default .05;
#'   a single number in (0, 1)).
#' @param seed Base seed (mandatory). Null calibration uses `seed`,
#'   grid point k uses `seed + k * 10000`. The global RNG state is
#'   restored on exit.
#' @param quiet Suppress progress output.
#' @return Object of class `tragweite_certificate`: `design`, `fixed`, `rule`
#'   (delta, recalibrated cut, in-sample FA, hold-out FA), `grid`
#'   (effect/power/mcse/failed fits), `failures`, `ceiling`, `target`.
#'
#'   On the two false-alarm numbers: the **in-sample FA** is by
#'   construction the largest grid value `<= alpha`, namely
#'   `floor(alpha * reps)/reps`, and is therefore no evidence that the
#'   recalibration hit its target -- it cannot be anything else. The
#'   **hold-out FA** calibrates the cut on the first half of the null
#'   worlds and measures it on the second: a genuine out-of-sample
#'   number, resting on half the calibration base and therefore
#'   pessimistic rather than optimistic. It costs no extra simulation.
#'   The split runs by position because the null worlds are
#'   independent draws -- their order carries no structure.
#' @seealso [design_ceiling()] for the retrospective ceiling under the
#'   Satterthwaite criterion, which needs no v2S sampler.
#' @export
ceiling_certificate <- function(n_therapists, caseload, n_sessions = 4,
                                icc = 0.10, rel_x = 0.8, rel_c = 0.8,
                                tau_xc_form = "product",
                                window_delta = 0.5,
                                fixed = list(tau = 0.3, tau_x = 0.2,
                                             tau_c = 0.2),
                                effect_grid = seq(0.10, 0.40, by = 0.05),
                                target = 0.80, delta_sd = 0.018,
                                reps = 200, alpha = 0.05, seed,
                                quiet = FALSE) {
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  stopifnot(is.numeric(effect_grid), !is.unsorted(effect_grid),
            all(effect_grid > 0))
  .check_alpha(alpha)
  .check_reps(reps)
  fixed <- .cert_fixed(fixed)
  .require_duet_s(
    what = "ceiling_certificate()",
    remedy = paste("For a detection ceiling that needs no v2S sampler,",
                   "use design_ceiling(), which runs against windkanal",
                   "0.3.0 in full."))
  rng_old <- .rng_state()
  on.exit(.rng_restore(rng_old), add = TRUE)
  first_error <- NULL

  one_world <- function(tau_xc, base_seed, i) {
    s <- windkanal::sim_stream(
      n_therapists = n_therapists, patients_per_therapist = caseload,
      n_sessions = n_sessions, icc = icc, z_level = "patient",
      tau = fixed$tau, tau_x = fixed$tau_x, tau_c = fixed$tau_c,
      tau_xc = tau_xc, tau_xc_form = tau_xc_form,
      fenster_delta = window_delta, seed = base_seed + i)
    if (is.null(s$therapist_c)) {          # null world: C is measured
      set.seed(base_seed + i + 1e7)        # in the field -> attach inert
      c_inert <- stats::rnorm(n_therapists)
      s$therapist_c <- c_inert[s$therapist_id]
    }
    if (rel_x < 1) {
      pid <- unique(s$patient_id)
      set.seed(base_seed + i + 2e7)
      ex <- stats::rnorm(length(pid))
      s$x <- sqrt(rel_x) * s$x +
        sqrt(1 - rel_x) * ex[match(s$patient_id, pid)]
    }
    if (rel_c < 1) {
      tid <- unique(s$therapist_id)
      set.seed(base_seed + i + 3e7)
      ec <- stats::rnorm(length(tid))
      s$therapist_c <- sqrt(rel_c) * s$therapist_c +
        sqrt(1 - rel_c) * ec[match(s$therapist_id, tid)]
    }
    .duet_p_sd(s, delta_sd)
  }

  # One failed fit must not throw away a long certificate run: count
  # and carry on, as tragweite() does.
  world_p_sd <- function(tau_xc, base_seed) {
    vapply(seq_len(reps), function(i) {
      tryCatch(one_world(tau_xc, base_seed, i),
               error = function(e) {
                 if (is.null(first_error)) {
                   first_error <<- conditionMessage(e)
                 }
                 NA_real_
               })
    }, numeric(1))
  }

  # ---- 1. Recalibration on our own stream ----
  if (!quiet) cat("Recalibration (", reps, " null worlds) ...\n", sep = "")
  p_null <- world_p_sd(0, seed)
  n_null <- sum(!is.na(p_null))
  if (n_null < 2) {
    stop("Recalibration with ", n_null, " of ", reps,
         " valid null worlds -- at least 2 required.",
         if (!is.null(first_error)) {
           paste0(" First error: ", first_error)
         },
         call. = FALSE)
  }
  cut_star <- .cert_cut(p_null, alpha)
  fa_calib <- mean(p_null > cut_star, na.rm = TRUE)
  ho <- .cert_holdout(p_null, alpha)
  if (!quiet) {
    cat(sprintf("  cut* = %.3f (FA in-sample %.3f, hold-out %.3f)\n",
                cut_star, fa_calib, ho$fa))
  }

  # ---- 2. Power across the grid ----
  g <- data.frame(effect = effect_grid, power = NA_real_,
                  mcse = NA_real_, n_valid = NA_integer_,
                  n_failed = NA_integer_)
  for (k in seq_along(effect_grid)) {
    p_sig <- world_p_sd(effect_grid[k], seed + k * 10000)
    valid <- sum(!is.na(p_sig))
    g$power[k] <- mean(p_sig > cut_star, na.rm = TRUE)
    g$mcse[k] <- sqrt(g$power[k] * (1 - g$power[k]) / max(valid, 1))
    g$n_valid[k] <- valid
    g$n_failed[k] <- reps - valid
    if (!quiet) {
      cat(sprintf("  tau_xc = %.2f: power %.3f [MCSE %.3f] (%d/%d worlds)\n",
                  effect_grid[k], g$power[k], g$mcse[k], valid, reps))
    }
  }
  failed_total <- (reps - n_null) + sum(g$n_failed)
  if (failed_total > 0) {
    warning(failed_total, " world fits failed -- first error: ",
            first_error, call. = FALSE)
  }
  idx <- which(g$power >= target)
  out <- list(
    design = list(n_therapists = n_therapists, caseload = caseload,
                  N = n_therapists * caseload, n_sessions = n_sessions,
                  icc = icc, rel_x = rel_x, rel_c = rel_c,
                  tau_xc_form = tau_xc_form, z_level = "patient"),
    fixed = fixed,
    rule = list(delta_sd = delta_sd, cut = cut_star,
                fa_in_sample = fa_calib, fa_holdout = ho$fa,
                fa_holdout_mcse = ho$mcse, holdout_cut = ho$cut,
                holdout_n_calib = ho$n_calib,
                holdout_n_test = ho$n_test, alpha = alpha,
                reps_calibration = reps, n_valid_calibration = n_null),
    failures = list(calibration = reps - n_null,
                    grid = g$n_failed,
                    first_error = first_error),
    grid = g, target = target,
    ceiling = if (length(idx)) g$effect[min(idx)] else NA_real_)
  class(out) <- "tragweite_certificate"
  out
}

# ---- The cut arithmetic, pulled out and therefore testable --------
# Strict cut convention: cut* = min{k : FA(k) <= alpha} on the grid of
# observed null values. That is the floor(alpha*n)+1-largest null
# observation -- an extreme order statistic whose base gets thin with
# small `reps`.
#' @keywords internal
.cert_cut <- function(p_null, alpha, warn = TRUE) {
  p <- p_null[!is.na(p_null)]
  if (!length(p)) stop("no valid null worlds", call. = FALSE)
  if (warn && length(p) < 1 / alpha) {
    warning("only ", length(p), " valid null worlds at alpha = ",
            alpha, ": the false-alarm target is then exactly 0 and a ",
            "SINGLE world fixes the cut. For a defensible certificate ",
            "use reps >= ", ceiling(5 / alpha), ".", call. = FALSE)
  }
  if (warn && length(unique(p)) < 3) {
    warning("degenerate null distribution (", length(unique(p)),
            " distinct values) -- the cut is not estimable this way.",
            call. = FALSE)
  }
  cuts <- sort(unique(c(0, p, 1)))
  fa <- vapply(cuts, function(k) mean(p > k), numeric(1))
  min(cuts[fa <= alpha])
}

# Out-of-sample FA without extra simulation: cut on the first half of
# the null worlds, measured on the second. The halved calibration base
# makes the number pessimistic rather than optimistic -- it replaces
# the tautological in-sample FA as evidence.
#' @keywords internal
.cert_holdout <- function(p_null, alpha) {
  p <- p_null[!is.na(p_null)]
  empty <- list(fa = NA_real_, mcse = NA_real_, cut = NA_real_,
                n_calib = NA_integer_, n_test = NA_integer_)
  if (length(p) < 4) return(empty)
  a <- seq_len(floor(length(p) / 2))
  b <- setdiff(seq_along(p), a)
  cut_a <- .cert_cut(p[a], alpha, warn = FALSE)
  fa <- mean(p[b] > cut_a)
  list(fa = fa, mcse = sqrt(fa * (1 - fa) / length(b)), cut = cut_a,
       n_calib = length(a), n_test = length(b))
}

# `fixed` fills in from the default and lets no unknown name fall
# through silently: without this a missing tau_x would travel into the
# generator as NULL and the run would abort with a message from deep
# inside the generator.
#' @keywords internal
.cert_fixed <- function(fixed) {
  default <- eval(formals(ceiling_certificate)$fixed)
  unknown <- setdiff(names(fixed), names(default))
  if (length(unknown)) {
    stop("`fixed` only knows ", paste(names(default), collapse = ", "),
         " -- unknown: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  utils::modifyList(default, fixed)
}

#' @export
print.tragweite_certificate <- function(x, ...) {
  d <- x$design
  cat("Detection certificate for a randomized pairing design\n")
  cat(sprintf("  Design    : %d therapists x caseload %d (N = %d), %d sessions, ICC %.2f\n",
              d$n_therapists, d$caseload, d$N, d$n_sessions, d$icc))
  cat(sprintf("  Batteries : rel_x = %.2f, rel_c = %.2f | form: %s\n",
              d$rel_x, d$rel_c, d$tau_xc_form))
  cat(sprintf("  World     : %s (fixed), assignment at %s level -- the ceiling holds for THIS configuration\n",
              paste(names(x$fixed), unlist(x$fixed), sep = " = ",
                    collapse = ", "), d$z_level))
  cat(sprintf("  Rule      : P{sd(h) > %.3f} > %.3f (cut recalibrated on THIS stream, %d null worlds)\n",
              x$rule$delta_sd, x$rule$cut,
              x$rule$n_valid_calibration))
  cat(sprintf("  Alarms    : in-sample %.3f -- by construction floor(alpha*n)/n, no evidence of calibration\n",
              x$rule$fa_in_sample))
  if (is.na(x$rule$fa_holdout)) {
    cat("              hold-out  : not computable (too few null worlds)\n")
  } else {
    cat(sprintf("              hold-out %.3f [MCSE %.3f] -- cut on %d worlds, measured on %d others (half the base, pessimistic)\n",
                x$rule$fa_holdout, x$rule$fa_holdout_mcse,
                x$rule$holdout_n_calib, x$rule$holdout_n_test))
  }
  if (is.na(x$ceiling)) {
    cat(sprintf("  Ceiling   : NOT REACHED -- no effect on the grid reaches %.0f%% power.\n",
                100 * x$target))
    cat("  Verdict   : this design cannot see the target hypothesis in the range examined.\n")
  } else {
    cat(sprintf("  Ceiling   : %.2f (smallest effect at >= %.0f%% power)\n",
                x$ceiling, 100 * x$target))
  }
  cat("  Power grid:\n")
  for (k in seq_len(nrow(x$grid))) {
    cat(sprintf("    tau_xc = %.2f : %.3f [MCSE %.3f] (%d worlds)\n",
                x$grid$effect[k], x$grid$power[k], x$grid$mcse[k],
                x$grid$n_valid[k]))
  }
  failed_total <- x$failures$calibration + sum(x$failures$grid)
  if (failed_total > 0) {
    cat(sprintf("  Failures  : %d world fits failed (%d in the calibration) -- first error: %s\n",
                failed_total, x$failures$calibration,
                x$failures$first_error))
  }
  cat("  Note      : freeze this certificate before data collection (it is preregistrable);\n")
  cat("              the hold-out number is the false-alarm evidence, not the in-sample one.\n")
  invisible(x)
}
