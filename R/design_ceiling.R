#' Retrospective design diagnosis: the detection ceiling of a design
#'
#' The third pillar: given the structure of an EXISTING study
#' (therapists, caseload, ICC, sessions), which effect sizes could it
#' ever have detected? `design_ceiling()` walks an ascending effect
#' grid with common random numbers (same base seed per cell) and
#' returns the smallest effect the design detects with the target
#' power -- its detection ceiling. Effects below that line were never
#' observable: a published null from such a design is a fact about
#' the design, not about the world ("tool null vs. world null").
#'
#' Named `design_ceiling` (not `ceiling`) to avoid masking
#' `base::ceiling`.
#'
#' @param effect `"ate"`, `"moderation"`, or `"matching"` as in
#'   [tragweite()]; the grid varies the corresponding parameter
#'   (`tau`, `tau_x`, `tau_xc`) while the others stay at `fixed`.
#' @param effect_grid Ascending vector of candidate effect sizes
#'   (default `seq(0.1, 0.6, by = 0.05)`).
#' @param n_therapists,caseload,icc,n_sessions,z_level The design
#'   whose ceiling is diagnosed (scalars).
#' @param fixed Named list of the non-varied effect parameters
#'   (default `list(tau = 0.3, tau_x = 0, tau_c = 0, tau_xc = 0)`).
#' @param target Target power defining the ceiling (default 0.80).
#' @param reps,alpha,seed As in [tragweite()]; `seed` mandatory.
#' @param quiet Suppress per-cell progress (default FALSE).
#' @return Object of class `tragweite_ceiling`: list with `ceiling`
#'   (smallest detectable effect, `NA` if the grid never reaches
#'   target), `grid` (data.frame effect/power/mcse), `design`,
#'   `target`, `criterion`.
#' @export
design_ceiling <- function(effect,
                           effect_grid = seq(0.1, 0.6, by = 0.05),
                           n_therapists, caseload, icc = 0.10,
                           n_sessions = 4, z_level = "patient",
                           fixed = list(tau = 0.3, tau_x = 0,
                                        tau_c = 0, tau_xc = 0),
                           target = 0.80, reps = 200, alpha = 0.05,
                           seed, quiet = FALSE) {
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  stopifnot(is.numeric(effect_grid), !is.unsorted(effect_grid),
            all(effect_grid > 0), target > 0, target < 1)
  ziel <- switch(effect, ate = "tau", moderation = "tau_x",
                 matching = "tau_xc",
                 stop("unknown effect", call. = FALSE))
  g <- data.frame(effect = effect_grid, power = NA_real_,
                  mcse = NA_real_)
  crit <- NULL
  for (k in seq_along(effect_grid)) {
    args <- fixed
    args[[ziel]] <- effect_grid[k]
    r <- do.call(tragweite, c(list(effect = effect),
      args, list(n_therapists = n_therapists, caseload = caseload,
                 n_sessions = n_sessions, icc = icc,
                 z_level = z_level, reps = reps, alpha = alpha,
                 seed = seed)))   # same seed = common random numbers across the grid
    g$power[k] <- r$power; g$mcse[k] <- r$mcse
    crit <- r$criterion
    if (!quiet) {
      cat(sprintf("  %s = %.2f: power %.3f [MCSE %.3f]\n",
                  ziel, effect_grid[k], r$power, r$mcse))
    }
  }
  idx <- which(g$power >= target)
  out <- list(
    ceiling = if (length(idx)) g$effect[min(idx)] else NA_real_,
    grid = g, target = target, effect = effect, criterion = crit,
    design = list(n_therapists = n_therapists, caseload = caseload,
                  N = n_therapists * caseload, icc = icc,
                  n_sessions = n_sessions, z_level = z_level))
  class(out) <- "tragweite_ceiling"
  out
}

#' @export
print.tragweite_ceiling <- function(x, ...) {
  d <- x$design
  cat("design_ceiling -- what this design could ever have detected\n")
  cat(sprintf("  Design    : %d therapists x caseload %d (N = %d), %d sessions, ICC %.2f\n",
              d$n_therapists, d$caseload, d$N, d$n_sessions, d$icc))
  cat(sprintf("  Effect    : %s, grid %.2f-%.2f\n", x$effect,
              min(x$grid$effect), max(x$grid$effect)))
  if (is.na(x$ceiling)) {
    cat(sprintf("  Ceiling   : NOT REACHED -- no effect on the grid reaches %.0f%% power.\n",
                100 * x$target))
    cat("  Verdict   : a null finding from this design is uninformative across the whole grid.\n")
  } else {
    cat(sprintf("  Ceiling   : %.2f (smallest effect at >= %.0f%% power)\n",
                x$ceiling, 100 * x$target))
    cat(sprintf("  Verdict   : effects below %.2f were never observable with this design.\n",
                x$ceiling))
  }
  invisible(x)
}
