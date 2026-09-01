#' Detection probability across a design grid
#'
#' Evaluates [tragweite()] over the cross product of design vectors.
#' All grid points share the same world seeds (common random numbers),
#' so *comparisons between designs* carry less Monte Carlo noise than
#' the individual power estimates suggest.
#'
#' @param effect,tau,tau_x,tau_c,tau_xc As in [tragweite()]; the
#'   effect-size arguments may be vectors and become grid dimensions.
#' @param n_therapists,caseload,icc Vectors of design values to cross.
#' @param n_sessions,z_level,reps,alpha,seed As in [tragweite()]
#'   (scalar).
#' @param quiet Suppress per-point progress messages (default FALSE).
#' @return A `data.frame` with one row per grid point: the design
#'   columns, `power`, `mcse`, `n_failed`. Class
#'   `c("tragweite_grid", "data.frame")`.
#' @section Limitations:
#' `design_grid()` and [min_design()] do **not** pass `tau_xc_form`,
#'   `window_delta`, `rel_x`, `rel_c`, `detector` and `duet_rule`
#'   through to [tragweite()] -- unlike [design_ceiling()]. Every grid
#'   point therefore runs the product form, perfectly measured
#'   attributes and the `"lmer"` detector, whatever the surrounding
#'   analysis assumes. To grid over reliability or congruence form,
#'   call [tragweite()] in your own loop.
#' @export
design_grid <- function(effect, tau = 0.5, tau_x = 0, tau_c = 0,
                        tau_xc = 0, n_therapists, caseload,
                        icc = 0.10, n_sessions = 4,
                        z_level = "patient", reps = 200,
                        alpha = 0.05, seed, quiet = FALSE) {
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  .check_alpha(alpha)
  .check_reps(reps)
  g <- expand.grid(tau = tau, tau_x = tau_x, tau_c = tau_c,
                   tau_xc = tau_xc, n_therapists = n_therapists,
                   caseload = caseload, icc = icc)
  g$power <- g$mcse <- g$n_failed <- NA_real_
  for (k in seq_len(nrow(g))) {
    r <- tragweite(effect = effect, tau = g$tau[k], tau_x = g$tau_x[k],
                   tau_c = g$tau_c[k], tau_xc = g$tau_xc[k],
                   n_therapists = g$n_therapists[k],
                   caseload = g$caseload[k], n_sessions = n_sessions,
                   icc = g$icc[k], z_level = z_level, reps = reps,
                   alpha = alpha, seed = seed)  # same seeds: CRN
    g$power[k] <- r$power; g$mcse[k] <- r$mcse
    g$n_failed[k] <- r$n_failed
    if (!quiet) {
      message(sprintf("[%d/%d] %d therapists x %d | power %.3f",
                      k, nrow(g), g$n_therapists[k], g$caseload[k],
                      g$power[k]))
    }
  }
  class(g) <- c("tragweite_grid", "data.frame")
  attr(g, "effect") <- effect
  g
}

#' Smallest design reaching a target detection probability
#'
#' Evaluates candidate designs in order of increasing total N and
#' returns the first whose estimated power reaches `target`. The MCSE
#' is reported alongside; near the boundary, increase `reps` before
#' committing to a design.
#'
#' @inheritParams design_grid
#' @param target Target detection probability (default 0.80).
#' @return A list with `design` (the chosen row, or `NULL` if no
#'   candidate reaches the target) and `grid` (all evaluated rows,
#'   ordered by N).
#' @inheritSection design_grid Limitations
#' @export
min_design <- function(effect, tau = 0.5, tau_x = 0, tau_c = 0,
                       tau_xc = 0, n_therapists, caseload,
                       icc = 0.10, n_sessions = 4,
                       z_level = "patient", target = 0.80,
                       reps = 200, alpha = 0.05, seed,
                       quiet = FALSE) {
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  .check_alpha(alpha)
  .check_reps(reps)
  g <- expand.grid(n_therapists = n_therapists, caseload = caseload)
  g$N <- g$n_therapists * g$caseload
  g <- g[order(g$N, g$n_therapists), ]
  g$power <- g$mcse <- NA_real_
  hit <- NULL
  for (k in seq_len(nrow(g))) {
    r <- tragweite(effect = effect, tau = tau, tau_x = tau_x,
                   tau_c = tau_c, tau_xc = tau_xc,
                   n_therapists = g$n_therapists[k],
                   caseload = g$caseload[k], n_sessions = n_sessions,
                   icc = icc, z_level = z_level, reps = reps,
                   alpha = alpha, seed = seed)
    g$power[k] <- r$power; g$mcse[k] <- r$mcse
    if (!quiet) {
      message(sprintf("N = %4d (%d x %d): power %.3f", g$N[k],
                      g$n_therapists[k], g$caseload[k], g$power[k]))
    }
    if (r$power >= target) { hit <- g[k, ]; break }
  }
  list(design = hit, grid = g[!is.na(g$power), ])
}
