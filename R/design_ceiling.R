#' Retrospective design diagnosis: the detection ceiling of a design
#'
#' The third pillar: given the structure of an EXISTING study
#' (therapists, caseload, ICC, sessions), which effect sizes could it
#' ever have detected? `design_ceiling()` walks an ascending effect
#' grid with common random numbers (same base seed per cell) and
#' returns the smallest effect the design detects with the target
#' power -- its detection ceiling. Effects below that line fall short
#' of the target convention; they are not unobservable, and the grid
#' states with what probability each of them was in fact detected. The
#' further an effect sits below the ceiling, the more a published null
#' from such a design is a fact about the design rather than about the
#' world ("tool null vs. world null") -- how far, the measured power
#' below the ceiling says.
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
#'   A partial `fixed` is filled in from this default.
#' @param target Target power defining the ceiling (default 0.80).
#' @param tau_xc_form,window_delta,rel_x,rel_c,detector,duet_rule
#'   As in [tragweite()] (congruence form, reliability on both sides,
#'   choice of detector). Note that `detector = "duet_s"` is not
#'   available against the public windkanal 0.3.0; see [tragweite()].
#' @param reps,alpha,seed As in [tragweite()]; `seed` mandatory.
#' @param quiet Suppress per-cell progress (default FALSE).
#' @return Object of class `tragweite_ceiling`: list with `ceiling`
#'   (smallest detectable effect, `NA` if the grid never reaches
#'   target), `grid` (data.frame effect/power/mcse), `fixed` (the
#'   parameters actually held fixed -- the varied one lives in `grid`),
#'   `design`, `target`, `effect`, `criterion`.
#' @export
design_ceiling <- function(effect,
                           effect_grid = seq(0.1, 0.6, by = 0.05),
                           n_therapists, caseload, icc = 0.10,
                           n_sessions = 4, z_level = "patient",
                           tau_xc_form = "product", window_delta = 0.5,
                           rel_x = 1, rel_c = 1,
                           detector = "lmer",
                           duet_rule = c(delta_sd = 0.018, cut_sd = 0.52),
                           fixed = list(tau = 0.3, tau_x = 0,
                                        tau_c = 0, tau_xc = 0),
                           target = 0.80, reps = 200, alpha = 0.05,
                           seed, quiet = FALSE) {
  if (missing(seed)) stop("`seed` is mandatory.", call. = FALSE)
  stopifnot(is.numeric(effect_grid), !is.unsorted(effect_grid),
            all(effect_grid > 0), target > 0, target < 1)
  .check_alpha(alpha)
  .check_reps(reps)
  # A partial `fixed` fills in from the DOCUMENTED default here, not
  # from tragweite()'s own defaults.
  fixed <- utils::modifyList(eval(formals(design_ceiling)$fixed), fixed)
  varied <- switch(effect, ate = "tau", moderation = "tau_x",
                   matching = "tau_xc",
                   stop("unknown effect", call. = FALSE))
  g <- data.frame(effect = effect_grid, power = NA_real_,
                  mcse = NA_real_)
  crit <- NULL
  for (k in seq_along(effect_grid)) {
    args <- fixed
    args[[varied]] <- effect_grid[k]
    r <- do.call(tragweite, c(list(effect = effect),
      args, list(n_therapists = n_therapists, caseload = caseload,
                 n_sessions = n_sessions, icc = icc,
                 z_level = z_level, tau_xc_form = tau_xc_form,
                 window_delta = window_delta, rel_x = rel_x,
                 rel_c = rel_c, detector = detector,
                 duet_rule = duet_rule, reps = reps, alpha = alpha,
                 seed = seed)))   # same seed = CRN across the grid
    g$power[k] <- r$power; g$mcse[k] <- r$mcse
    crit <- r$criterion
    if (!quiet) {
      cat(sprintf("  %s = %.2f: power %.3f [MCSE %.3f]\n",
                  varied, effect_grid[k], r$power, r$mcse))
    }
  }
  idx <- which(g$power >= target)
  out <- list(
    ceiling = if (length(idx)) g$effect[min(idx)] else NA_real_,
    grid = g, target = target, effect = effect, criterion = crit,
    fixed = fixed[setdiff(names(fixed), varied)],  # without the varied one
    design = list(n_therapists = n_therapists, caseload = caseload,
                  N = n_therapists * caseload, icc = icc,
                  n_sessions = n_sessions, z_level = z_level,
                  rel_x = rel_x, rel_c = rel_c,
                  tau_xc_form = tau_xc_form, detector = detector))
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
  if (!is.null(x$fixed)) {
    cat(sprintf("  Fixed     : %s\n",
                paste(names(x$fixed), unlist(x$fixed), sep = " = ",
                      collapse = ", ")))
  }
  if (!is.null(d$rel_x)) {
    cat(sprintf("              rel_x = %.2f, rel_c = %.2f, form %s, detector %s -- the ceiling holds for THIS world\n",
                d$rel_x, d$rel_c, d$tau_xc_form, d$detector))
  }
  if (is.na(x$ceiling)) {
    top <- which.max(x$grid$effect)
    cat(sprintf("  Ceiling   : NOT REACHED -- no effect on the grid reaches %.0f%% power.\n",
                100 * x$target))
    cat(sprintf("  Verdict   : even the largest effect on the grid (%.2f) reaches only %.2f power -- a null from this design constrains the effect weakly at best.\n",
                x$grid$effect[top], x$grid$power[top]))
  } else {
    below <- x$grid[x$grid$effect < x$ceiling, , drop = FALSE]
    cat(sprintf("  Ceiling   : %.2f (smallest effect at >= %.0f%% power)\n",
                x$ceiling, 100 * x$target))
    if (nrow(below)) {
      last <- below[nrow(below), ]
      cat(sprintf("  Verdict   : effects below %.2f fall short of the %.0f%% convention -- not unobservable: the grid measures %.2f detection probability at effect %.2f.\n",
                  x$ceiling, 100 * x$target, last$power,
                  last$effect))
    } else {
      cat(sprintf("  Verdict   : the grid reaches the %.0f%% convention at its lowest point (%.2f); smaller effects were not examined.\n",
                  100 * x$target, x$ceiling))
    }
  }
  invisible(x)
}
