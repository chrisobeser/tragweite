#' Calibrated effect-size scenarios (anchored, with estimand caveats)
#'
#' Scenario entries anchor tragweite runs in the published literature
#' instead of round guesses. Every scenario carries its source and --
#' deliberately prominent -- an *estimand note*: published effect
#' sizes measure different quantities (strategy-level ATEs, allocation
#' contrasts, slope differences) and are **not** automatically
#' data-generating parameters. Where no defensible translation exists
#' yet, the scenario provides the anchor and says so, rather than
#' inventing a mapping. (The formal translation layer is the package's
#' first roadmap pillar.)
#'
#' Available scenarios:
#' * `"field_realistic"`: tau = 0.2, tau_x = 0.2 -- the windkanal
#'   realistic-effect configuration; a direct DGM specification, no translation
#'   involved.
#' * `"personalization_floor"`: anchor d = 0.14 (95% CI 0.08-0.20),
#'   the risk-of-bias-adjusted meta-analytic floor for personalized
#'   vs. standardized care (Nye, Delgadillo & Barkham, 2023, J Consult
#'   Clin Psychol, doi:10.1037/ccp0000820). Strategy-level ATE; no
#'   automatic DGM mapping -- use as a sensitivity anchor.
#' * `"matching_rct"`: anchor d = 0.50-0.75, prospective
#'   therapist-matching RCT (Constantino et al., 2021, JAMA
#'   Psychiatry, doi:10.1001/jamapsychiatry.2021.1221). Weekly-slope
#'   contrasts under measurement-feedback matching; no automatic DGM
#'   mapping.
#'
#' @param name Scenario name (see above).
#' @return A list with `name`, `anchor`, `source`, `estimand_note`,
#'   and `suggested_args` (a list ready for [tragweite()], or `NULL`
#'   where an honest translation does not exist yet).
#' @export
scenario <- function(name = c("field_realistic",
                              "personalization_floor",
                              "matching_rct")) {
  name <- match.arg(name)
  s <- switch(name,
    field_realistic = list(
      name = name,
      anchor = "tau = 0.2, tau_x = 0.2 (outcome-SD units)",
      source = "windkanal realistic-effect configuration, anchored to routine-care meta-analytic ranges",
      estimand_note = "Direct DGM specification; no translation involved.",
      suggested_args = list(tau = 0.2, tau_x = 0.2)),
    personalization_floor = list(
      name = name,
      anchor = "d = 0.14 [0.08, 0.20], risk-of-bias-adjusted",
      source = "Nye, Delgadillo & Barkham (2023), J Consult Clin Psychol, doi:10.1037/ccp0000820",
      estimand_note = paste("Strategy-level ATE of personalized vs.",
        "standardized care -- NOT a moderation amplitude. No",
        "defensible automatic mapping to tau_x exists yet; use as a",
        "sensitivity anchor and report the assumption you choose."),
      suggested_args = NULL),
    matching_rct = list(
      name = name,
      anchor = "d = 0.50-0.75 across outcomes",
      source = "Constantino et al. (2021), JAMA Psychiatry, doi:10.1001/jamapsychiatry.2021.1221",
      estimand_note = paste("Weekly-slope contrasts from a",
        "measurement-feedback matching RCT -- NOT a tau_xc amplitude.",
        "No defensible automatic mapping exists yet; use as an upper",
        "sensitivity anchor."),
      suggested_args = NULL))
  class(s) <- "tragweite_scenario"
  s
}

#' @export
print.tragweite_scenario <- function(x, ...) {
  cat(sprintf("Scenario '%s'\n  Anchor : %s\n  Source : %s\n  Note   : %s\n",
              x$name, x$anchor, x$source, x$estimand_note))
  if (!is.null(x$suggested_args)) {
    cat(sprintf("  Ready-to-use args: %s\n",
                paste(names(x$suggested_args), x$suggested_args,
                      sep = " = ", collapse = ", ")))
  } else {
    cat("  Ready-to-use args: none (see note) -- anchor only.\n")
  }
  invisible(x)
}
