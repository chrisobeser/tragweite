# Tests of the detection certificate. The expensive part is the
# dyadic sampler -- it is mocked here, together with the availability
# guard, because the public windkanal 0.3.0 does not ship it at all.
# What is then checked is exactly what can be wrong without a sampler:
# the cut arithmetic, the two false-alarm numbers, and the behaviour
# when fits fail. The sampler physics itself belongs to windkanal.

# The documentation text of tragweite(): in the source tree the Rd
# file, under R CMD check the installed Rd database.
doc_text <- function() {
  path <- test_path("..", "..", "man", "tragweite.Rd")
  if (file.exists(path)) {
    return(paste(readLines(path, warn = FALSE), collapse = " "))
  }
  db <- tryCatch(tools::Rd_db("tragweite"), error = function(e) NULL)
  if (is.null(db[["tragweite.Rd"]])) return(NULL)
  paste(utils::capture.output(print(db[["tragweite.Rd"]])),
        collapse = " ")
}

test_that("the in-sample FA is alpha by construction", {
  # cut* is the smallest cut with FA <= alpha, and fa_in_sample is the
  # FA exactly there -- so always floor(alpha*n)/n, whatever the null
  # distribution looks like. The number cannot be anything else and
  # therefore evidences nothing.
  set.seed(20260811)
  for (n in c(60, 200)) {
    for (a in c(0.05, 0.10)) {
      for (draw in 1:3) {
        p <- stats::runif(n)
        expect_equal(mean(p > .cert_cut(p, a)), floor(a * n) / n)
      }
    }
  }
})

test_that("the hold-out FA does not follow alpha -- it measures", {
  # 12 null values, ascending: cut on the first half (0.01-0.06),
  # measured on the second (0.07-0.12). By hand: the FA target on A is
  # floor(.25*6) = 1 -> cut_A = 0.05; all 6 B values lie above it.
  p <- seq(0.01, 0.12, by = 0.01)
  ho <- .cert_holdout(p, 0.25)
  expect_equal(ho$cut, 0.05)
  expect_equal(ho$fa, 1)
  expect_equal(ho$n_calib, 6L)
  expect_equal(ho$n_test, 6L)
  expect_equal(ho$mcse, 0)
  # in-sample on the same vector: exactly alpha (the tautology)
  expect_equal(mean(p > .cert_cut(p, 0.25)), 0.25)
  # too few worlds: no number rather than an invented one
  expect_true(is.na(.cert_holdout(c(0.1, 0.2, 0.3), 0.05)$fa))
})

test_that("the cut rule warns where a single world fixes it", {
  set.seed(4242)
  expect_warning(.cert_cut(stats::runif(10), 0.05), "SINGLE world")
  expect_warning(.cert_cut(rep(0, 60), 0.05), "degenerate null")
  expect_warning(.cert_cut(rep(1, 60), 0.05), "degenerate null")
  expect_equal(suppressWarnings(.cert_cut(rep(0, 60), 0.05)), 0)
  expect_equal(suppressWarnings(.cert_cut(rep(1, 60), 0.05)), 1)
  # positive case: enough worlds, real distribution -> no warning
  expect_silent(.cert_cut(stats::runif(200), 0.05))
  expect_error(.cert_cut(c(NA_real_, NA_real_), 0.05), "no valid")
})

test_that("fixed is checked and filled in, not passed through", {
  expect_error(ceiling_certificate(n_therapists = 4, caseload = 3,
                                   reps = 4, fixed = list(taux = 0.2),
                                   seed = 1, quiet = TRUE),
               "unknown: taux")
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) 0.5)
  z <- suppressWarnings(
    ceiling_certificate(n_therapists = 4, caseload = 3, n_sessions = 2,
                        reps = 6, effect_grid = 0.4,
                        fixed = list(tau = 0.5), seed = 1,
                        quiet = TRUE))
  expect_equal(z$fixed$tau, 0.5)     # passed in
  expect_equal(z$fixed$tau_x, 0.2)   # filled in from the default
  expect_equal(z$fixed$tau_c, 0.2)
})

test_that("the certificate reports both FA numbers honestly", {
  calls <- 0
  values <- c(seq(0.01, 0.12, by = 0.01), rep(0.5, 12))
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) {
    calls <<- calls + 1
    values[calls]
  })
  z <- ceiling_certificate(n_therapists = 4, caseload = 3,
                           n_sessions = 2, reps = 12, effect_grid = 0.4,
                           alpha = 0.25, seed = 1, quiet = TRUE)
  expect_equal(z$rule$cut, 0.09)
  expect_equal(z$rule$fa_in_sample, 0.25)   # = alpha, tautological
  expect_equal(z$rule$fa_holdout, 1)        # measures instead of following
  expect_equal(z$rule$holdout_n_calib, 6L)
  expect_equal(z$grid$power, 1)
  expect_equal(z$grid$n_valid, 12L)
  txt <- paste(capture.output(print(z)), collapse = " ")
  expect_match(txt, "in-sample 0.250", fixed = TRUE)
  expect_match(txt, "by construction", fixed = TRUE)
  expect_match(txt, "hold-out 1.000", fixed = TRUE)
  expect_match(txt, "tau = 0.3, tau_x = 0.2, tau_c = 0.2", fixed = TRUE)
})

test_that("failed fits do not throw the certificate away", {
  calls <- 0
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) {
    calls <<- calls + 1
    if (calls %% 5 == 0) stop("MCMC failure")
    0.1 + calls / 100
  })
  expect_warning(
    z <- ceiling_certificate(n_therapists = 4, caseload = 3,
                             n_sessions = 2, reps = 10,
                             effect_grid = 0.4, alpha = 0.25, seed = 1,
                             quiet = TRUE),
    "world fits failed")
  expect_equal(z$failures$calibration, 2)
  expect_equal(z$grid$n_failed, 2L)
  expect_equal(z$grid$n_valid, 8L)
  expect_true(is.finite(z$grid$power))
  expect_match(z$failures$first_error, "MCMC failure")
  # the MCSE uses the valid worlds, not reps
  expect_equal(z$grid$mcse,
               sqrt(z$grid$power * (1 - z$grid$power) / 8))
  expect_match(paste(capture.output(print(z)), collapse = " "),
               "4 world fits failed", fixed = TRUE)
})

test_that("a total calibration failure stops the run loudly", {
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) {
    stop("MCMC failure")
  })
  err <- tryCatch(
    ceiling_certificate(n_therapists = 4, caseload = 3, n_sessions = 2,
                        reps = 6, effect_grid = 0.4, seed = 1,
                        quiet = TRUE),
    error = function(e) conditionMessage(e))
  expect_match(err, "0 of 6 valid null worlds")
  expect_match(err, "MCMC failure")            # the reason travels with it
})

test_that("too few worlds stop without an empty reason", {
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) 0.5)
  err <- tryCatch(
    ceiling_certificate(n_therapists = 4, caseload = 3, n_sessions = 2,
                        reps = 1, effect_grid = 0.4, seed = 1,
                        quiet = TRUE),
    error = function(e) conditionMessage(e))
  expect_match(err, "at least 2 required")
  expect_false(grepl("First error", err, fixed = TRUE))
})

test_that("the certificate restores the RNG state too", {
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.duet_p_sd = function(s, delta_sd, ...) 0.5)
  set.seed(2026)
  before <- stats::runif(1)
  set.seed(2026)
  suppressWarnings(
    ceiling_certificate(n_therapists = 4, caseload = 3, n_sessions = 2,
                        reps = 6, effect_grid = 0.4, seed = 9,
                        quiet = TRUE))
  expect_identical(stats::runif(1), before)
})

test_that("no reference to an argument that does not exist", {
  expect_false("recalibrate" %in% names(formals(ceiling_certificate)))
  src <- doc_text()
  skip_if(is.null(src), "documentation source not found in the test tree")
  expect_false(grepl("recalibrate = TRUE", src, fixed = TRUE))
  expect_match(src, "ceiling_certificate", fixed = TRUE)
})
