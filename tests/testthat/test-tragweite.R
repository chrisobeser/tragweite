test_that("null design: false-positive rate near alpha", {
  r <- tragweite("ate", tau = 0, n_therapists = 10, caseload = 5,
                 reps = 60, seed = 100)
  expect_lt(r$power, 0.15)
  expect_equal(r$n_valid + r$n_failed, 60L)
})

test_that("strong ATE in decent design: high power", {
  r <- tragweite("ate", tau = 0.8, n_therapists = 20, caseload = 10,
                 reps = 40, seed = 200)
  expect_gt(r$power, 0.8)
})

test_that("deterministic under identical seed", {
  a <- tragweite("moderation", tau = 0.5, tau_x = 0.4,
                 n_therapists = 10, caseload = 8, reps = 20, seed = 7)
  b <- tragweite("moderation", tau = 0.5, tau_x = 0.4,
                 n_therapists = 10, caseload = 8, reps = 20, seed = 7)
  expect_identical(a$power, b$power)
})

test_that("matching effect uses the triple interaction and detects a large tau_xc", {
  r <- tragweite("matching", tau = 0.3, tau_x = 0.2, tau_c = 0.2,
                 tau_xc = 0.9, n_therapists = 30, caseload = 10,
                 reps = 25, seed = 300)
  expect_gt(r$power, 0.5)
})

test_that("matching null is calibratable: inert measured C, no silent NaN", {
  # With tau_c = tau_xc = 0 the generator draws no therapist
  # attribute; the triple test would degenerate against a constant,
  # every fit would fail by subscript error, and the error used to be
  # swallowed into a NaN power. Now: an inert measured attribute, all
  # fits valid, false-alarm rate near alpha.
  r <- tragweite("matching", tau = 0.3, tau_c = 0, tau_xc = 0,
                 n_therapists = 10, caseload = 5, n_sessions = 2,
                 reps = 30, seed = 500)
  expect_equal(r$n_failed, 0)
  expect_true(is.finite(r$power))
  expect_lt(r$power, 0.25)                    # null discipline (rough)
  # deterministic: the inert C stream depends only on the seed
  r2 <- tragweite("matching", tau = 0.3, tau_c = 0, tau_xc = 0,
                  n_therapists = 10, caseload = 5, n_sessions = 2,
                  reps = 30, seed = 500)
  expect_identical(r$power, r2$power)
})

test_that(".detect aborts loudly when C is missing instead of testing c = 0", {
  s <- windkanal::sim_stream(n_therapists = 4, patients_per_therapist = 3,
                             n_sessions = 2, tau = 0.3, seed = 1)
  expect_error(.detect(s, "matching", 0.05), "therapist attribute")
})

test_that("min_design returns smallest hitting design", {
  m <- min_design("ate", tau = 0.8, n_therapists = c(6, 20),
                  caseload = c(5, 10), reps = 30, seed = 400,
                  quiet = TRUE)
  expect_false(is.null(m$design))
  expect_true(m$design$power >= 0.8)
})

test_that("scenario registry is honest about missing translations", {
  s <- scenario("personalization_floor")
  expect_null(s$suggested_args)
  expect_match(s$estimand_note, "NOT a moderation amplitude")
  f <- scenario("field_realistic")
  expect_equal(f$suggested_args$tau, 0.2)
})

test_that("design_ceiling: structure, monotonicity, plausible ceiling", {
  d <- design_ceiling("ate", effect_grid = c(0.2, 0.8),
                      n_therapists = 12, caseload = 8,
                      fixed = list(tau = 0, tau_x = 0, tau_c = 0,
                                   tau_xc = 0),
                      reps = 40, seed = 3, quiet = TRUE)
  expect_s3_class(d, "tragweite_ceiling")
  expect_equal(nrow(d$grid), 2L)
  expect_true(is.na(d$ceiling) || d$ceiling %in% c(0.2, 0.8))
  expect_gt(d$grid$power[2], d$grid$power[1])  # more effect, more power
  expect_error(design_ceiling("ate", n_therapists = 5, caseload = 5),
               "mandatory")
  expect_output(print(d), "Ceiling")
})

# ---- Argument guards ---------------------------------------------

test_that("alpha is validated -- the typo alpha = 5 stops the run", {
  # Before: alpha = 5 returned power 1.000 with a positive verdict at
  # a TRUE tau of 0, and alpha = -1 returned 0.000 -- both silently.
  expect_error(tragweite("ate", tau = 0, n_therapists = 10, caseload = 5,
                         n_sessions = 2, reps = 20, alpha = 5, seed = 1),
               "alpha")
  expect_error(tragweite("ate", n_therapists = 6, caseload = 4,
                         reps = 5, alpha = -1, seed = 1), "alpha")
  expect_error(tragweite("ate", n_therapists = 6, caseload = 4,
                         reps = 5, alpha = c(0.05, 0.1), seed = 1),
               "alpha")
  expect_error(tragweite("ate", n_therapists = 6, caseload = 4,
                         reps = 5, alpha = NA, seed = 1), "alpha")
  expect_error(tragweite("ate", n_therapists = 6, caseload = 4,
                         reps = 5, alpha = "0.05", seed = 1), "alpha")
  # the same guard in every wrapper that passes alpha through
  expect_error(design_grid("ate", n_therapists = 6, caseload = 4,
                           reps = 5, alpha = 5, seed = 1, quiet = TRUE),
               "alpha")
  expect_error(min_design("ate", n_therapists = 6, caseload = 4,
                          reps = 5, alpha = 5, seed = 1, quiet = TRUE),
               "alpha")
  expect_error(design_ceiling("ate", n_therapists = 6, caseload = 4,
                              reps = 5, alpha = 5, seed = 1,
                              quiet = TRUE), "alpha")
  expect_error(ceiling_certificate(n_therapists = 6, caseload = 4,
                                   reps = 5, alpha = 5, seed = 1,
                                   quiet = TRUE), "alpha")
  # positive case: a valid non-default alpha runs through
  r <- tragweite("ate", tau = 0.8, n_therapists = 8, caseload = 5,
                 n_sessions = 2, reps = 5, alpha = 0.10, seed = 11)
  expect_equal(r$alpha, 0.10)
})

test_that("reps = 0 names the reason instead of leaving it empty", {
  # Before: "All 0 worlds failed to fit -- first error: "
  err <- tryCatch(tragweite("ate", n_therapists = 6, caseload = 4,
                            reps = 0, seed = 1),
                  error = function(e) conditionMessage(e))
  expect_match(err, "reps")
  expect_false(grepl("first error", err, fixed = TRUE))
})

# ---- The ceiling verdict claims a convention, not impossibility ----

test_that("the ceiling verdict does not claim unobservability", {
  # A real diagnosed case rebuilt as an object: ceiling .45, but .670
  # measured detection at effect .40 -- "never observable" was
  # refutable from the function's own table.
  g <- data.frame(
    effect = seq(0.10, 0.60, by = 0.05),
    power = c(.080, .150, .240, .355, .450, .555, .670, .815, .885,
              .925, .975),
    mcse = c(.0192, .0252, .0302, .0338, .0352, .0351, .0332, .0275,
             .0226, .0186, .0110))
  x <- structure(
    list(ceiling = 0.45, grid = g, target = 0.80, effect = "ate",
         criterion = "Satterthwaite", fixed = list(tau = 0),
         design = list(n_therapists = 20, caseload = 10, N = 200,
                       n_sessions = 4, icc = 0.10,
                       z_level = "patient")),
    class = "tragweite_ceiling")
  txt <- paste(capture.output(print(x)), collapse = " ")
  expect_false(grepl("never observable", txt, fixed = TRUE))
  expect_match(txt, "fall short of the 80% convention", fixed = TRUE)
  expect_match(txt, "0.67", fixed = TRUE)   # the measured rate below
  expect_match(txt, "0.40", fixed = TRUE)
  # NOT-REACHED branch: no blanket "uninformative" claim
  y <- x; y$ceiling <- NA_real_; y$grid <- g[1:6, ]
  txt_na <- paste(capture.output(print(y)), collapse = " ")
  expect_false(grepl("uninformative", txt_na, fixed = TRUE))
  expect_match(txt_na, "reaches only 0.56", fixed = TRUE)
  # ceiling at the lowest grid point: no claim about the unmeasured
  z <- x; z$ceiling <- 0.10
  expect_match(paste(capture.output(print(z)), collapse = " "),
               "not examined", fixed = TRUE)
})

test_that("a real run prints the honest verdict too", {
  d <- design_ceiling("ate", effect_grid = c(0.2, 0.8),
                      n_therapists = 8, caseload = 5, n_sessions = 2,
                      fixed = list(tau = 0, tau_x = 0, tau_c = 0,
                                   tau_xc = 0),
                      reps = 10, seed = 3, quiet = TRUE)
  txt <- paste(capture.output(print(d)), collapse = " ")
  expect_false(grepl("never observable", txt, fixed = TRUE))
  # the fixed world belongs in the output, not only in the object --
  # the varied parameter (here tau) lives in the grid, not there
  expect_match(txt, "tau_x = 0, tau_c = 0, tau_xc = 0", fixed = TRUE)
  expect_false(grepl("Fixed     : tau = ", txt, fixed = TRUE))
  expect_match(txt, "holds for THIS world", fixed = TRUE)
})

test_that("a partial fixed falls back on the DOCUMENTED default", {
  # Before, a partial `fixed` silently fell through to tragweite()'s
  # own tau = 0.5, although the documentation promises tau = 0.3.
  d <- design_ceiling("moderation", effect_grid = 0.5, n_therapists = 6,
                      caseload = 4, n_sessions = 2,
                      fixed = list(tau_x = 0), reps = 4, seed = 21,
                      quiet = TRUE)
  expect_equal(d$fixed$tau, 0.3)
  expect_equal(d$fixed$tau_c, 0)
  expect_equal(d$fixed$tau_xc, 0)
  # a complete fixed is left untouched
  v <- design_ceiling("moderation", effect_grid = 0.5, n_therapists = 6,
                      caseload = 4, n_sessions = 2,
                      fixed = list(tau = 0, tau_x = 0, tau_c = 0,
                                   tau_xc = 0),
                      reps = 4, seed = 21, quiet = TRUE)
  expect_equal(v$fixed$tau, 0)
})

test_that("the global RNG state survives a run", {
  set.seed(999)
  before <- runif(1)
  set.seed(999)
  invisible(tragweite("ate", tau = 0.3, n_therapists = 6, caseload = 4,
                      n_sessions = 2, reps = 3, seed = 500))
  expect_identical(runif(1), before)
  # negative case: where there was NO state, none is left behind
  if (exists(".Random.seed", envir = globalenv())) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(tragweite("ate", tau = 0.3, n_therapists = 6, caseload = 4,
                      n_sessions = 2, reps = 3, seed = 500))
  expect_false(exists(".Random.seed", envir = globalenv()))
  set.seed(1)
})

# ---- The v2S detector guard (windkanal 0.4) -----------------------

test_that("duet_s stops honestly while windkanal ships no v2S sampler", {
  # The public windkanal 0.3.0 exports fit_cate_dyade_v2p but not
  # fit_cate_dyade_v2s. The package must say so at the point of call
  # rather than fail deep inside a world loop.
  skip_if(.duet_s_available(),
          "windkanal ships fit_cate_dyade_v2s -- guard not applicable")
  expect_error(
    tragweite("matching", tau = 0.3, tau_x = 0.2, tau_c = 0.2,
              tau_xc = 0.3, n_therapists = 4, caseload = 3,
              n_sessions = 2, reps = 2, detector = "duet_s", seed = 1),
    "fit_cate_dyade_v2s")
  err <- tryCatch(
    tragweite("matching", tau_xc = 0.3, n_therapists = 4, caseload = 3,
              reps = 2, detector = "duet_s", seed = 1),
    error = function(e) conditionMessage(e))
  expect_match(err, "windkanal 0.4", fixed = TRUE)
  expect_match(err, "windkanal 0.3.0 does not ship it", fixed = TRUE)
  expect_match(err, "detector = 'lmer'", fixed = TRUE)
  # The certificate rests on the same sampler and stops the same way,
  # but must NOT quote `detector` at the user: ceiling_certificate()
  # has no such argument, so that remedy would be unfollowable.
  cert_err <- tryCatch(
    ceiling_certificate(n_therapists = 4, caseload = 3, reps = 4,
                        seed = 1, quiet = TRUE),
    error = function(e) conditionMessage(e))
  expect_match(cert_err, "fit_cate_dyade_v2s", fixed = TRUE)
  expect_match(cert_err, "ceiling_certificate() requires", fixed = TRUE)
  expect_match(cert_err, "use design_ceiling()", fixed = TRUE)
  expect_false(grepl("detector", cert_err, fixed = TRUE))
  # and the guard is reached through design_ceiling() as well
  expect_error(design_ceiling("matching", effect_grid = 0.3,
                              n_therapists = 4, caseload = 3,
                              detector = "duet_s", reps = 2, seed = 1,
                              quiet = TRUE),
               "fit_cate_dyade_v2s")
})

test_that("the guard fires before any simulation work", {
  skip_if(.duet_s_available(),
          "windkanal ships fit_cate_dyade_v2s -- guard not applicable")
  # reps = 1e6 would take hours if the guard were placed after the
  # world loop; it must return immediately.
  t0 <- Sys.time()
  expect_error(
    tragweite("matching", tau_xc = 0.3, n_therapists = 40,
              caseload = 20, reps = 1e6, detector = "duet_s", seed = 1),
    "fit_cate_dyade_v2s")
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 5)
})

test_that("duet_s is refused for non-matching effects before the guard", {
  expect_error(
    tragweite("ate", n_therapists = 4, caseload = 3, reps = 2,
              detector = "duet_s", seed = 1),
    "built for effect")
})

test_that("the sampler seed is constant and visible in the argument", {
  # Deliberately NOT tied to the user seed: the sampler stream stays
  # fixed across worlds so that a power difference is a design
  # difference. Conservative = behaviour unchanged, but the choice is
  # named.
  expect_equal(eval(formals(.duet_p_sd)$mcmc_seed), 1)
  expect_true(all(c("mcmc_seed", "nburn", "nsim") %in%
                    names(formals(.duet_p_sd))))
})

test_that("alpha does not steer the duet_s detector, and it says so", {
  local_mocked_bindings(.require_duet_s = function(...) invisible(TRUE))
  local_mocked_bindings(.detect_duet_s = function(s, duet_rule) TRUE)
  expect_warning(
    r <- tragweite("matching", tau = 0.3, tau_x = 0.2, tau_c = 0.2,
                   tau_xc = 0.3, n_therapists = 4, caseload = 3,
                   n_sessions = 2, reps = 2, alpha = 0.01,
                   detector = "duet_s", seed = 1),
    "alpha does not act")
  expect_true(is.na(r$alpha))        # alpha does not act -> no value
  txt <- paste(capture.output(print(r)), collapse = " ")
  expect_match(txt, "false-alarm rate is not nominal", fixed = TRUE)
  expect_match(txt, ".083", fixed = TRUE)
  expect_match(txt, ".067", fixed = TRUE)
  expect_match(txt, "P{sd(h) > 0.018} > 0.520", fixed = TRUE)
  # negative case: the lmer path carries its alpha and no FA note
  r2 <- tragweite("ate", tau = 0.5, n_therapists = 6, caseload = 4,
                  n_sessions = 2, reps = 3, alpha = 0.01, seed = 1)
  expect_equal(r2$alpha, 0.01)
  expect_null(r2$fa_note)
  expect_false(grepl("false-alarm",
                     paste(capture.output(print(r2)), collapse = " "),
                     fixed = TRUE))
})

test_that("design_grid: rows, CRN identity, attribute", {
  g <- design_grid("ate", tau = c(0.3, 0.8), n_therapists = c(6, 10),
                   caseload = 4, n_sessions = 2, reps = 8, seed = 77,
                   quiet = TRUE)
  expect_s3_class(g, "tragweite_grid")
  expect_equal(nrow(g), 4L)
  expect_equal(attr(g, "effect"), "ate")
  expect_true(all(g$n_failed == 0))
  # CRN: the same cell computed on its own gives EXACTLY the same power
  ref <- tragweite("ate", tau = 0.8, n_therapists = 10, caseload = 4,
                   n_sessions = 2, reps = 8, seed = 77)
  cell <- g$power[g$tau == 0.8 & g$n_therapists == 10]
  expect_identical(cell, ref$power)
  # negative case: a different cell is not the same number
  expect_false(identical(g$power[g$tau == 0.3 & g$n_therapists == 6],
                         ref$power))
})
