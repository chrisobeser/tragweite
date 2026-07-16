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

test_that("matching null is calibratable: inert measured C, no silent NaN", {
  # With tau_c = tau_xc = 0 the generator draws no therapist
  # attribute; the triple test would degenerate against a constant.
  # tragweite() attaches an inert measured attribute instead, which
  # is what a null world looks like from the analyst's chair.
  r <- tragweite("matching", tau = 0.3, tau_c = 0, tau_xc = 0,
                 n_therapists = 10, caseload = 5, n_sessions = 2,
                 reps = 30, seed = 500)
  expect_equal(r$n_failed, 0)
  expect_true(is.finite(r$power))
  expect_lt(r$power, 0.25)
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

test_that("design_ceiling: structure, monotonicity, plausible ceiling", {
  d <- design_ceiling("ate", effect_grid = c(0.2, 0.8),
                      n_therapists = 12, caseload = 8,
                      fixed = list(tau = 0, tau_x = 0, tau_c = 0,
                                   tau_xc = 0),
                      reps = 40, seed = 3, quiet = TRUE)
  expect_s3_class(d, "tragweite_ceiling")
  expect_equal(nrow(d$grid), 2L)
  expect_true(is.na(d$ceiling) || d$ceiling %in% c(0.2, 0.8))
  expect_gt(d$grid$power[2], d$grid$power[1])
  expect_error(design_ceiling("ate", n_therapists = 5, caseload = 5),
               "mandatory")
  expect_output(print(d), "Ceiling")
})
