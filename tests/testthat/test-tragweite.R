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
