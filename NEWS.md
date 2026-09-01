# tragweite 0.2.0 (2026-09-01)

Second feature release. It adds a preregistrable detection
certificate for pairing designs, four new levers on the simulation
core, and a round of hardening fixes across the package.

## New function

- `ceiling_certificate()` -- a preregistrable detection certificate
  for a randomized pairing design, in two steps: the rule's cut is
  **recalibrated on your own world stream** (null worlds of the same
  design, strict cut convention `cut* = min{k : FA(k) <= alpha}`),
  then power is measured across the effect grid with that cut.
  Recalibration is mandatory, not a switch: the frozen amplitude rule
  is stream-sensitive (its false-alarm rate has been observed to move
  from .06 to .15 across world streams), so a certificate carrying a
  cut over from elsewhere would be self-deception. **This function
  requires a windkanal sampler that is not public yet -- see the
  caveat below.**

## Extensions to existing functions

- `tragweite()` and `design_ceiling()` pass four new levers through to
  the generator and to the analysis:
  - `tau_xc_form` / `window_delta`: functional form of the dyadic
    term.
  - `rel_x` / `rel_c`: reliability of the OBSERVED patient /
    therapist attribute. Below 1 the analysis sees an attenuated
    version from its own seed stream while the truth in the generator
    is untouched -- measurement error on both sides, which attenuates
    a matching effect twice.
  - `detector`: `"lmer"` (default) or `"duet_s"`, an amplitude rule on
    a dyadic sampler's posterior, for `effect = "matching"` only.
  - `duet_rule`: the frozen rule `c(delta_sd, cut_sd)`.
- `design_ceiling()` passes all four through; `design_grid()` and
  `min_design()` do **not**, and say so in an `@section Limitations`
  in `?design_grid` and `?min_design`.

## Caveat: the duet_s detector is not usable against windkanal 0.3.0

`detector = "duet_s"` and `ceiling_certificate()` call
`windkanal::fit_cate_dyade_v2s()`. The public windkanal 0.3.0 ships
`fit_cate_dyade_v2p` but **not** `fit_cate_dyade_v2s`; that sampler is
planned for windkanal 0.4. Both entry points therefore stop
immediately with an explanatory error naming the missing function and
the version that will provide it, instead of failing deep inside a
world loop. This is stated in `?tragweite`, in
`?ceiling_certificate`, and in the README. Everything else in the
package -- including `design_ceiling()` with its default
`detector = "lmer"` -- runs against the public windkanal 0.3.0 in
full, and the test suite verifies both the working path and the guard.

## Fixes and hardening

- `alpha` is validated in every entry point (`tragweite()`,
  `design_grid()`, `min_design()`, `design_ceiling()`,
  `ceiling_certificate()`): a single number in (0, 1). The typo
  `alpha = 5` used to return power 1.000 with a positive verdict at a
  TRUE effect of zero, silently. `reps` is checked the same way, which
  also fixes the empty-reason error message at `reps = 0`.
- `design_ceiling()` no longer prints "effects below X were never
  observable". The ceiling is a statement about the 80% convention,
  not about possibility, and the verdict now names the **measured**
  detection probability below the ceiling. A third branch covers the
  ceiling sitting at the edge of the grid, and the
  ceiling-not-reached branch names the power at the top of the grid
  instead of calling the null "uninformative". The detection-ceiling
  vignette carried the same retired claim and has been brought in
  line.
- The detection-ceiling vignette's session-count comparison has been
  corrected. It asserted that measurement moves the ceiling further
  than recruitment does; its own output does not show that. At
  `reps = 60` the 4-session and 8-session designs land on the same
  ceiling (0.40) and differ by less than a combined Monte Carlo
  standard error at every grid point, with the sign of the difference
  changing across cells. The section now reports that result and
  explains what resolution would be needed to settle the question,
  instead of asserting a pattern the run does not support.
- `ceiling_certificate()` reports an out-of-sample false-alarm rate:
  the cut is calibrated on the first half of the null worlds and
  measured on the second. The in-sample number is
  `floor(alpha * reps)/reps` by construction and is now labeled as
  such -- it cannot be anything else and therefore evidences nothing.
  Failed fits are counted per grid point instead of aborting a long
  run, and `fixed` is validated and filled in from the default.
- `detector = "duet_s"` no longer carries an `alpha` that does not
  act: a non-default `alpha` warns, the returned `alpha` is `NA`, and
  the print states the rule's measured false-alarm rate (.083
  calibration stream, .067 fresh stream; standard design, seed 4242)
  rather than implying a nominal .05.
- The global RNG state is restored on exit from `tragweite()` and
  `ceiling_certificate()`, so a call no longer moves the caller's own
  stream. The fixed world (`fixed`) travels in the return object and
  into the print of `design_ceiling()` and `ceiling_certificate()`.
- The sampler seed behind `duet_s` stays constant (common random
  numbers on the estimator side) but is now a named argument of the
  internal `.duet_p_sd()` rather than buried in the call.
- Documentation ghost removed: the help page recommended a
  `recalibrate = TRUE` argument that does not exist and never did.
- Tests: 31 blocks, 135 assertions, 0 failures against the public
  windkanal 0.3.0. (Two blocks exercise the v2S guard and are skipped
  automatically on an installation that does ship the v2S sampler.)

## Packaging

- `DESCRIPTION` requires `windkanal (>= 0.3.0)` and adds `utils`;
  `Remotes` points at the public windkanal repository.
- `URL` and `BugReports` fields added.
- The `scenario()` literature anchors (Nye et al. 2023, Constantino
  et al. 2021) are unchanged and remain what they were: anchors with
  an estimand note and **no** automatic translation.

# tragweite 0.1.0

First feature release.

- New: `design_ceiling()`, retrospective design diagnosis. Given the
  structure of an existing study (therapists, caseload, ICC,
  sessions), it walks an ascending effect grid under common random
  numbers and returns the smallest effect the design detects at a
  target power -- its detection ceiling. Named `design_ceiling`, not
  `ceiling`, so that it does not mask `base::ceiling()`.
- New vignette "The detection ceiling: what a design could ever have
  shown", including a session-count comparison (corrected in 0.2.0,
  see above).
- Fixed: matching null worlds are now calibratable. With
  `tau_c = tau_xc = 0` the generator draws no therapist attribute,
  so the triple test degenerated against a constant, every fit
  failed, and the failure was silently swallowed into a `NaN` power.
  `tragweite()` now attaches an inert *measured* therapist attribute
  from its own seed stream (the null world an analyst actually
  faces), fit failures carry the first error message as a warning,
  and an all-worlds-failed run stops loudly.

# tragweite 0.0.0.9000 (initial development version)

- Core working: `tragweite()` (simulation-based power for ATE,
  moderation, and matching contrasts under therapist nesting, with
  Satterthwaite detection criteria), `design_grid()` (common random
  numbers across designs), `min_design()`, `scenario()` (anchored
  scenario library that refuses undefensible translations).
- Validation: exact-seed reproduction of the windkanal
  reference (power 0.6300, 300/300 decisions identical);
  fresh-seed cross-check within Monte Carlo error.
- CI: R CMD check on every push (windkanal via remotes).
- Hex logo.
