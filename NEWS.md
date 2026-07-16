# tragweite 0.1.0

First feature release.

- New: `design_ceiling()`, retrospective design diagnosis. Given the
  structure of an existing study (therapists, caseload, ICC,
  sessions), it walks an ascending effect grid under common random
  numbers and returns the smallest effect the design detects at a
  target power -- its detection ceiling. Effects below that line
  were never observable: a null from such a design is a fact about
  the design, not about the world.
- New vignette "The detection ceiling: what a design could ever have
  shown", including the measurement-beats-recruitment comparison.
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
- Validation: exact-seed reproduction of the audited windkanal
  hardening reference (power 0.6300, 300/300 decisions identical);
  fresh-seed cross-check within Monte Carlo error.
- CI: R CMD check on every push (windkanal via remotes).
- Hex logo.
