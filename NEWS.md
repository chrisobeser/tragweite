# tragweite (development version)

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
