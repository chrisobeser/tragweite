# Contributing to tragweite

Thank you for your interest. The package is young; contributions are welcome.

* **Bugs / questions:** open a GitHub issue with a minimal reproducible
  example (a `tragweite()` call with its seed is usually enough).
* **New scenarios:** every anchor needs a source (reference with DOI)
  and an estimand note that states what the published number is and
  whether a defensible mapping to data-generating parameters exists.
  Uncited anchors and silent translations are not accepted.
* **New detection criteria:** a criterion enters the package only with
  a validation reference: either a windkanal validation cell or an
  equivalent simulation study showing its false-alarm rate and power
  behavior under nesting.
* **Style:** match the existing code; comments state constraints, not
  narration. Run the test suite before submitting; CI must stay green.
* **Reproducibility discipline:** `seed` stays mandatory in every
  user-facing function, and grid comparisons keep common random
  numbers.
