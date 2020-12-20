---
layout: page
title: Updates
---

The most recent package updates:

# pavo 2.7.0 (development)

## MINOR FEATURES AND BUG FIXES

# pavo 2.6.1 (2020-12-21)

## MINOR FEATURES AND BUG FIXES

* Minor fix to a vignette to avoid installation issues. 
* `options()` and `par()` are now always set locally, including in vignettes and
examples, as to prevent spillover of these changes in the user session

# pavo 2.6.0

## MINOR FEATURES AND BUG FIXES

* `bootcoldist()` and `adjacent()` now use a random number generator that 
generates statistically sound values, even when ran in parallel. The output of
these functions is thus expected to slightly change, even if you set the seed
before.
* The alphashape3d package has been moved to `Suggests`, which means it will
not be installed automatically when you install pavo from CRAN and that you
will need to install it yourself if you need it (for `vol(type = "alpha")`,
`vol(type = "alpha")` and `voloverlap(type = "alpha")`)

---

# lightr 1.4 (development)

## Minor changes and bug fixes

* `lr_parse_generic()` now makes sure that the data is ordered by increasing 
wavelengths, which fixes a bug reported by @itamshab

# lightr 1.3 (2020-07-02)

## Minor changes

* (mostly internal) `compute_processed()` function is now named
`lr_compute_processed()`
* disable hash tests on Solaris (the output is still checked by other tests)