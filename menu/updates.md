---
layout: page
title: Updates
---

The most recent package updates:

## pavo 2.5.0 (development)

### New features and significant changes

* Add ability to compute colour volume by using alphashapes instead of convex 
  hulls. The functions `vol()`, `tcsvol()` and `voloverlap()` gain a new 
  argument `type = c("convex", "alpha")` to decide how you want to compute the
  colour volume. Please refer to the vignette `vignette("pavo-5-alphashapes", 
  package = "pavo")` for more information. As a result, the argument order
  in these 3 function has changed. Check the documentation to update your
  scripts accordingly.

### Minor changes

* Added a continuous measure of hue to the output of the categorical model of Troje (1993)
* `teal` example dataset columns have been renamed to add an additional zero in
front of single digit numbers, so that column names now sort in the correct
order by default.
* Some very small negative values in the built-in visual system data have 
been corrected .
* `as.rspec()` is now more lenient for wavelength trimming when `interp = FALSE`
and now works even if the specified `lims` do not correspond to actual wl values
from the input object.

## pavo 2.4.0 (08/02/2020)

### NEW FEATURES AND SIGNIFICANT CHANGES

* Fixed a bug introduced in version 2.3.0 that gave wrong values for S1UV and 
S1V in `summary.rspec()`.
* several `plot()` functions for colspace (`triplot()`, `tcsplot()`, 
`tetraplot()`) gain a new `gamut` argument to plot the maximum gamut for a given
visual system and illuminant. `summary.colspace()` also now returns the maximum
colour volume for a given visual system and illuminant that you can use to
compare to the realised volume by a given dataset. More information in PR #180.
* parallel processing now relies on the `future` package, which offers windows 
and high performance computing (HPC) environments support. The progress bar is
produced by the `progressr` package and can be customised as well. As a 
consequence, the `cores` argument in `getspec()`, `adjacent()` and `classify()`
has been deprecated.

### MINOR FEATURES AND BUG FIXES

* fixed a plotting bug introduced in version 2.3.0 where it was required to run
`projplot()` twice for the background grid to be displayed.
* fixed a bug in `summary.colspace()` where `NULL` was returned instead of 
`summary.data.frame()` for non-tcs colourspaces.
* fix partial matching warnings in examples and in `bootcooldist()`
* the package has a new website at `pavo.colrverse.com`
* fixed a bug in `coldist()`that prevented the calculation of achromatic contrast
when using custom quantum catch data 

---

## lightr 1.2 (development)

### Minor changes

* `spec_ID` extraction from Avantes exported files (`ttt` and `trt`) is now
more robust, meaning it should work for more files.

## lightr 1.1 (02/05/2020)

### New features and major changes

* `date` column in metadata is now always formatted as ISO 8601.
* `lightr` can now import AvaSoft8 files (test files provided by M.D. Shawkey 
and L. Swierk), via the functions `lr_parse_rfl8()`/`lr_parse_raw8()`.
* `lightr` can correctly imports `TRM` files from AvaSoft 6.0 (previously it 
only supported files from AvaSoft 7.0).
* `lightr` can now import binary `.spc` files (via the `lr_parse_spc()` parser).
This format is used by OceanInsight and CRAIC.

### Minor changes

* new test suite on a different locale (in this case `fr_FR.UTF-8`) to ensure
parsing is locale-independent.
* warnings on CRAN build system for platforms that don't support markdown 2 
have been fixed.
* new, stricter tests for various file formats.
* `jdx` files saved in a locale that uses `,` as the decimal separator are now
parsed correctly.
* Avantes exported files in non-English locales (`ttt` and `trt` files) are now 
parsed correctly again (this was a regression compared to pavo's 
`pavo::getspec()`). Thanks to A. Fargevieille for reporting the issue and 
providing a test file.