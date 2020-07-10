---
layout: page
title: Updates
---

The most recent package updates:

## pavo 2.5.0 (development)

### NEW FEATURES AND SIGNIFICANT CHANGES

* Add ability to compute colour volume by using alphashapes instead of convex 
  hulls. The functions `vol()`, `tcsvol()` and `voloverlap()` gain a new 
  argument `type = c("convex", "alpha")` to decide how you want to compute the
  colour volume. Please refer to the vignette `vignette("pavo-5-alphashapes", 
  package = "pavo")` for more information. As a result, the argument order
  in these 3 function has changed. Check the documentation to update your
  scripts accordingly. The function `summary.colspace()` also gains an 
  additional column that returns that colour volume computed with an alpha-shape
  of parameter alpha* in the case of `tcs` objects.

### MINOR FEATURES AND BUG FIXES

* Maximum quantum catches computation (`data.maxqcatches` attribute) now works
for segment "visual model" as well. As a side effect, this removes a warning
that occurred when users ran `vismodel(..., visual = "segment")`.
* `sensmodel()` now accepts the argument `sensnames`, for specifying the names
of the resulting sensitivity curves on-the-fly.
* CIE models now accept data created outside of `vismodel()`, by allowing users to 
specify the illuminant and viewer sensitivity function used when estimating XYZ values 
(via `illum` and `visual` arguments in `colspace()`).
* `bootcoldist()` is now laxer in its argument checks and accept objects that 
are neither `vismodel` or `colspace` objects. This means you can now use this 
function on quantum catches dataframe that you obtained outside of pavo, such
as the MICA toolbox.
* `summary.colspace()` now prints a more explicit error when the `by` argument
value is not a multiple of the number of rows in the colspace object (i.e., the
number of spectra)
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

## lightr 1.3 (development)

### Minor changes

* `spec_ID` extraction from Avantes exported files (`ttt` and `trt`) is now
more robust, meaning it should work for more files.

## lightr 1.2 (24/06/2020)

## Minor changes

* fixed tests on platform with no long-doubles ('noLD') 
* restored tests on 32bits machines
* `spec_ID` extraction from Avantes exported files (`ttt` and `trt`) is now
more robust, meaning it should work for more files.