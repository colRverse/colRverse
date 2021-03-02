---
layout: page
title: Packages
---

# [pavo](http://pavo.colrverse.com)  
**Organize, visualise, and analyse spectral and spatial colour data** 

`pavo` offers a flexible and integrated workflow for working with spectral and spatial colour data. It includes functions that take advantage of new data classes to work seamlessly from importing raw spectra and images, to publication-quality visualisations, and analyses via a suite of analytical methods and visual models.  
[Learn more...](http://pavo.colrverse.com)

**Install:** `install.packages('pavo')`, or the development version, via `remotes::install_github("rmaia/pavo")`

**Cite:** Maia R, Gruson H, Endler JA, White TE (2019) pavo 2: new tools for the spectral and spatial analysis of colour in R.  [_Methods in Ecology and Evolution_](http://dx.doi.org/10.1111/2041-210X.13174) 10, 1097-1107. 

---

# [lightr](http://lightr.colrverse.com) 
**Import spectral data and metadata**  

`lightr` offers a unified, user-friendly interface for reading UV-VIS reflectance, transmittance, and/or absorbance spectral files and associated metadata from a suite of proprietary (and generally unfriendly) file formats, across all systems.  
[Learn more...](http://lightr.colrverse.com)

**Install:** `install.packages('lightr')` 

**Cite:** Gruson H, White TE, Maia R (2019) lightr: import spectral data and metadata in R. [_Journal of Open Source Software_](https://doi.org/10.21105/joss.01857). 

---

# [acuityview](https://doi.org/10.1111/2041-210X.12911) 
**Represent the effects of visual acuity**

`AcuityView` offers an intuitive method for representing scenes as they may appear to viewers with less acute visual systems, either as a valuable end unto itself or as part of a broader analysis.

**Install:** The `AcuityView` 2.0 algorithm is fully implemented in `pavo` and can be accessed via the `procimg()` function, while the original package remains available via CRAN (`install.packages('AcuityView')`). Either way, please remember to cite and read the accompanying publication for full details.

**Cite:** Caves E, Johnsen S (2018) AcuityView: An r package for portraying the effects of visual acuity on scenes observed by an animal. [_Methods in Ecology and Evolution_](https://doi.org/10.1111/2041-210X.12911) 9:793-797.
