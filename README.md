# CausalInf_Workshop_2025
This repository contains materials for the Causal Inference Workshop presented at the 2025 NCI Biostatistics Branch Workshop series, presented by Drs. Alex Keil and Maria Kamenetsky.

The contents of this repository include:

- *cross_sectional_gformula.Rmd*: RMarkdown document containing the R code needed to conduct a cross-sectional analysis using the parametric g-formula. Functions are included within the Rmd document. The HTML document contains the rendered version.
- *longitudinal_gformula.Rmd*: RMarkdown document containing the R code needed to conduct a longitudinal analysis using the parametric g-formula. Functions are included within the Rmd document. The HTML document contains the rendered version.
- *data/coalplant.csv*: data set used for the cross-sectional g-formula data example.
- *data/miners.csv*: data set used for the longitudinal g-formula data example.

**Installation**

Materials for the workshop were developed using `R` version 4.3.1 ("Beagle Scouts"). The latest version of `R` can be downloaded [here](https://cran.r-project.org/). We also recommend using RStudio, which can be downloaded [here](https://posit.co/downloads/).

The packages required for this workshop are:

- `ggplot2`
- `dplyr`
- `survival`
- `boot`
- `ggcorrplot`

To install a package, use `install.packages("packagename")`.




