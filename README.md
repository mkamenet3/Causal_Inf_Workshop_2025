# Causal_Inf_Workshop_2025
This repository contains materials for the Causal Inference Workshop presented at the February 2025 NCI Biostatistics Branch Workshop series, presented by Drs. Alex Keil and Maria Kamenetsky.

The contents of this repository include:

- **cross_sectional_gformula.html**: RMarkdown document containing the R code needed to conduct a cross-sectional analysis using the parametric g-formula. Functions are included within the Rmd document. The Rmd document contains the R code.
- **longitudinal_gformula.html**: RMarkdown document containing the R code needed to conduct a longitudinal analysis using the parametric g-formula. Functions are included within the Rmd document. The Rmd document contains the R code.
- **data/coalplant.csv**: data set used for the cross-sectional g-formula data example.
- **data/miners.csv**: data set used for the longitudinal g-formula data example, as well as the inverse probability weighting data example.

Files with HTML slides:

- **slides_cross_sectional_gformula.html**: Slides HTML document. The Rmd document contains the R code for the slides
- **slides_longitudinal_gformula.html**: Slides HTML document. The Rmd document contains the R code for the slides
  
Finally, though not covered in slides, we have included a tutorial on *inverse probability weighting (IPW)* as well. This tutorial uses the `miners.csv` data set. These files are ***ipw.html** (and respective Rmd file). 

*To render the HTML files, it is recommended you download them locally onto your computer as they will not render on GitHub*


**Installation**

Materials for the workshop were developed using `R` version 4.3.1 ("Beagle Scouts"). The latest version of `R` can be downloaded [here](https://cran.r-project.org/). We also recommend using RStudio, which can be downloaded [here](https://posit.co/downloads/).

The packages required for this workshop are:

- `ggplot2`
- `dplyr`
- `survival`
- `boot`
- `ggcorrplot`
- `Hmisc`
- `zoo`

To install a package, use `install.packages("packagename")`.




