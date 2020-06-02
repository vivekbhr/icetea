# icetea

[![Build Status](https://travis-ci.org/vivekbhr/icetea.svg?branch=master)](https://travis-ci.org/vivekbhr/icetea)
[![Build status](https://ci.appveyor.com/api/projects/status/58od90th3qg5o7f0/branch/master?svg=true)](https://ci.appveyor.com/project/vivekbhr/icetea/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/vivekbhr/icetea/badge.svg?branch=master)](https://coveralls.io/github/vivekbhr/icetea?branch=master)

*An R package for analysis of data produced by transcript 5' profiling methods like RAMPAGE and MAPCap.*

![icetea_logo](https://raw.githubusercontent.com/vivekbhr/icetea/master/icetea_front.png)

**The icetea R package for analysis of TSS profiling data** allows users to processes data from multiplexed,
5’-profiling techniques such as RAMPAGE and the recently developed MAPCap protocol. TSS detection and
differential TSS analysis can be performed using replicates on any 5’-profiling dataset.
**Left panel :**  Typical analysis steps for MAPCap data that can be performed using icetea.
**Right panel :** Showing some of the quality control and visualization outputs from the package.
Proportion of sequencing reads used for each step (Top), comparison of TSS accuracy (w.r.t. annotated TSS)
between samples (middle), and MA-plots from differential TSS analysis (Bottom).

Additionally, analysis of RNA-binding protein locations via the [FLASH protocol](https://academic.oup.com/nar/article/48/3/e15/5658447) can also be preformed via icetea.


## Installing icetea

### Stable version

Stable release of icetea is available via [bioconductor](https://bioconductor.org/packages/release/bioc/html/icetea.html).

```{r}
## first install BiocManager
install.packages("BiocManager")

## then install icetea
BiocManager::install("icetea")
```
### Latest version
The latest (development) version of icetea can be installed using Devtools. Devtools can be installed from CRAN.

```{r}
## first install devtools
install.packages("devtools")

## then install icetea
devtools::install_github("vivekbhr/icetea")
```

## Documentation

Please visit the icetea [website](https://vivekbhr.github.io/icetea/) for the package documentation and vignette.

## Citation

If you use icetea in your analysis, please cite:

Bhardwaj, V., Semplicio, G., Erdogdu, N. U., Manke, T. & Akhtar, A. MAPCap allows high-resolution detection and differential expression analysis of transcription start sites. [Nat. Commun. 10, 3219 (2019)](https://www.nature.com/articles/s41467-019-11115-x)
