**NOTE: This package is not actively developed, nor is it complete.**

**For an actively developed (and more efficient) package for Bayesian estimation of non-Gaussian SVARs, see [bsvar](https://github.com/jetroant/bsvar).**

**This package ([rbsvar](https://github.com/jetroant/rbsvar)) implements all the methods discussed in [Anttonen, Lanne and Luoto (2023)](https://onlinelibrary.wiley.com/doi/full/10.1002/jae.3019).**


# rbsvar 

Bayesian estimation of statistically identified "robust" structural vectorautoregressive models. 

## Getting Started

### Prerequisites

Make sure you have the latest version of R installed on your computer. On top of that, see the operating system specific further prerequisities below for the development version of the package to work on your computer.


#### Windows: 
Rtools must be installed on your computer. See: https://cran.r-project.org/bin/windows/Rtools/

#### Mac: 
Make sure you have Xcode installed. If not, it can be found from the App Store OR it can be installed in the shell by: 

```
xcode-select --install
```

After installing Xcode, a few more steps might still be necessary. For comprehensive instructions, see: https://thecoatlessprofessor.com/programming/cpp/openmp-in-r-on-os-x/#after-3-4-0

In short, (i) the official gfortran 6.1 build (see: https://gcc.gnu.org/wiki/GFortranBinaries#MacOS-11) and (ii) clang (compiler, see: https://uofi.box.com/v/r-macos-clang-pkg) may need to be downloaded and installed.

#### Linux:
Everything should be just fine. Just make sure you have everything Rcpp needs.

### Installing the package

If you do not have devtools installed, install it in R by:

```
install.packages("devtools")
```

After devtools is installed, install and load the package by:

```
devtools::install_github("jetroant/rbsvar")
library(rbsvar)
```
